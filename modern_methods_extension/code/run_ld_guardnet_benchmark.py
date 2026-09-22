#!/usr/bin/env python
"""Synthetic, leakage-aware benchmark for an LD-graph input anomaly detector.

The model is intentionally positioned as a triage layer for content-level faults
that deterministic metadata contracts cannot authenticate. It never overrides a
contract failure and it is not a biological association or colocalisation model.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import platform
import random
import time
from dataclasses import asdict, dataclass
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
import numpy as np
import pandas as pd
import sklearn
import torch
from sklearn.calibration import calibration_curve
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import (
    accuracy_score,
    brier_score_loss,
    confusion_matrix,
    f1_score,
    log_loss,
    roc_auc_score,
)
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from torch import nn
from torch.utils.data import DataLoader, Dataset


CLASSES = ("valid", "order_mismatch", "sign_corruption", "ld_mismatch")
ARCHITECTURES = ("one_shared", "two_shared", "shared_plus_specific", "distinct_only")
PANELS = ("EUR", "EAS", "AFR", "SAS", "AMR")


@dataclass(frozen=True)
class Config:
    seed: int = 20260920
    train_replicates: int = 220
    validation_replicates: int = 60
    test_replicates: int = 100
    max_epochs: int = 100
    patience: int = 15
    batch_size: int = 128
    learning_rate: float = 1e-3
    weight_decay: float = 1e-4
    hidden_dim: int = 48
    dropout: float = 0.20
    training_seeds: tuple[int, ...] = (11, 29, 47)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--ld-dir", required=True, type=Path)
    parser.add_argument("--out-dir", required=True, type=Path)
    parser.add_argument("--quick", action="store_true", help="Run a small smoke benchmark.")
    return parser.parse_args()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def load_ld_matrices(ld_dir: Path) -> tuple[list[str], dict[str, np.ndarray], dict[str, str]]:
    variants: list[str] | None = None
    matrices: dict[str, np.ndarray] = {}
    hashes: dict[str, str] = {}
    for panel in PANELS:
        path = ld_dir / f"common_ld_{panel}.tsv"
        frame = pd.read_csv(path, sep="\t")
        current = frame.iloc[:, 0].astype(str).tolist()
        if variants is None:
            variants = current
        elif current != variants:
            raise ValueError(f"variant order differs for {panel}")
        matrix = frame.iloc[:, 1:].to_numpy(dtype=np.float64)
        matrix = nearest_correlation(matrix)
        matrices[panel] = matrix.astype(np.float32)
        hashes[panel] = sha256_file(path)
    if variants is None:
        raise ValueError("no LD matrices loaded")
    return variants, matrices, hashes


def nearest_correlation(matrix: np.ndarray, floor: float = 1e-8) -> np.ndarray:
    symmetric = (matrix + matrix.T) / 2
    values, vectors = np.linalg.eigh(symmetric)
    rebuilt = (vectors * np.maximum(values, floor)) @ vectors.T
    scale = np.sqrt(np.diag(rebuilt))
    rebuilt = rebuilt / np.outer(scale, scale)
    rebuilt = (rebuilt + rebuilt.T) / 2
    np.fill_diagonal(rebuilt, 1.0)
    return rebuilt


def draw_correlated(rng: np.random.Generator, matrix: np.ndarray) -> np.ndarray:
    values, vectors = np.linalg.eigh(matrix.astype(np.float64))
    return vectors @ (np.sqrt(np.maximum(values, 0)) * rng.standard_normal(len(values)))


def causal_configuration(replicate: int, architecture: str, p: int) -> tuple[list[int], list[int]]:
    first = (replicate * 37 - 1) % p
    second = (first + p // 3) % p
    third = (first + (2 * p) // 3) % p
    if architecture == "one_shared":
        return [first], [first]
    if architecture == "two_shared":
        return [first, second], [first, second]
    if architecture == "shared_plus_specific":
        return [first, second], [first, third]
    return [first], [second]


def simulate_z(
    rng: np.random.Generator,
    matrix: np.ndarray,
    sample_size: int,
    causal: list[int],
    total_effect: float,
) -> np.ndarray:
    effects = np.repeat(total_effect / math.sqrt(len(causal)), len(causal))
    mean = math.sqrt(sample_size) * (matrix[:, causal] @ effects)
    return mean + draw_correlated(rng, matrix)


def standardize(vector: np.ndarray) -> np.ndarray:
    scale = vector.std(ddof=0)
    if not np.isfinite(scale) or scale <= 0:
        scale = 1.0
    return (vector - vector.mean()) / scale


def node_features(z1: np.ndarray, z2: np.ndarray) -> np.ndarray:
    first = standardize(z1)
    second = standardize(z2)
    return np.column_stack(
        [first, second, np.abs(first), np.abs(second), np.sign(first * second)]
    ).astype(np.float32)


def apply_fault(
    rng: np.random.Generator,
    label: str,
    z1: np.ndarray,
    z2: np.ndarray,
    data_panel: str,
    available_fit_panels: tuple[str, ...],
) -> tuple[np.ndarray, np.ndarray, str]:
    altered = z2.copy()
    fit_panel = data_panel
    if label == "order_mismatch":
        step = 53
        permutation = ((np.arange(len(altered)) + 1) * step - 1) % len(altered)
        altered = altered[permutation]
    elif label == "sign_corruption":
        offset = int(rng.integers(0, 10))
        altered[offset::10] *= -1
    elif label == "ld_mismatch":
        alternatives = [panel for panel in available_fit_panels if panel != data_panel]
        fit_panel = alternatives[int(rng.integers(0, len(alternatives)))]
    return z1, altered, fit_panel


def build_split(
    *,
    matrices: dict[str, np.ndarray],
    replicate_start: int,
    replicates: int,
    data_panels: tuple[str, ...],
    fit_panels: tuple[str, ...],
    architectures: tuple[str, ...],
    sample_sizes: tuple[int, ...],
    seed: int,
) -> dict[str, np.ndarray | list[dict[str, object]]]:
    xs: list[np.ndarray] = []
    fit_panel_indices: list[int] = []
    ys: list[int] = []
    handcrafted: list[np.ndarray] = []
    metadata: list[dict[str, object]] = []
    p = next(iter(matrices.values())).shape[0]
    panel_to_index = {panel: index for index, panel in enumerate(PANELS)}

    for local_rep in range(replicates):
        replicate = replicate_start + local_rep
        for data_panel in data_panels:
            for architecture in architectures:
                for sample_size in sample_sizes:
                    rng = np.random.default_rng(
                        seed + replicate * 10007 + panel_to_index[data_panel] * 1009
                        + ARCHITECTURES.index(architecture) * 101 + sample_size
                    )
                    causal1, causal2 = causal_configuration(replicate, architecture, p)
                    data_ld = matrices[data_panel]
                    z1 = simulate_z(rng, data_ld, sample_size, causal1, 0.085)
                    z2 = simulate_z(rng, data_ld, sample_size, causal2, 0.075)
                    for class_index, label in enumerate(CLASSES):
                        local_rng = np.random.default_rng(rng.integers(0, 2**32 - 1))
                        fz1, fz2, fit_panel = apply_fault(
                            local_rng, label, z1, z2, data_panel, fit_panels
                        )
                        features = node_features(fz1, fz2)
                        xs.append(features)
                        fit_panel_indices.append(panel_to_index[fit_panel])
                        ys.append(class_index)
                        handcrafted.append(
                            aggregate_baseline_features(features, matrices[fit_panel])
                        )
                        metadata.append(
                            {
                                "replicate": replicate,
                                "data_panel": data_panel,
                                "fit_panel": fit_panel,
                                "architecture": architecture,
                                "sample_size": sample_size,
                                "label": label,
                            }
                        )
    return {
        "x": np.stack(xs),
        "panel": np.asarray(fit_panel_indices, dtype=np.int64),
        "y": np.asarray(ys, dtype=np.int64),
        "baseline": np.stack(handcrafted),
        "metadata": metadata,
    }


def aggregate_baseline_features(x: np.ndarray, matrix: np.ndarray) -> np.ndarray:
    signed = matrix.copy().astype(np.float64)
    np.fill_diagonal(signed, 0)
    denominator = np.sum(np.abs(signed), axis=1, keepdims=True)
    denominator[denominator == 0] = 1
    operator = signed / denominator
    z1, z2 = x[:, 0], x[:, 1]
    smooth1, smooth2 = operator @ z1, operator @ z2
    return np.asarray(
        [
            np.mean(np.abs(z1 - smooth1)),
            np.mean(np.abs(z2 - smooth2)),
            np.sqrt(np.mean((z1 - smooth1) ** 2)),
            np.sqrt(np.mean((z2 - smooth2) ** 2)),
            np.corrcoef(z1, z2)[0, 1],
            np.mean(np.sign(z1 * z2)),
            np.quantile(np.abs(z1), 0.95),
            np.quantile(np.abs(z2), 0.95),
            np.mean(np.abs(matrix[np.triu_indices_from(matrix, 1)])),
        ],
        dtype=np.float32,
    )


def signed_operators(matrices: dict[str, np.ndarray]) -> tuple[torch.Tensor, torch.Tensor]:
    positive, negative = [], []
    for panel in PANELS:
        matrix = matrices[panel].astype(np.float64).copy()
        np.fill_diagonal(matrix, 0)
        pos = np.maximum(matrix, 0)
        neg = np.maximum(-matrix, 0)
        for operator in (pos, neg):
            operator += np.eye(len(operator))
            degree = operator.sum(axis=1)
            inv = 1 / np.sqrt(np.maximum(degree, 1e-12))
            operator[:] = inv[:, None] * operator * inv[None, :]
        positive.append(pos.astype(np.float32))
        negative.append(neg.astype(np.float32))
    return torch.from_numpy(np.stack(positive)), torch.from_numpy(np.stack(negative))


class GraphDataset(Dataset):
    def __init__(self, split: dict[str, object]):
        self.x = torch.from_numpy(split["x"])
        self.panel = torch.from_numpy(split["panel"])
        self.y = torch.from_numpy(split["y"])

    def __len__(self) -> int:
        return len(self.y)

    def __getitem__(self, index: int) -> tuple[torch.Tensor, torch.Tensor, torch.Tensor]:
        return self.x[index], self.panel[index], self.y[index]


class SignedLDGuardNet(nn.Module):
    def __init__(self, input_dim: int, hidden_dim: int, dropout: float):
        super().__init__()
        self.layer1 = nn.Sequential(
            nn.Linear(input_dim * 3, hidden_dim), nn.LayerNorm(hidden_dim), nn.GELU(),
            nn.Dropout(dropout),
        )
        self.layer2 = nn.Sequential(
            nn.Linear(hidden_dim * 3, hidden_dim), nn.LayerNorm(hidden_dim), nn.GELU(),
            nn.Dropout(dropout),
        )
        self.classifier = nn.Sequential(
            nn.Linear(hidden_dim * 2, hidden_dim), nn.GELU(), nn.Dropout(dropout),
            nn.Linear(hidden_dim, len(CLASSES)),
        )

    def forward(
        self, x: torch.Tensor, panel: torch.Tensor,
        positive: torch.Tensor, negative: torch.Tensor,
    ) -> torch.Tensor:
        apos, aneg = positive[panel], negative[panel]
        h = self.layer1(torch.cat([x, torch.bmm(apos, x), torch.bmm(aneg, x)], dim=-1))
        h = self.layer2(torch.cat([h, torch.bmm(apos, h), torch.bmm(aneg, h)], dim=-1))
        pooled = torch.cat([h.mean(dim=1), h.amax(dim=1)], dim=-1)
        return self.classifier(pooled)


def set_all_seeds(seed: int) -> None:
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    torch.use_deterministic_algorithms(True)


def probability_metrics(y: np.ndarray, probability: np.ndarray) -> dict[str, float]:
    prediction = probability.argmax(axis=1)
    anomaly_true = (y != 0).astype(int)
    anomaly_probability = 1 - probability[:, 0]
    output = {
        "accuracy": accuracy_score(y, prediction),
        "macro_f1": f1_score(y, prediction, average="macro"),
        "macro_ovr_auc": roc_auc_score(y, probability, multi_class="ovr", average="macro"),
        "log_loss": log_loss(y, probability, labels=np.arange(len(CLASSES))),
        "anomaly_auc": roc_auc_score(anomaly_true, anomaly_probability),
        "anomaly_brier": brier_score_loss(anomaly_true, anomaly_probability),
    }
    for index, label in enumerate(CLASSES):
        output[f"auc_{label}"] = roc_auc_score((y == index).astype(int), probability[:, index])
        output[f"recall_{label}"] = np.mean(prediction[y == index] == index)
    return output


@torch.no_grad()
def predict_network(
    model: nn.Module,
    split: dict[str, object],
    positive: torch.Tensor,
    negative: torch.Tensor,
    batch_size: int,
) -> np.ndarray:
    model.eval()
    loader = DataLoader(GraphDataset(split), batch_size=batch_size, shuffle=False)
    probabilities = []
    for x, panel, _ in loader:
        probabilities.append(torch.softmax(model(x, panel, positive, negative), dim=1).cpu().numpy())
    return np.vstack(probabilities)


def train_network(
    train: dict[str, object],
    validation: dict[str, object],
    positive: torch.Tensor,
    negative: torch.Tensor,
    config: Config,
    training_seed: int,
) -> tuple[nn.Module, dict[str, object]]:
    set_all_seeds(config.seed + training_seed)
    model = SignedLDGuardNet(train["x"].shape[-1], config.hidden_dim, config.dropout)
    optimizer = torch.optim.AdamW(
        model.parameters(), lr=config.learning_rate, weight_decay=config.weight_decay
    )
    loss_function = nn.CrossEntropyLoss()
    generator = torch.Generator().manual_seed(config.seed + training_seed)
    loader = DataLoader(
        GraphDataset(train), batch_size=config.batch_size, shuffle=True, generator=generator
    )
    best_state, best_auc, best_epoch = None, -np.inf, -1
    history = []
    for epoch in range(config.max_epochs):
        model.train()
        losses = []
        for x, panel, y in loader:
            optimizer.zero_grad(set_to_none=True)
            loss = loss_function(model(x, panel, positive, negative), y)
            loss.backward()
            optimizer.step()
            losses.append(float(loss.detach()))
        validation_probability = predict_network(
            model, validation, positive, negative, config.batch_size
        )
        validation_auc = roc_auc_score(
            validation["y"], validation_probability, multi_class="ovr", average="macro"
        )
        history.append({"epoch": epoch + 1, "loss": np.mean(losses), "validation_auc": validation_auc})
        if validation_auc > best_auc + 1e-5:
            best_auc, best_epoch = validation_auc, epoch + 1
            best_state = {key: value.detach().clone() for key, value in model.state_dict().items()}
        elif epoch + 1 - best_epoch >= config.patience:
            break
    if best_state is None:
        raise RuntimeError("training produced no checkpoint")
    model.load_state_dict(best_state)
    return model, {"best_epoch": best_epoch, "best_validation_auc": best_auc, "history": history}


def bootstrap_interval(
    y: np.ndarray, probability: np.ndarray, metric: str, seed: int,
    repeats: int = 2000, groups: np.ndarray | None = None,
) -> tuple[float, float]:
    if groups is None or len(groups) != len(y):
        raise ValueError("explicit outer-replicate groups are required")
    members = [np.flatnonzero(groups == group) for group in np.unique(groups)]
    rng = np.random.default_rng(seed)
    values = []
    for _ in range(repeats):
        indices = np.concatenate([members[i] for i in rng.integers(0, len(members), len(members))])
        if len(np.unique(y[indices])) < len(CLASSES):
            continue
        current = probability_metrics(y[indices], probability[indices])[metric]
        values.append(current)
    return tuple(np.quantile(values, [0.025, 0.975]))


def save_figures(
    out_dir: Path,
    y: np.ndarray,
    probability: np.ndarray,
    baseline_probability: np.ndarray,
) -> None:
    prediction = probability.argmax(axis=1)
    matrix = confusion_matrix(y, prediction, normalize="true")
    anomaly_true = (y != 0).astype(int)
    anomaly_probability = 1 - probability[:, 0]
    baseline_anomaly_probability = 1 - baseline_probability[:, 0]
    fraction, mean_prediction = calibration_curve(
        anomaly_true, anomaly_probability, n_bins=10, strategy="quantile"
    )
    baseline_fraction, baseline_mean = calibration_curve(
        anomaly_true, baseline_anomaly_probability, n_bins=10, strategy="quantile"
    )

    figure, axes = plt.subplots(1, 3, figsize=(10.6, 3.35), constrained_layout=True)
    axes[0].set_axis_off()
    axes[0].set_title("a", loc="left", fontweight="bold")
    stages = [
        (0.82, "Node features", "z1, z2, |z1|, |z2|, sign"),
        (0.61, "Signed LD propagation", "positive / negative channels"),
        (0.40, "Two message blocks", "normalization, GELU, dropout"),
        (0.19, "Graph readout", "mean + max pooling; 4 states"),
    ]
    colours = ("#D9EAF7", "#E8DFF5", "#DDEEDB", "#FCE5CD")
    for index, ((ypos, title, detail), colour) in enumerate(zip(stages, colours)):
        box = FancyBboxPatch(
            (0.08, ypos - 0.07), 0.84, 0.14,
            boxstyle="round,pad=0.012,rounding_size=0.02",
            linewidth=0.8, edgecolor="#3A3A3A", facecolor=colour,
            transform=axes[0].transAxes,
        )
        axes[0].add_patch(box)
        axes[0].text(0.50, ypos + 0.018, title, ha="center", va="center",
                     fontsize=8.5, fontweight="bold", transform=axes[0].transAxes)
        axes[0].text(0.50, ypos - 0.033, detail, ha="center", va="center",
                     fontsize=7.2, transform=axes[0].transAxes)
        if index < len(stages) - 1:
            axes[0].annotate(
                "", xy=(0.50, ypos - 0.105), xytext=(0.50, ypos - 0.075),
                xycoords=axes[0].transAxes, textcoords=axes[0].transAxes,
                arrowprops={"arrowstyle": "-|>", "lw": 0.8, "color": "#555555"},
            )

    image = axes[1].imshow(matrix, vmin=0, vmax=1, cmap="Blues")
    axes[1].set_xticks(range(len(CLASSES)), CLASSES, rotation=35, ha="right", fontsize=8)
    axes[1].set_yticks(range(len(CLASSES)), CLASSES, fontsize=8)
    axes[1].set_xlabel("Predicted state")
    axes[1].set_ylabel("True state")
    axes[1].set_title("b", loc="left", fontweight="bold")
    for row in range(len(CLASSES)):
        for column in range(len(CLASSES)):
            axes[1].text(column, row, f"{matrix[row, column]:.2f}", ha="center", va="center",
                         color="white" if matrix[row, column] > 0.55 else "black", fontsize=8)
    figure.colorbar(image, ax=axes[1], fraction=0.046, pad=0.04)

    axes[2].plot([0, 1], [0, 1], "--", color="#777777", linewidth=1)
    axes[2].plot(mean_prediction, fraction, "o-", color="#1F77B4", label="LD-GuardNet")
    axes[2].plot(baseline_mean, baseline_fraction, "s-", color="#D95F02", label="Logistic baseline")
    axes[2].set(xlabel="Predicted anomaly probability", ylabel="Observed anomaly proportion",
                xlim=(0, 1), ylim=(0, 1))
    axes[2].set_title("c", loc="left", fontweight="bold")
    axes[2].legend(frameon=False, fontsize=8)
    for spine in ("top", "right"):
        axes[2].spines[spine].set_visible(False)

    figure.savefig(out_dir / "figure_ld_guardnet_ood.png", dpi=600, facecolor="white")
    figure.savefig(out_dir / "figure_ld_guardnet_ood.svg", facecolor="white")
    plt.close(figure)


def main() -> None:
    args = parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)
    config = Config()
    if args.quick:
        config = Config(
            train_replicates=24, validation_replicates=12, test_replicates=16,
            max_epochs=20, patience=6, training_seeds=(11,), batch_size=64,
        )
    started = time.time()
    variants, matrices, input_hashes = load_ld_matrices(args.ld_dir)
    positive, negative = signed_operators(matrices)

    train = build_split(
        matrices=matrices, replicate_start=1, replicates=config.train_replicates,
        data_panels=("EUR", "EAS", "AFR"), fit_panels=("EUR", "EAS", "AFR"),
        architectures=("one_shared", "two_shared"), sample_sizes=(5000,), seed=config.seed,
    )
    validation = build_split(
        matrices=matrices, replicate_start=10001, replicates=config.validation_replicates,
        data_panels=("EUR", "EAS", "AFR"), fit_panels=("EUR", "EAS", "AFR"),
        architectures=("one_shared", "two_shared"), sample_sizes=(5000,), seed=config.seed,
    )
    test = build_split(
        matrices=matrices, replicate_start=20001, replicates=config.test_replicates,
        data_panels=("SAS", "AMR"), fit_panels=("SAS", "AMR"),
        architectures=("shared_plus_specific", "distinct_only"),
        sample_sizes=(2000, 20000), seed=config.seed,
    )

    baseline = make_pipeline(
        StandardScaler(), LogisticRegression(max_iter=3000, random_state=config.seed)
    )
    baseline.fit(train["baseline"], train["y"])
    baseline_probability = baseline.predict_proba(test["baseline"])
    baseline_metrics = probability_metrics(test["y"], baseline_probability)

    seed_probabilities, run_rows, histories = [], [], {}
    for training_seed in config.training_seeds:
        model, training = train_network(
            train, validation, positive, negative, config, training_seed
        )
        validation_probability = predict_network(
            model, validation, positive, negative, config.batch_size
        )
        test_probability = predict_network(model, test, positive, negative, config.batch_size)
        seed_probabilities.append(test_probability)
        run_metrics = probability_metrics(test["y"], test_probability)
        run_rows.append({"training_seed": training_seed, **run_metrics, **{k: training[k] for k in ("best_epoch", "best_validation_auc")}})
        histories[str(training_seed)] = training["history"]
        torch.save(
            {"state_dict": model.state_dict(), "config": asdict(config), "classes": CLASSES},
            args.out_dir / f"ld_guardnet_seed_{training_seed}.pt",
        )

    ensemble_probability = np.mean(seed_probabilities, axis=0)
    ensemble_metrics = probability_metrics(test["y"], ensemble_probability)
    success = (
        ensemble_metrics["macro_ovr_auc"] >= 0.80
        and min(ensemble_metrics[f"auc_{label}"] for label in CLASSES) >= 0.70
        and ensemble_metrics["macro_f1"] >= 0.60
        and ensemble_metrics["macro_ovr_auc"] - baseline_metrics["macro_ovr_auc"] >= 0.03
    )
    ci = {
        metric: bootstrap_interval(test["y"], ensemble_probability, metric, config.seed + index,
                                   groups=np.asarray([row['replicate'] for row in test['metadata']]))
        for index, metric in enumerate(("accuracy", "macro_f1", "macro_ovr_auc", "anomaly_auc"))
    }

    pd.DataFrame(run_rows).to_csv(args.out_dir / "training_seed_metrics.tsv", sep="\t", index=False)
    summary_rows = []
    for model_name, metrics in (("LD-GuardNet ensemble", ensemble_metrics), ("Logistic baseline", baseline_metrics)):
        summary_rows.append({"model": model_name, **metrics})
    pd.DataFrame(summary_rows).to_csv(args.out_dir / "model_summary.tsv", sep="\t", index=False)

    prediction_frame = pd.DataFrame(test["metadata"])
    prediction_frame["true_class"] = test["y"]
    for index, label in enumerate(CLASSES):
        prediction_frame[f"prob_{label}"] = ensemble_probability[:, index]
        prediction_frame[f"baseline_prob_{label}"] = baseline_probability[:, index]
    prediction_frame.to_csv(args.out_dir / "ood_predictions.tsv", sep="\t", index=False)
    save_figures(args.out_dir, test["y"], ensemble_probability, baseline_probability)

    report = {
        "schema_version": "1.0",
        "scope": "synthetic LD-graph anomaly triage; not biological inference",
        "model_name": "LD-GuardNet",
        "architecture": (
            "two signed LD message-passing blocks; node mean/max pooling; four-state classifier"
        ),
        "classes": CLASSES,
        "predeclared_inclusion_rule": (
            "Formal-run OOD macro one-vs-rest AUROC >=0.80, every class AUROC >=0.70, "
            "macro F1 >=0.60, and macro AUROC improvement over the logistic baseline >=0.03; "
            "locked after code smoke testing and before the formal run"
        ),
        "inclusion_rule_passed": bool(success),
        "config": asdict(config),
        "split_contract": {
            "train": "EUR/EAS/AFR; N=5000; one_shared and two_shared architectures",
            "validation": "new seeds under the training regimes",
            "ood_test": "unseen SAS/AMR panels; N=2000/20000; shared_plus_specific and distinct_only architectures",
        },
        "n_variants": len(variants),
        "n_train": len(train["y"]),
        "n_validation": len(validation["y"]),
        "n_ood_test": len(test["y"]),
        "ensemble_metrics": ensemble_metrics,
        "bootstrap_95_ci": {key: list(value) for key, value in ci.items()},
        "baseline_metrics": baseline_metrics,
        "input_sha256": input_hashes,
        "runtime_seconds": time.time() - started,
        "environment": {
            "python": platform.python_version(),
            "platform": platform.platform(),
            "numpy": np.__version__,
            "pandas": pd.__version__,
            "sklearn": sklearn.__version__,
            "torch": torch.__version__,
            "cuda_available": torch.cuda.is_available(),
        },
        "limitations": [
            "All examples are synthetic and use one 100-variant genomic region.",
            "The OOD test changes ancestry panel, architecture, sample size and seed, but not genomic region.",
            "The detector is probabilistic triage and cannot replace deterministic contract checks.",
            "No association, colocalisation, causal, biological or clinical claim is supported.",
        ],
    }
    (args.out_dir / "benchmark_report.json").write_text(
        json.dumps(report, indent=2), encoding="utf-8"
    )
    (args.out_dir / "training_history.json").write_text(
        json.dumps(histories, indent=2), encoding="utf-8"
    )
    print(json.dumps({"passed": success, "ensemble": ensemble_metrics, "baseline": baseline_metrics}, indent=2))


if __name__ == "__main__":
    main()
