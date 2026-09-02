"""Record numerical gates for a stored LD matrix without altering it."""
from __future__ import annotations
import argparse, json
from pathlib import Path
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument('--ld-npz', required=True)
parser.add_argument('--out', required=True)
args = parser.parse_args()
x = np.load(args.ld_npz, allow_pickle=False)['ld']
e = np.linalg.eigvalsh(x)
ae = np.abs(e)
largest = float(np.max(ae))
smallest = float(np.min(ae))
tol = float(x.shape[0] * np.finfo(float).eps * largest)
rank = int(np.count_nonzero(ae > tol))
rank_deficient = rank < x.shape[0]
result = {
    'n_variants': int(x.shape[0]),
    'symmetric_max_abs_error': float(np.max(np.abs(x-x.T))),
    'diagonal_max_abs_error_from_one': float(np.max(np.abs(np.diag(x)-1))),
    'finite': bool(np.isfinite(x).all()),
    'min_eigenvalue': float(e[0]),
    'max_eigenvalue': float(e[-1]),
    'smallest_absolute_eigenvalue': smallest,
    'numerical_rank': rank,
    'rank_tolerance': tol,
    'rank_deficient': rank_deficient,
    'condition_number_2': None if rank_deficient else largest / max(smallest, np.finfo(float).tiny),
    'condition_number_is_infinite': rank_deficient,
}
Path(args.out).write_text(json.dumps(result, indent=2)+'\n', encoding='utf-8')
print(json.dumps(result, sort_keys=True))
