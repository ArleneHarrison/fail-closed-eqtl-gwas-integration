"""Write 600-dpi PNG metadata and refresh sibling figure-manifest hashes."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

from PIL import Image


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("manifests", nargs="+", type=Path)
    args = parser.parse_args()
    for manifest_path in args.manifests:
        payload = json.loads(manifest_path.read_text(encoding="utf-8"))
        png_path = Path(payload["outputs"]["png"]["path"])
        temporary = png_path.with_suffix(".dpi-stamped.png")
        with Image.open(png_path) as image:
            image.save(temporary, format="PNG", dpi=(600, 600), optimize=False)
        temporary.replace(png_path)
        payload["outputs"]["png"]["sha256"] = sha256(png_path)
        payload["outputs"]["png"]["dpi"] = 600
        manifest_path.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
        with Image.open(png_path) as image:
            actual = image.info.get("dpi")
        if actual is None or min(actual) < 599:
            raise RuntimeError(f"DPI metadata was not written: {png_path}")
        print(f"{png_path}: dpi={actual}, sha256={sha256(png_path)}")


if __name__ == "__main__":
    main()
