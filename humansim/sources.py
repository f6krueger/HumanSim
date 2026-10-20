"""Citation sidecar. The model runs without it; it is loaded for reports and integrity checks."""

from __future__ import annotations

from pathlib import Path

import yaml

from .registry import DATA_DIR


def load_sources(data_dir: str | Path | None = None) -> dict[str, dict]:
    root = Path(data_dir) if data_dir else DATA_DIR
    path = root / "sources" / "sources.yaml"
    if not path.exists():
        return {}
    doc = yaml.safe_load(path.read_text(encoding="utf-8")) or {}
    return dict(doc.get("sources") or {})
