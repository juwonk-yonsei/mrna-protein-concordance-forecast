import os

# Stage 1 caps CPU at 16. Stage 2 sets B_CPU_CAP=none before import.
_cap = os.environ.get("B_CPU_CAP", "16")
if _cap != "none":
    os.environ["OMP_NUM_THREADS"] = _cap
    os.environ["OPENBLAS_NUM_THREADS"] = _cap
    os.environ["MKL_NUM_THREADS"] = _cap
    os.environ["NUMEXPR_MAX_THREADS"] = _cap

import hashlib
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parents[1]


def load_paths():
    with open(ROOT / "config" / "paths.yaml") as f:
        cfg = yaml.safe_load(f)
    out = {"root": ROOT}
    for key, value in cfg.items():
        if key.endswith("_dir"):
            out[key] = ROOT / value
        else:
            out[key] = value
    return out


def load_plan():
    with open(ROOT / "config" / "analysis_plan_B1.yaml") as f:
        return yaml.safe_load(f)


def sha256_file(path, chunk=1 << 20):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        while True:
            b = f.read(chunk)
            if not b:
                break
            h.update(b)
    return h.hexdigest()


def free_gb(path=None):
    st = os.statvfs(path or ROOT)
    return st.f_bavail * st.f_frsize / 1e9
