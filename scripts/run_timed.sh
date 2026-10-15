#!/bin/bash
# Usage: scripts/run_timed.sh STEP_NAME command...
set -euo pipefail
# Stage 1 default is 16. Stage 2 passes B_CPU_CAP=none and does not export a cap.
if [ "${B_CPU_CAP:-16}" != "none" ]; then
  export OMP_NUM_THREADS="${B_CPU_CAP:-16}"
  export OPENBLAS_NUM_THREADS="${B_CPU_CAP:-16}"
  export MKL_NUM_THREADS="${B_CPU_CAP:-16}"
  export NUMEXPR_MAX_THREADS="${B_CPU_CAP:-16}"
fi
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
STEP="$1"
shift
mkdir -p "$ROOT/logs"
TMP="$(mktemp)"
/usr/bin/time -v -o "$TMP" "$@"
"$ROOT/.venv/bin/python" - "$ROOT" "$STEP" "$TMP" "$*" << 'PY'
import sys
from datetime import datetime, timezone
from pathlib import Path
root, step, tmp, cmd = sys.argv[1:]
text = Path(tmp).read_text(errors="replace")
elapsed = ""
rss = ""
for line in text.splitlines():
    if "Elapsed (wall clock) time" in line:
        elapsed = line.split("):", 1)[-1].strip()
    if "Maximum resident set size" in line:
        rss = line.split(":", 1)[-1].strip()
out = Path(root) / "logs" / "resource_usage.tsv"
if not out.exists():
    out.write_text("step\telapsed_wall\tmax_rss_kb\tdatetime\tcommand\n")
stamp = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
with out.open("a") as f:
    f.write(f"{step}\t{elapsed}\t{rss}\t{stamp}\t{cmd}\n")
Path(tmp).unlink(missing_ok=True)
PY
