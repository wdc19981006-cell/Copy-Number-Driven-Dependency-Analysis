"""Backward-compatible entry point for the generalized repository index guard."""
import json
from git_size_guard import inspect_index

report = inspect_index()
print(json.dumps(report, indent=2))
if report['issues']:
    raise SystemExit('Staged audit failed; do not commit/push')
