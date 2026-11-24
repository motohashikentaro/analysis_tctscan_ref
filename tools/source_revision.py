"""Identify the implementation saved with analysis conditions."""

import hashlib
from pathlib import Path
import sys

root = Path(__file__).resolve().parent.parent
paths = [root / "Makefile", Path(__file__).resolve()]
for directory in ("include", "src", "app"):
    paths.extend(path for path in (root / directory).iterdir() if path.is_file())

digest = hashlib.sha256()
for path in sorted(paths):
    digest.update(str(path.relative_to(root)).encode() + b"\0")
    digest.update(path.read_bytes() + b"\0")

content = '#define TCTSCAN_SOURCE_SHA256 "' + digest.hexdigest() + '"\n'
target = Path(sys.argv[1])
target.parent.mkdir(parents=True, exist_ok=True)
if not target.exists() or target.read_text() != content:
    target.write_text(content)
