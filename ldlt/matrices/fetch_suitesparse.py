#!/usr/bin/env python3
"""Download SuiteSparse Matrix Collection .mtx fixtures.

Reads matrices/manifest.txt (one `group/name` per line). Fetches each as a
tar.gz from https://suitesparse-collection-website.herokuapp.com/MM/<group>/<name>.tar.gz,
extracts the .mtx into matrices/<name>.mtx. Idempotent.
"""
import os, sys, tarfile, urllib.request, pathlib

ROOT = pathlib.Path(__file__).parent
MANIFEST = ROOT / "manifest.txt"
BASE = "https://suitesparse-collection-website.herokuapp.com/MM"

def fetch(group: str, name: str) -> None:
    out = ROOT / f"{name}.mtx"
    if out.exists():
        print(f"[skip] {out.name}")
        return
    url = f"{BASE}/{group}/{name}.tar.gz"
    print(f"[get ] {url}")
    tmp = ROOT / f"{name}.tar.gz"
    urllib.request.urlretrieve(url, tmp)
    with tarfile.open(tmp) as tf:
        for m in tf.getmembers():
            if m.name.endswith(f"/{name}.mtx"):
                m.name = f"{name}.mtx"
                tf.extract(m, ROOT)
                break
    tmp.unlink()
    print(f"[ok  ] {out}")

def main() -> int:
    if not MANIFEST.exists():
        print("manifest.txt missing", file=sys.stderr); return 1
    for line in MANIFEST.read_text().splitlines():
        s = line.strip()
        if not s or s.startswith("#"): continue
        if "/" not in s: continue
        group, name = s.split("/", 1)
        try:
            fetch(group, name)
        except Exception as e:
            print(f"[fail] {s}: {e}", file=sys.stderr)
    return 0

if __name__ == "__main__":
    sys.exit(main())
