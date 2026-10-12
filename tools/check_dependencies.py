#!/usr/bin/env python3
"""Check the installed Debian package versions against the development lock."""

import json
import pathlib
import subprocess
import sys


def main() -> int:
    lock_path = pathlib.Path(sys.argv[1])
    lock = json.loads(lock_path.read_text())
    failures = []
    for package, expected in lock["packages"].items():
        result = subprocess.run(
            ["dpkg-query", "-W", "-f=${Version}", package],
            text=True,
            capture_output=True,
            check=False,
        )
        actual = result.stdout.strip() if result.returncode == 0 else "not-installed"
        if actual != expected["version"]:
            failures.append(f"{package}: expected {expected['version']}, found {actual}")
    if failures:
        print("dependency lock mismatch:", file=sys.stderr)
        for failure in failures:
            print(f"  {failure}", file=sys.stderr)
        return 1
    print(f"dependency lock: {lock['platform']} packages match")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
