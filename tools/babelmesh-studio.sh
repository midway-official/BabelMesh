#!/usr/bin/env bash
set -euo pipefail
task_root=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)
task_binary="$task_root/build/qt/bin/babelmesh-gui"
if [[ ! -x "$task_binary" ]]; then
  printf '%s\n' 'Build the Qt application first: make gui' >&2
  exit 1
fi
exec "$task_binary" "$@"
