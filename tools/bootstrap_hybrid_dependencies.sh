#!/usr/bin/env bash
# Rootless development-header/runtime cache for Ubuntu 22.04.
set -euo pipefail
task_root=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)
task_cache="$task_root/.cache/dev-packages"
command -v apt-get >/dev/null
command -v dpkg-deb >/dev/null
for task_family in gmsh cgns; do
  if [[ "$task_family" = gmsh ]]; then
    task_header=gmsh.h; task_development=libgmsh-dev; task_runtime=libgmsh4.8
  else
    task_header=cgnslib.h; task_development=libcgns-dev; task_runtime=libcgns3.4
  fi
  if [[ -f "/usr/include/$task_header" ]]; then
    printf '%s\n' "$task_family: system development header available"
    continue
  fi
  if [[ -f "$task_cache/$task_family/usr/include/$task_header" ]]; then
    printf '%s\n' "$task_family: existing local development cache available"
    continue
  fi
  mkdir -p "$task_cache/$task_family"
  (
    cd "$task_cache/$task_family"
    apt-get download "$task_development" "$task_runtime"
    for task_archive in ./*.deb; do dpkg-deb -x "$task_archive" .; done
  )
done
printf '%s\n' 'Gmsh/CGNS cache ready. OCCT, Qt 6, Eigen and JSON remain system prerequisites.'
