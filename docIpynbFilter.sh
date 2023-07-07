#!/bin/bash

target="${1/Doc\//}"
dir="$(dirname $1)"
base="$(basename $1)"

target="Doc/html/$dir_$base"

cp "$1" "$target"

jupyter nbconvert --to markdown "$target"

echo [TOC]
echo

cat "${target/.ipynb/.md}" \
  | sed 's/```python/```{.py}/g' \
  | sed 's/@/@ /g' \
  | sed 's/^##/#/g' \
  | sed 's/\\(\(.*\)\\)/\\f$\1\\f$/g;s/\\\[/\\f[/g;s/\\\]/\\f]/g;s/$$\(.\+\)$\$/\\f[\1\\f]/g;s/$\([^$]\+\)\$/\\f$\1\\f$/g' \
  | tail -n +2
