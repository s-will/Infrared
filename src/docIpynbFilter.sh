#!/bin/bash


target="${1/Doc\//}"
dir="$(dirname "$target")"
base="$(basename "$target")"

target="Doc/html/$dir_$base"


cp "$1" "$target"

jupyter nbconvert --to markdown "$target"


head -n 1 "${target/.ipynb/.md}"

echo [TOC]
echo

cat "${target/.ipynb/.md}" \
  | sed 's/```python/```{.py}/g' \
  | sed 's/@/@ /g' \
  | sed 's/^##/#/g' \
  | sed 's/\\(\(.*\)\\)/\\f$\1\\f$/g;s/\\\[/\\f[/g;s/\\\]/\\f]/g;s/$$\(.\+\)$\$/\\f[\1\\f]/g;s/$\([^$]\+\)\$/\\f$\1\\f$/g' \
  | tail -n +2
