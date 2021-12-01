#!/bin/bash

target="${1/Doc/Doc\/html}"

cp "$1" "$target"

jupyter nbconvert --to markdown "$target"

echo [TOC]
echo
cat "${target/.ipynb/.md}" | sed 's/```python/```{.py}/g' | sed 's/^##/#/g' | tail -n +2
