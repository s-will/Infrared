#!/bin/bash

target="${1/Doc/Doc\/html}"

cp "$1" "$target"

jupyter nbconvert --to markdown "$target"

cat "${target/.ipynb/.md}" | sed 's/```python/```{.py}/g'
