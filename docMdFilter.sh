#!/bin/bash

cat $1 | sed 's/\\(\(.*\)\\)/\\f$\1\\f$/g;s/\\\[/\\f[/g;s/\\\]/\\f]/g;s/$$\(.\+\)$\$/\\f[\1\\f]/g;s/$\([^$]\+\)\$/\\f$\1\\f$/g'
