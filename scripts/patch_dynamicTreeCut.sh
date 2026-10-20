#!/usr/bin/env bash
# Vendors dynamicTreeCut 0.1.1 (Python port of the R Dynamic Tree Cut) into scripts/analysis/ and applies the small
# fixes needed for NumPy 2 / recent SciPy (np.in1d -> np.isin; integer casting of rankdata output used as an index).
# The PyPI wheel build fails on recent Python, so the source distribution is used directly.
set -euo pipefail
HERE=$(cd "$(dirname "$0")" && pwd)
TMP=$(mktemp -d)
pip download dynamicTreeCut==0.1.1 --no-deps --no-binary :all: -d "$TMP" >/dev/null
tar -xzf "$TMP"/dynamicTreeCut-0.1.1.tar.gz -C "$TMP"
rm -rf "$HERE/analysis/dynamicTreeCut"
cp -r "$TMP/dynamicTreeCut-0.1.1/dynamicTreeCut" "$HERE/analysis/dynamicTreeCut"
patch -d "$HERE/analysis/dynamicTreeCut" -p1 < "$HERE/dynamicTreeCut_numpy2.patch"
echo "dynamicTreeCut vendored and patched in $HERE/analysis/dynamicTreeCut"
