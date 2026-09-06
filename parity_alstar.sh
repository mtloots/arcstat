#!/bin/sh
# Parity harness for the spacing-corrected band arc length and the scale estimators.
#
# The R package arcstat and the Python package arcstat-py compile the SAME C source, arcdistc.c, so
# every value below must agree to the last printed digit. The comparison is byte-identical: the two
# sides print with the same format and the files are diffed, not compared numerically. The data are
# supplied as a fixed literal vector rather than drawn, so no generator has to be matched across the
# two front ends.
#
# Usage: sh parity_alstar.sh
set -e
OUT=${TMPDIR:-/tmp}/parity_alstar
mkdir -p "$OUT"
PKG="$(cd "$(dirname "$0")" && pwd)"
sh "$PKG/parity_sync_guard.sh" "$PKG/arcstat/src" "$PKG/arcstat-py" || exit 1

# fixed data, written once and read by both sides
Rscript -e 'set.seed(20260906); cat(sprintf("%.17g", rnorm(4000)), sep="\n")' > "$OUT/x.txt"

Rscript -e '
suppressMessages(devtools::load_all("'"$PKG"'/arcstat", quiet=TRUE))
x <- scan("'"$OUT"'/x.txt", quiet=TRUE)
for (b in list(c(0.05,0.95), c(0.10,0.90), c(0.25,0.75)))
  for (s in c(0.25, 1, 5))
    cat(sprintf("star %.2f %.2f %.2f %.12f\n", b[1], b[2], s, al_band_model_star(s, b[1], b[2], 4000)))
for (b in list(c(0.05,0.95), c(0.10,0.90), c(0.25,0.75))) {
  cat(sprintf("scale %.2f %.2f %.12f\n", b[1], b[2], al_scale(x, b[1], b[2])))
  cat(sprintf("raw   %.2f %.2f %.12f\n", b[1], b[2], al_scale_raw(x, b[1], b[2])))
}' > "$OUT/r.txt"

python3 -c "
import sys; sys.path.insert(0, '$PKG/arcstat-py/src')
import arcstat as A
x = [float(l) for l in open('$OUT/x.txt')]
for b in [(0.05,0.95),(0.10,0.90),(0.25,0.75)]:
    for s in [0.25,1.0,5.0]:
        print('star %.2f %.2f %.2f %.12f' % (b[0], b[1], s, A.al_band_model_star(s,b[0],b[1],4000)))
for b in [(0.05,0.95),(0.10,0.90),(0.25,0.75)]:
    print('scale %.2f %.2f %.12f' % (b[0], b[1], A.al_scale(x,b[0],b[1])))
    print('raw   %.2f %.2f %.12f' % (b[0], b[1], A.al_scale_raw(x,b[0],b[1])))
" > "$OUT/py.txt"

if diff -q "$OUT/r.txt" "$OUT/py.txt" >/dev/null; then
  echo "VERDICT: PARITY, $(wc -l < "$OUT/r.txt" | tr -d ' ') values byte-identical"
else
  echo "VERDICT: MISMATCH"; diff "$OUT/r.txt" "$OUT/py.txt"; exit 1
fi
