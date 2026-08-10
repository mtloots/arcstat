#!/bin/sh
# Shared pre-condition for every parity harness.
#
# Usage:  sh parity_sync_guard.sh <r-package-dir> <python-package-dir>
#
# Each Python front carries its OWN copy of the shared C back-end and links its own library.
# A parity run made against a stale copy therefore compares the new C on one side with the old
# C on the other, and it will still report agreement on every value the edit happened not to
# touch. That is the worst possible failure mode for this harness: it is silent, and it is most
# silent exactly when the edit is narrow. So no parity verdict is worth reading until the
# computational sources are known to agree.
#
# init.c is excluded by design. It is R's own symbol-registration glue, it is not compiled into
# the Python library, and it carries no numerics; the two copies legitimately differ.
set -e
RDIR="$1"; PYDIR="$2"
if [ ! -d "$RDIR" ] || [ ! -d "$PYDIR" ]; then
  echo "VERDICT: GUARD ERROR -- cannot find $RDIR or $PYDIR"; exit 1
fi
drift=""
for rc in "$RDIR"/*.c "$RDIR"/*.h; do
  [ -e "$rc" ] || continue
  b=$(basename "$rc")
  [ "$b" = "init.c" ] && continue
  pc=$(find "$PYDIR" -name "$b" -not -path "*Rcheck*" 2>/dev/null | head -1)
  [ -n "$pc" ] || continue
  cmp -s "$rc" "$pc" || drift="$drift $b"
done
if [ -n "$drift" ]; then
  echo "VERDICT: OUT OF SYNC --$drift"
  echo "         the Python front's copy differs from the R package source."
  echo "         Copy across and rebuild the Python library before trusting any parity result."
  exit 1
fi
