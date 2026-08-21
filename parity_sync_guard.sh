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
# The same argument applies to the COMPILED library. The Python front compiles once and caches
# libarcstat, so sources that match can still be fronted by a library built before the last edit;
# if a signature changed, ctypes will push the wrong number of arguments and segfault, and if only
# a numeric constant changed it will quietly report the old value. Any library older than the
# sources it was built from is therefore removed here, and the next import rebuilds it.
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
# `read -d ""` is a bash extension, and this script declares #!/bin/sh. On macOS sh is bash in
# posix mode so it worked locally, but on Ubuntu sh is dash and the harness died with
# "read: Illegal option -d" -- which CI found on 21 Aug 2026. find -exec is POSIX and still safe
# for the paths with spaces the original -print0 was chosen to handle, because find passes each
# path as a single argument with no word splitting.
# The library glob is lib*, not libarcstat*: this guard is shared, and parity_ellstat.sh points
# it at the ellstat front, whose library is libellstat.
find "$PYDIR" \( -name "lib*.so" -o -name "lib*.dylib" \) -not -path "*Rcheck*" \
     -exec sh -c '
       lib="$1"; RDIR="$2"
       for rc in "$RDIR"/*.c "$RDIR"/*.h; do
         [ -e "$rc" ] || continue
         if [ "$rc" -nt "$lib" ]; then
           echo "NOTE: removing $(basename "$lib"), older than $(basename "$rc"); import will rebuild it"
           rm -f "$lib"; break
         fi
       done
     ' sh {} "$RDIR" \; 2>/dev/null

if [ -n "$drift" ]; then
  echo "VERDICT: OUT OF SYNC --$drift"
  echo "         the Python front's copy differs from the R package source."
  echo "         Copy across and rebuild the Python library before trusting any parity result."
  exit 1
fi
