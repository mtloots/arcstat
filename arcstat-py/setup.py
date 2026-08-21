"""Build the shared arc-length C back-end at INSTALL time.

The package previously shelled out to `cc` on FIRST IMPORT and wrote the resulting library into
its own site-packages directory. That works on the author's machine and fails everywhere else:
it needs a compiler at import rather than install, so the failure surfaces as a crash instead of
an install error; it writes into site-packages, which is read-only in containers and many system
installs; and it hardcoded `cc`, which stock Windows does not have.

Declaring the sources as an Extension hands the compile to pip, inside build isolation, at the
point where a missing compiler is a clear install-time error -- and lets real platform wheels be
built, after which users need no compiler at all.

init.c is excluded deliberately: it is R's symbol-registration glue, it references R headers, and
it is not part of the Python front.
"""
import glob
import os

from setuptools import Extension, setup

# Paths must be RELATIVE to setup.py: setuptools rejects absolute source paths outright
# ("setup script specifies an absolute path"), and a relative path is what ends up in the sdist.
SRC = os.path.join("src", "arcstat")

sources = sorted(
    s for s in glob.glob(os.path.join(SRC, "*.c"))
    if os.path.basename(s) != "init.c"
)

# Named with a leading underscore so it is never mistaken for an importable module: it carries no
# PyInit, and is loaded with ctypes.CDLL by _core.py.
ext = Extension(
    name="arcstat._libarcstat",
    sources=sources,
    include_dirs=[SRC],
    libraries=[] if os.name == "nt" else ["m"],
    optional=False,
)

setup(ext_modules=[ext])
