"""Shared output helpers (no numerical code).  Everything a script prints is also written to
ancillary/output/<name>; machine-readable results go to ancillary/output/<name>.json and are
read by ../verify.py.  SciPy is required by the scripts that need it (there is no fallback)."""
import atexit, json, os, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "output")
os.makedirs(OUT, exist_ok=True)


class _Tee:
    def __init__(self, path):
        self.f = open(path, "w")
        self.o = sys.stdout

    def write(self, s):
        self.o.write(s)
        if not self.f.closed:
            self.f.write(s)

    def flush(self):
        self.o.flush()
        if not self.f.closed:
            self.f.flush()


def tee_stdout(name, with_scipy=True):
    """Redirect stdout also to output/<name>; first line records the library versions."""
    t = _Tee(os.path.join(OUT, name))
    atexit.register(t.f.close)
    sys.stdout = t
    v = f"numpy {np.__version__}"
    if with_scipy:
        import scipy
        v += f"; scipy {scipy.__version__}"
    print(f"# {v}")


def save_json(name, d):
    with open(os.path.join(OUT, name), "w") as fh:
        json.dump(d, fh, indent=1, sort_keys=True, default=float)
        fh.write("\n")
