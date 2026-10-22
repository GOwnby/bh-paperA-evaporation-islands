#!/usr/bin/env python3
"""verify.py - run every generator of this repository and compare the computed numbers with the numbers quoted
in main.tex.   Usage: python3 verify.py [--quick|--full]     (normally called through `bash verify.sh`)

  --quick  (default)  all generators (symbolic, mpmath, scipy WKB, time-domain at 3 grids, direct shooting, figure)
  --full              additionally: finer time-domain grid (Nx=24000) and halved RK4 steps in the direct shooting

Outputs (all under verification/):  logs/<script>.log (stdout of each script), logs/<script>.stderr (only if non-empty),
checks.json, REPORT.md (PASS/FAIL table + environment), timings.txt.  The scripts also refresh ancillary/output/.
Exit code: 0 all checks passed; 1 at least one check FAILED; 3 no failure but some checks SKIPPED because a Python
module is missing (the summary line then says INCOMPLETE).  Expected values: verification/expected.json.
"""
import argparse, datetime, importlib.util, json, math, os, platform, re, shutil, subprocess, sys, time

ROOT = os.path.dirname(os.path.abspath(__file__))
ANC = os.path.join(ROOT, "ancillary")
OUT = os.path.join(ANC, "output")
VER = os.path.join(ROOT, "verification")
LOGS = os.path.join(VER, "logs")
REPO = os.path.basename(ROOT)
T_START = 0.0

# name, script, args, required modules, rough duration, full-only, json/log it produces
JOBS = [
    ("check_symbolic",        "check_symbolic.py",        [],              ["sympy"],             "~30 s",   False),
    ("check_analytic",        "check_analytic.py",        [],              ["mpmath", "numpy"],   "~3 s",    False),
    ("wkb",                   "wkb.py",                   [],              ["scipy"],             "~5 s",    False),
    ("crosscheck_wkb_mpmath", "crosscheck_wkb_mpmath.py", [],              ["scipy", "mpmath"],   "~10 s",   False),
    ("timedomain_nx6000",     "timedomain.py",            [],              ["scipy"],             "~30 s",   False),
    ("timedomain_nx3000",     "timedomain.py",            ["--nx", "3000"], ["scipy"],            "~10 s",   False),
    ("timedomain_nx12000",    "timedomain.py",            ["--nx", "12000"], ["scipy"],           "~2.5 min", False),
    ("make_fig_sublimation",  "make_fig_sublimation.py",  [],              ["matplotlib", "numpy"], "~3 s",  False),
    ("compare_wkb_original",  "compare_wkb_original.py",  [],              ["scipy"],             "~40 s",   False),
    ("timedomain_nx24000",    "timedomain.py",            ["--nx", "24000"], ["scipy"],           "~10 min", True),
    ("fit_timedomain",        "fit_timedomain.py",        None,            ["scipy"],             "~3 s",    False),   # args set below
    ("qnm_shooting",          "qnm_shooting.py",          [],              ["scipy"],             "~1.5 min", False),
    ("qnm_shooting_refine2",  "qnm_shooting.py",          ["2"],           ["scipy"],             "~5 min",  True),
]


PRODUCER = {"symbolic_checks.json": "check_symbolic", "analytic_checks.json": "check_analytic", "wkb.json": "wkb",
            "crosscheck_wkb_mpmath.json": "crosscheck_wkb_mpmath", "timedomain_fit.json": "fit_timedomain",
            "qnm_shooting.json": "qnm_shooting", "sublimation_curve.json": "make_fig_sublimation",
            "compare_wkb_original.json": "compare_wkb_original"}


def have(mod):
    return importlib.util.find_spec(mod) is not None


def versions():
    v = {"python": sys.version.split()[0], "platform": platform.platform()}
    for m in ["numpy", "scipy", "sympy", "mpmath", "matplotlib"]:
        try:
            v[m] = __import__(m).__version__
        except Exception:
            v[m] = "NOT INSTALLED"
    try:
        v["gcc"] = subprocess.run(["gcc", "--version"], capture_output=True, text=True).stdout.splitlines()[0] + " (not used by this repository)"
    except Exception:
        v["gcc"] = "not found (not used by this repository)"
    return v


def run_job(name, script, args, full):
    cmd = [sys.executable, script] + args
    t0 = time.time()
    p = subprocess.run(cmd, cwd=ANC, capture_output=True, text=True)
    dt = time.time() - t0
    with open(os.path.join(LOGS, name + ".log"), "w") as fh:
        fh.write(p.stdout)
    errp = os.path.join(LOGS, name + ".stderr")
    if p.stderr.strip():
        with open(errp, "w") as fh:
            fh.write(p.stderr)
    elif os.path.exists(errp):
        open(errp, "w").close()                    # keep stale content out; (empty file = no stderr)
    return p.returncode, dt


def load_json(fn):
    p = os.path.join(OUT, fn)
    if not os.path.exists(p) or os.path.getmtime(p) < T_START:
        return None
    with open(p) as fh:
        return json.load(fh)


def fetch(cache, src):
    fn, *path = src
    d = cache.get(fn)
    if d is None:
        raise KeyError(f"{fn} not available")
    for k in path:
        d = d[k]
    return d


def table_rows_in_tex():
    rows = []
    with open(os.path.join(ROOT, "main.tex")) as fh:
        for line in fh:
            if re.match(r"^\$0\.\d\d\$ & \$\d+\$ & \$", line.strip()):
                rows.append(re.sub(r"\s+", " ", line.strip()))
    return rows


def main():
    ap = argparse.ArgumentParser()
    g = ap.add_mutually_exclusive_group()
    g.add_argument("--quick", action="store_true")
    g.add_argument("--full", action="store_true")
    a = ap.parse_args()
    full = a.full
    mode = "full" if full else "quick"
    os.makedirs(LOGS, exist_ok=True)
    # stale results must not leak into this run: old time-domain signals are removed if possible, and a JSON result counts only
    # if it was written during this run (modification time check in load_json; no file deletion needed)
    global T_START
    T_START = time.time()
    for fn in os.listdir(OUT) if os.path.isdir(OUT) else []:
        if fn.startswith("td_l") and fn.endswith(".npz"):
            try:
                os.remove(os.path.join(OUT, fn))
            except OSError:
                pass

    jobstat = {}
    timings = []
    nxs = ["3000", "6000", "12000"] + (["24000"] if full else [])
    for name, script, args, mods, est, fullonly in JOBS:
        if fullonly and not full:
            continue
        if args is None:
            args = nxs
        missing = [m for m in mods if not have(m)]
        if missing:
            jobstat[name] = ("SKIPPED", f"module(s) missing: {', '.join(missing)}")
            print(f"[{name}] SKIPPED (missing: {', '.join(missing)})", flush=True)
            continue
        print(f"[{name}] running python3 {script} {' '.join(args)}  (expected {est}) ...", flush=True)
        rc, dt = run_job(name, script, args, full)
        timings.append((name, dt))
        jobstat[name] = ("OK" if rc == 0 else "FAILED", f"exit code {rc}")
        print(f"[{name}] {'ok' if rc == 0 else 'FAILED (exit code %d, see verification/logs/%s.log/.stderr)' % (rc, name)}  {dt:.1f} s", flush=True)

    cache = {fn: load_json(fn) for fn in ["symbolic_checks.json", "analytic_checks.json", "wkb.json", "crosscheck_wkb_mpmath.json", "timedomain_fit.json",
                                          "qnm_shooting.json", "qnm_shooting_refine2.json", "sublimation_curve.json", "compare_wkb_original.json"]}
    # derived quantities
    D = {}
    try:
        x = cache["crosscheck_wkb_mpmath.json"]["rows"]
        D["wkb_mp_max_rel"] = max(max(r["rel_wI"], r["rel_Q"], r["rel_tau"], r["rel_barrier"]) for r in x)
        D["wkb_bs_max_resid"] = max(abs(r["BS_residual"]) for r in x)
    except Exception:
        pass
    try:
        t = cache["timedomain_fit.json"]; runs = t["runs"]
        D["td_wr_reldiff"] = abs(runs["6000"]["wR"] - t["wkb"]["wR"])/t["wkb"]["wR"]
        D["td_decay_grid_change"] = max(abs(runs[k]["decay"]/runs["6000"]["decay"] - 1) for k in ("3000", "12000"))
        D["td_wr_grid_change"] = max(abs(runs[k]["wR"]/runs["6000"]["wR"] - 1) for k in ("3000", "12000"))
        D["td_decay_lsq_vs_env"] = max(abs(runs[k]["decay_lsq"]/runs[k]["decay"] - 1) for k in ("3000", "6000", "12000"))
        if "24000" in runs:
            D["td_decay_change_24000"] = abs(runs["24000"]["decay"]/runs["12000"]["decay"] - 1)
    except Exception:
        pass
    try:
        q = cache["qnm_shooting.json"]; t = cache["timedomain_fit.json"]
        D["qnm_td_wr_diff"] = abs(q["direct"]["0.90_2"]["re"] - t["runs"]["6000"]["wR"])
        D["qnm_td_decay_reldiff"] = abs(abs(q["direct"]["0.90_2"]["im"]) - t["runs"]["6000"]["decay"])/t["runs"]["6000"]["decay"]
        rat = [q["direct"][k]["ratio_to_wkb"] for k in ("0.90_1", "0.90_2", "0.90_4", "0.95_1", "0.95_2", "0.99_1")]
        D["qnm_ratio_min"], D["qnm_ratio_max"] = min(rat), max(rat)
        D["qnm_rm_robust"] = max(max(abs(q["robustness_rm1.8"][k]["re"] - q["direct"][k]["re"]), abs(q["robustness_rm1.8"][k]["im"] - q["direct"][k]["im"]))
                                 for k in q["robustness_rm1.8"])
        if cache.get("qnm_shooting_refine2.json"):
            r2 = cache["qnm_shooting_refine2.json"]
            D["qnm_refine_diff"] = max(abs(complex(r2["direct"][k]["re"], r2["direct"][k]["im"]) - complex(q["direct"][k]["re"], q["direct"][k]["im"])) for k in q["direct"])
    except Exception:
        pass
    cache["derived"] = D
    tex_rows = table_rows_in_tex()

    expected = json.load(open(os.path.join(VER, "expected.json")))
    results = []
    for e in expected:
        if e.get("mode") == "full" and not full:
            continue
        r = dict(id=e["id"], where=e["where"], quantity=e["quantity"], reason=e.get("reason", ""), cmp=e["cmp"])
        try:
            v = fetch(cache, e["src"])
        except Exception as ex:
            # a missing result is a SKIP only if the producing script was skipped (missing module); otherwise it is a FAILURE
            producer = PRODUCER.get(e["src"][0])
            skipped = (jobstat.get(producer, ("",))[0] == "SKIPPED") if producer else any(s == "SKIPPED" for s, _ in jobstat.values())
            r.update(status="SKIPPED" if skipped else "FAIL", computed=None, paper=e.get("paper"), tol="",
                     note=f"result not available ({type(ex).__name__}: {ex})")
            results.append(r); continue
        if e.get("absval"):
            v = abs(v)
        if "paper_ref" in e:
            try:
                e = dict(e, paper=fetch(cache, e["paper_ref"]))
            except Exception:
                r.update(status="FAIL", computed=v, paper=None, tol="", note="reference value not available"); results.append(r); continue
        c = e["cmp"]; tol = e.get("tol"); paper = e.get("paper")
        if c == "abs":
            ok = abs(v - paper) <= tol; tols = f"abs {tol:g}"
        elif c == "rel":
            ok = abs(v - paper) <= tol*abs(paper); tols = f"rel {tol:g}"
        elif c == "log10":
            ok = abs(math.log10(v) - math.log10(paper)) <= tol; tols = f"|dlog10| {tol:g}"
        elif c == "range":
            lo, hi = e.get("lo"), e.get("hi")
            ok = (lo is None or v >= lo) and (hi is None or v <= hi); tols = f"in [{lo}, {hi}]"
        elif c == "true":
            ok = (v is True); tols = "exact"
        elif c == "table_row":
            norm = re.sub(r"\s+", " ", v.strip())
            ok = norm in tex_rows; tols = "exact string"; paper = "row present in main.tex"
        else:
            raise ValueError(c)
        r.update(status="PASS" if ok else "FAIL", computed=v, paper=paper, tol=tols, note="")
        results.append(r)

    # job failures that produced no checks still count
    failed_jobs = [n for n, (s, _) in jobstat.items() if s == "FAILED"]
    n_pass = sum(r["status"] == "PASS" for r in results)
    n_fail = sum(r["status"] == "FAIL" for r in results) + len(failed_jobs)
    n_skip = sum(r["status"] == "SKIPPED" for r in results)
    n_tot = len(results)

    def fmt(x):
        if x is None: return ""
        if isinstance(x, bool): return str(x)
        if isinstance(x, float): return f"{x:.6g}"
        return str(x)

    v = versions()
    L = []
    L.append(f"# Verification report - {REPO}\n")
    L.append(f"Mode: **{mode}**.  Date (UTC): {datetime.datetime.now(datetime.timezone.utc).strftime('%Y-%m-%d')}.\n")
    L.append("## Summary\n")
    L.append(f"* checks passed: **{n_pass}/{n_tot}**, failed: {n_fail}, skipped: {n_skip}")
    L.append("* scripts: " + ", ".join(f"{n}: {s}" for n, (s, _) in jobstat.items()) + "\n")
    L.append("## Environment\n")
    for k in ["python", "numpy", "scipy", "sympy", "mpmath", "matplotlib", "gcc", "platform"]:
        L.append(f"* {k}: {v[k]}")
    L.append("")
    L.append("## Checks (paper value vs value computed by the scripts in `ancillary/`)\n")
    L.append("| ID | Paper location | Quantity | Paper value | Computed | Tolerance (reason) | Status |")
    L.append("|---|---|---|---|---|---|---|")
    for r in results:
        tol = r["tol"] + (f" ({r['reason']})" if r["reason"] else "")
        L.append(f"| {r['id']} | {r['where']} | {r['quantity']} | {fmt(r['paper'])} | {fmt(r['computed'])} | {tol} | {r['status']}{(' - ' + r['note']) if r['note'] else ''} |")
    L.append("")
    L.append("Symbolic (S01-S19) rows are exact sympy derivations, see `ancillary/output/symbolic_checks.txt` for the derived expressions.")
    open(os.path.join(VER, "REPORT.md"), "w").write("\n".join(L) + "\n")
    json.dump(dict(mode=mode, passed=n_pass, failed=n_fail, skipped=n_skip, total=n_tot, jobs={n: s for n, (s, _) in jobstat.items()}, checks=results),
              open(os.path.join(VER, "checks.json"), "w"), indent=1, default=float)
    with open(os.path.join(VER, "timings.txt"), "w") as fh:
        fh.write(f"wall-clock seconds per script ({mode}), this machine; sum = {sum(d for _, d in timings):.0f} s\n")
        for n, d in timings:
            fh.write(f"{n:28s} {d:8.1f}\n")

    for r in results:
        if r["status"] == "FAIL":
            print(f"FAIL {r['id']}: {r['quantity']} paper={fmt(r['paper'])} computed={fmt(r['computed'])} tol={r['tol']}")
    for n in failed_jobs:
        print(f"FAILED script: {n}")
    if n_fail:
        print(f"RESULT {REPO}: FAILED {n_fail} of {n_tot} checks ({mode})")
        return 1
    if n_skip:
        print(f"RESULT {REPO}: INCOMPLETE - {n_pass}/{n_tot} checks passed, {n_skip} SKIPPED (missing modules: "
              + ", ".join(sorted({m for n, (s, why) in jobstat.items() if s == 'SKIPPED' for m in why.split(': ')[1].split(', ')})) + f") ({mode})")
        return 3
    print(f"RESULT {REPO}: PASSED {n_pass}/{n_tot} checks ({mode})")
    return 0


if __name__ == "__main__":
    sys.exit(main())
