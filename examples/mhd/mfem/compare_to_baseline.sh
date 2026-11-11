#!/usr/bin/env python3
"""
Compare a new lmmhd run against a stored baseline.

    python3 scripts/compare_baseline.py BASELINE_DIR [NEW_DIR] [options]

NEW_DIR defaults to the current directory.  In each directory the script
looks for simulation.log and console.log (whichever exist), and a ParaView
output folder (any sub-directory containing Cycle*/ folders, e.g. lmmhd/ or
data_*/lmmhd/).

Logs are compared line by line: each line is split into its text (with
numbers replaced by '#') and its numbers.  The text must match exactly; the
numbers must agree to within --rtol / --atol.  Wall-clock values ("time
elapsed", "Total simulation time") are ignored.

If meshio is installed (pip install meshio), the fields in the last common
ParaView cycle are compared too.  Both runs must use the same number of MPI
ranks for that to be meaningful.

Exit status: 0 = match, 1 = differences found, 2 = usage / missing files.
"""

import argparse
import math
import re
import sys
from pathlib import Path

NUM = re.compile(r"[-+]?(?:\d+\.\d*|\.\d+|\d+)(?:[eE][-+]?\d+)?|[-+]?(?:nan|inf)\b",
                 re.IGNORECASE)

# Wall-clock numbers that legitimately change between runs.
IGNORE = [
    re.compile(r"time elapsed\s*=\s*\S+"),
    re.compile(r"Total simulation time:\s*\S+"),
]

# Blank lines are not compared.
SKIP_LINE = re.compile(r"^\s*$")

# ScalingTimer::report() prints a wall-clock table at the end of the run
# ("=== Timings: N ranks, M DOFs ===" followed by one row per section).
# Everything from that header onwards is timing, so stop reading there.
TIMINGS_START = re.compile(r"^=== Timings:")

LOG_NAMES = ["simulation.log", "console.log"]


def split_line(line):
    for pat in IGNORE:
        line = pat.sub("<ignored>", line)
    nums = [float(x) for x in NUM.findall(line)]
    text = NUM.sub("#", line).strip()
    return text, nums


def load_log(path, ignore=()):
    out = []
    for lineno, raw in enumerate(path.read_text(errors="replace").splitlines(), 1):
        if TIMINGS_START.search(raw):
            break
        if SKIP_LINE.search(raw):
            continue
        if any(pat.search(raw) for pat in ignore):
            continue
        text, nums = split_line(raw)
        out.append((lineno, raw.rstrip(), text, nums))
    return out


def rel_diff(a, b, atol):
    if math.isnan(a) or math.isnan(b):
        return 0.0 if (math.isnan(a) and math.isnan(b)) else math.inf
    d = abs(a - b)
    if d <= atol:
        return 0.0
    return d / max(abs(a), abs(b))


def compare_log(base_path, new_path, rtol, atol, max_report, ignore=()):
    base, new = load_log(base_path, ignore), load_log(new_path, ignore)
    print(f"\n== {base_path.name}: {len(base)} vs {len(new)} comparable lines")

    n_bad = 0
    worst = (0.0, None)
    for (bl, braw, btext, bnums), (nl, nraw, ntext, nnums) in zip(base, new):
        if btext != ntext or len(bnums) != len(nnums):
            print(f"  STRUCTURE DIVERGES (baseline line {bl}, new line {nl}) -- stopping this file")
            print(f"    baseline: {braw}")
            print(f"    new     : {nraw}")
            print("    (usually a different FGMRES iteration count or an extra/missing print)")
            return False, worst
        for a, b in zip(bnums, nnums):
            r = rel_diff(a, b, atol)
            if r > worst[0]:
                worst = (r, (bl, braw, nraw))
            if r > rtol:
                n_bad += 1
                if n_bad <= max_report:
                    print(f"  line {bl}: {a!r} vs {b!r}  (rel diff {r:.3e})")
                    print(f"    baseline: {braw}")
                    print(f"    new     : {nraw}")
                break

    ok = n_bad == 0
    if len(base) != len(new):
        ok = False
        longer = "baseline" if len(base) > len(new) else "new run"
        print(f"  LENGTH DIFFERS: {longer} has {abs(len(base) - len(new))} extra lines "
              "(run stopped early/late?)")

    if n_bad > max_report:
        print(f"  ... {n_bad - max_report} more mismatching lines")
    print(f"  max relative difference: {worst[0]:.3e}")
    if worst[1] and worst[0] > 0:
        print(f"    at baseline line {worst[1][0]}: {worst[1][1].strip()}")
    print("  OK" if ok else "  DIFFERENT")
    return ok, worst


def find_paraview_root(d, exclude=None):
    exclude = exclude.resolve() if exclude else None
    hits = set()
    for p in d.rglob("Cycle*"):
        if not p.is_dir():
            continue
        if exclude and exclude in p.resolve().parents:
            continue  # don't pick up the baseline when NEW_DIR contains it
        hits.add(p.parent)
    return sorted(hits)


def compare_paraview(base_dir, new_dir, rtol, atol):
    try:
        import meshio  # noqa: F401
        import numpy as np
    except ImportError:
        print("\n== ParaView output: skipped (pip install meshio numpy to enable)")
        return True

    broots, nroots = find_paraview_root(base_dir), find_paraview_root(new_dir, exclude=base_dir)
    if not broots or not nroots:
        print("\n== ParaView output: not found in one of the directories, skipped")
        return True
    if len(nroots) > 1:
        print(f"\n== ParaView output: several candidates in {new_dir}:")
        for r in nroots:
            print(f"     {r}")
        print("   using the newest; pass the run directory explicitly to be sure")
        nroots.sort(key=lambda p: p.stat().st_mtime)
    broot, nroot = broots[-1], nroots[-1]

    bcyc = {p.name for p in broot.glob("Cycle*")}
    ncyc = {p.name for p in nroot.glob("Cycle*")}
    common = sorted(bcyc & ncyc)
    if bcyc != ncyc:
        print(f"\n== ParaView: cycle sets differ ({len(bcyc)} vs {len(ncyc)})")
    if not common:
        print("== ParaView: no common cycles")
        return False
    cyc = common[-1]

    bfiles = sorted((broot / cyc).glob("*.vtu"))
    nfiles = sorted((nroot / cyc).glob("*.vtu"))
    print(f"\n== ParaView {broot.name}/{cyc}: {len(bfiles)} vs {len(nfiles)} .vtu files")
    if [f.name for f in bfiles] != [f.name for f in nfiles]:
        print("  different proc files -- were the runs on different numbers of ranks?")
        return False

    ok = True
    for bf, nf in zip(bfiles, nfiles):
        bm, nm = meshio.read(bf), meshio.read(nf)
        if bm.points.shape != nm.points.shape or not np.allclose(bm.points, nm.points,
                                                                 rtol=rtol, atol=atol):
            print(f"  {bf.name}: mesh points differ")
            ok = False
            continue
        for name in sorted(set(bm.point_data) | set(nm.point_data)):
            if name not in bm.point_data or name not in nm.point_data:
                print(f"  {bf.name}: field '{name}' only in one run")
                ok = False
                continue
            a = np.asarray(bm.point_data[name], dtype=float)
            b = np.asarray(nm.point_data[name], dtype=float)
            scale = max(np.abs(a).max(), np.abs(b).max(), atol)
            err = np.abs(a - b).max() / scale
            flag = "" if err <= rtol else "   <-- DIFFERENT"
            if err > rtol:
                ok = False
            print(f"  {bf.name:20s} {name:22s} max|diff|/max|field| = {err:.3e}{flag}")
    print("  OK" if ok else "  DIFFERENT")
    return ok


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("baseline", type=Path)
    ap.add_argument("new", type=Path, nargs="?", default=Path("."))
    ap.add_argument("--rtol", type=float, default=1e-8,
                    help="relative tolerance (default 1e-8; use 0 to demand bitwise-identical text values)")
    ap.add_argument("--atol", type=float, default=1e-12,
                    help="differences below this are ignored (default 1e-12, for ~0 quantities like ||div J||)")
    ap.add_argument("--max-report", type=int, default=5)
    ap.add_argument("--no-vtu", action="store_true", help="skip ParaView comparison")
    ap.add_argument("--ignore", action="append", default=[], metavar="REGEX",
                    help="skip log lines matching REGEX (repeatable), e.g. "
                         "--ignore 'FGMRES iteration' --ignore 'Pass :'")
    ap.add_argument("--logs", nargs="+", default=LOG_NAMES, metavar="NAME",
                    help="which log files to compare (default: simulation.log console.log). "
                         "Use --logs simulation.log for MOOSE runs, whose console output differs.")
    ap.add_argument("--results-only", action="store_true",
                    help="shorthand: ignore Krylov iteration histories, compare only "
                         "per-step results (useful across different rank counts)")
    args = ap.parse_args()

    if not args.baseline.is_dir() or not args.new.is_dir():
        print("baseline and new must be directories", file=sys.stderr)
        return 2

    # Inputs first: comparing runs with different parameters is meaningless.
    bp, np_ = args.baseline / "params.in", args.new / "params.in"
    if bp.exists() and np_.exists() and bp.read_text() != np_.read_text():
        print("WARNING: params.in differs between baseline and new run "
              "(command-line overrides are not checked)")

    ignore = [re.compile(r) for r in args.ignore]
    if args.results_only:
        ignore += [re.compile(r"FGMRES iteration\s*:"), re.compile(r"^\s*Pass\s*:")]

    all_ok = True
    found = False
    for name in args.logs:
        b, n = args.baseline / name, args.new / name
        if b.exists() and n.exists():
            found = True
            ok, _ = compare_log(b, n, args.rtol, args.atol, args.max_report, ignore)
            all_ok &= ok
        elif b.exists() or n.exists():
            print(f"\n== {name}: present in only one directory, skipped")
    if not found:
        print("no logs to compare", file=sys.stderr)
        return 2

    if not args.no_vtu:
        all_ok &= compare_paraview(args.baseline, args.new, args.rtol, args.atol)

    print("\nRESULT:", "MATCH" if all_ok else "DIFFERENT")
    return 0 if all_ok else 1


if __name__ == "__main__":
    sys.exit(main())
