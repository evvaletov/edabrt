#!/usr/bin/env python3
"""Regression guard for the edabrt aberration formulas.

Property checks that would have caught the n1>3 x-plane bug:
  1. First-order symplecticity: (x|x)(a|a) - (x|a)(a|x) == 1 in every field-index regime.
  2. Continuity of the second-order coefficients across the n1=3 branch boundary.
  3. Golden values for the n1=4 hyperbolic case.

Compiles edabrt.c and runs the binary; needs only python3 and a C compiler. Exits 1 on failure.
"""
import os, re, subprocess, sys, tempfile

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def build():
    out = os.path.join(tempfile.mkdtemp(), "edabrt")
    subprocess.run(["cc", "-O2", "-o", out, os.path.join(ROOT, "edabrt.c"), "-lm"], check=True)
    return out


def run(binp, r, ang, n1, n2):
    txt = subprocess.run([binp, str(r), str(ang), str(n1), str(n2)],
                         capture_output=True, text=True).stdout
    res, plane = {}, None
    for ln in txt.splitlines():
        s = ln.strip()
        if s.startswith("(x|"):
            plane = "x"
        elif s.startswith("(a|"):
            plane = "a"
        f = ln.split()
        if plane and len(f) == 5 and re.match(r"-?\d", f[1]) and f[2] in ("1", "2"):
            res[(plane, f[3], f[4])] = float(f[1])
    return res


def c(d, plane, e1, e2):
    return d.get((plane, str(e1), str(e2)), 0.0)  # edabrt omits coefficients that are exactly 0


def main():
    binp = build()
    fails = []

    # 1. Symplecticity of the first-order map across all regimes (n1<3, n1=3, n1>3).
    for n1 in (0.5, 1.0, 2.0, 2.5, 3.0, 3.5, 4.0, 6.0):
        d = run(binp, 1.0, 60.0, n1, 0.0)
        det = c(d, "x", 1, 0) * c(d, "a", 0, 1) - c(d, "x", 0, 1) * c(d, "a", 1, 0)
        if abs(det - 1.0) > 1e-9:
            fails.append(f"symplecticity n1={n1}: det={det:.12f} != 1")

    # 2. Continuity across n1=3: both sides must match the n1=3 branch value. The bug
    #    jumped by ~26 here, while the smooth variation over dn1=0.001 is ~0.01.
    for n2 in (0.0, 2.0, -0.5):
        mid = run(binp, 1.0, 60.0, 3.0, n2)
        for tag, side in (("2.999", run(binp, 1.0, 60.0, 2.999, n2)),
                          ("3.001", run(binp, 1.0, 60.0, 3.001, n2))):
            for plane in ("x", "a"):
                for e1, e2 in ((2, 0), (1, 1), (0, 2)):
                    m, v = c(mid, plane, e1, e2), c(side, plane, e1, e2)
                    if abs(v - m) > 0.05:
                        fails.append(f"discontinuity ({plane}|{e1}{e2}) n2={n2} n1={tag}: {v:.4f} vs n1=3 {m:.4f}")

    # 3. Golden values, n1=4 hyperbolic deflector (r=1, 60 deg, n2=0).
    g = run(binp, 1.0, 60.0, 4.0, 0.0)
    for (plane, e1, e2), want in {("x", 2, 0): 16.089065361694130,
                                  ("x", 0, 1): 1.249367050523975,
                                  ("a", 2, 0): 31.488350359942770}.items():
        got = c(g, plane, e1, e2)
        if abs(got - want) > 1e-9:
            fails.append(f"golden ({plane}|{e1}{e2}): {got} != {want}")

    if fails:
        print("FAIL:\n  " + "\n  ".join(fails))
        sys.exit(1)
    print("PASS: symplecticity, n1=3 continuity, and n1=4 golden values OK")


if __name__ == "__main__":
    main()
