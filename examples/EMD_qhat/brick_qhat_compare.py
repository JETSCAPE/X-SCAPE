#!/usr/bin/env python3
"""Brick energy loss for the different qhat(T) choices (pure Python, no numpy).

For each brick run: the energy lost by the leading parton and within R = 0.4
of it (holes subtracted), relative to the vacuum run with the same gun, with
the statistical error of the mean. Writes pgfplots tables to figs_brick/ and
prints a summary. Run sets, in bricks/ (made by make_configs.py bricks and
run_bricks.sh; 100 GeV quark unless stated, L = 4 fm):

  tscan/   T scan, lambda_t = 1: tune (MATTER Type 5 + LBT tables),
           emd (MATTER Type 8 + MyEloss), emdv (MATTER Type 9 + MyEloss)
  lscan/   lambda_t = 4, 9, 16 (Type 9 + MyEloss) and the conformal control
           emdc_l9 (Q(T) frozen at T = 0.3 GeV)
  diag/    path length L = 1-8 fm at T = 0.3 GeV, tune vs emdv lambda_t = 9
  htlscan/ JETSCAPE's HTL qhat in the same chain (MATTER Type 5 + MyEloss
           qhat_model=1, both x K): K = 1, 3 vs T, and K = 3 vs L
  opt1/    route B (MATTER Type 9 + LBT with EMD radiative qhat),
           lambda_t = 1-8, quark and gluon, T = 0.3 GeV
"""
import math
import os

from brick_ana import events, kin

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "figs_brick")
os.makedirs(OUT, exist_ok=True)


def per_event(path, R=0.4):
    """Per-event leading-parton pT and cone pT (as in brick_ana.analyse)."""
    lead, cone = [], []
    for ev in events(path):
        parts = [(st, *kin(px, py, pz)) for pid, st, E, px, py, pz in ev]
        pos = [p for p in parts if p[0] >= 0]
        l = max(pos, key=lambda p: p[1])
        c = 0.
        for st, pt, eta, phi in parts:
            dphi = abs((phi - l[3] + math.pi) % (2 * math.pi) - math.pi)
            if math.hypot(eta - l[2], dphi) < R:
                c += pt if st >= 0 else (-pt if st == -1 else 0.)
        lead.append(l[1]); cone.append(c)
    return lead, cone


def mean_err(x):
    n = len(x)
    m = sum(x) / n
    v = sum((a - m) ** 2 for a in x) / (n - 1)
    return m, math.sqrt(v / n)


CACHE = {}


def stats(path):
    if path not in CACHE:
        lead, cone = per_event(path)
        CACHE[path] = (mean_err(lead), mean_err(cone), len(lead))
    return CACHE[path]


def loss(path, vac):
    """(dlead, err, dcone, err, n) relative to the vacuum run."""
    (ml, el), (mc, ec), n = stats(path)
    (vl, evl), (vc, evc), _ = stats(vac)
    return vl - ml, math.hypot(el, evl), vc - mc, math.hypot(ec, evc), n


def f(d, name):
    return os.path.join(HERE, "bricks", d, name + "_final_state_partons.dat")


def write(name, header, rows):
    with open(os.path.join(OUT, name), "w") as fh:
        fh.write(" ".join(header) + "\n")
        for r in rows:
            fh.write(" ".join(f"{x:.5g}" for x in r) + "\n")


VAC = f("tscan", "vac")
TS = [0.18, 0.2, 0.25, 0.3, 0.4, 0.5]

# ---------------------------------------------------------------- 1. T scan, lambda_t = 1
print("# T scan, 100 GeV quark, L = 4 fm, lambda_t = 1 [GeV lost: lead, cone]")
for m in ("tune", "emd", "emdv"):
    rows = []
    for T in TS:
        p = f("tscan", f"{m}_T{T:g}")
        if os.path.exists(p):
            dl, el, dc, ec, n = loss(p, VAC)
            rows.append((T, dl, el, dc, ec))
            print(f"  {m:5s} T={T:<5g} lead {dl:6.2f}+-{el:4.2f}  cone {dc:6.2f}+-{ec:4.2f}  n={n}")
    write(f"tscan_{m}.dat", ["T", "dlead", "elead", "dcone", "econe"], rows)

# ---------------------------------------------------------------- 2. lambda_t scan + conformal
print("\n# lambda_t scan (Type 9 + MyEloss) and conformal control")
for m in ("emdv_l4", "emdv_l9", "emdv_l16", "emdc_l9"):
    rows = []
    for T in (0.2, 0.3, 0.4, 0.5):
        p = f("lscan", f"{m}_T{T:g}")
        if os.path.exists(p):
            dl, el, dc, ec, n = loss(p, VAC)
            rows.append((T, dl, el, dc, ec))
            print(f"  {m:9s} T={T:<4g} lead {dl:6.2f}+-{el:4.2f}  cone {dc:6.2f}+-{ec:4.2f}")
    write(f"lscan_{m}.dat", ["T", "dlead", "elead", "dcone", "econe"], rows)

# ratio to the tune at the same T (lambda_t = 9 and conformal)
for m in ("emdv_l9", "emdc_l9"):
    rows = []
    for T in (0.2, 0.3, 0.4, 0.5):
        p, t = f("lscan", f"{m}_T{T:g}"), f("tscan", f"tune_T{T:g}")
        if os.path.exists(p) and os.path.exists(t):
            dl, el, dc, ec, _ = loss(p, VAC)
            tl, etl, tc, etc, _ = loss(t, VAC)
            rows.append((T, dl / tl, dl / tl * math.hypot(el / dl, etl / tl),
                         dc / tc, dc / tc * math.hypot(ec / dc, etc / tc)))
    write(f"ratio_{m}.dat", ["T", "rlead", "erlead", "rcone", "ercone"], rows)

# ---------------------------------------------------------------- 3. path length
print("\n# path length, T = 0.3 GeV: tune vs emdv lambda_t = 9")
for m, l4 in (("tune", f("tscan", "tune_T0.3")), ("emdv", f("lscan", "emdv_l9_T0.3"))):
    rows = []
    for L in (1, 2, 3, 4, 6, 8):
        p = l4 if L == 4 else f("diag", f"L{L}_{m}")
        if os.path.exists(p):
            dl, el, dc, ec, n = loss(p, VAC)
            rows.append((L, dl, el, dc, ec))
            print(f"  {m:5s} L={L}  lead {dl:6.2f}+-{el:4.2f}  cone {dc:6.2f}+-{ec:4.2f}  n={n}")
    write(f"lscan_L_{m}.dat", ["L", "dlead", "elead", "dcone", "econe"], rows)

# ---------------------------------------------------------------- 4. route B, EMD qhat in LBT
print("\n# route B: MATTER Type 9 + LBT with EMD radiative qhat, T = 0.3 GeV, L = 4 fm")
for pid, lab in ((1, "quark"), (21, "gluon")):
    vac = f("opt1", f"vac_pid{pid}")
    rows = []
    for lam in (1, 2, 4, 8):
        p = f("opt1", f"lbtemd_l{lam}_pid{pid}")
        dl, el, dc, ec, n = loss(p, vac)
        rows.append((lam, dl, el, dc, ec))
        print(f"  {lab} lambda_t={lam}  lead {dl:6.2f}+-{el:4.2f}  cone {dc:6.2f}+-{ec:4.2f}")
    write(f"opt1_{lab}.dat", ["lam", "dlead", "elead", "dcone", "econe"], rows)
    dl, el, dc, ec, _ = loss(f("opt1", f"tune_pid{pid}"), vac)
    print(f"  {lab} tune          lead {dl:6.2f}+-{el:4.2f}  cone {dc:6.2f}+-{ec:4.2f}")
    write(f"opt1_{lab}_tune.dat", ["lam", "dlead", "elead", "dcone", "econe"],
          [(0.7, dl, el, dc, ec), (11, dl, el, dc, ec)])

# ---------------------------------------------------------------- 4b. HTL qhat in the same kernel
# htlscan/: MATTER Type 5 (HTL x virtuality factor) + MyEloss qhat_model=1 (HTL),
# both times qhat_K = K. At lambda_t = 1 EMD and HTL have the same size, so
# HTL K = 1 pairs with EMD lambda_t = 1 and K = 3 with lambda_t = 9 (sqrt 9 = 3).
print("\n# HTL qhat in MATTER + MyEloss")
for K in (1, 3):
    rows = []
    for T in (0.2, 0.3, 0.4, 0.5):
        dl, el, dc, ec, n = loss(f("htlscan", f"htl_k{K}_T{T:g}"), VAC)
        rows.append((T, dl, el, dc, ec))
        print(f"  HTL K={K} T={T:<4g} lead {dl:6.2f}+-{el:4.2f}  cone {dc:6.2f}+-{ec:4.2f}  n={n}")
    write(f"htl_k{K}.dat", ["T", "dlead", "elead", "dcone", "econe"], rows)
rows = []
for L in (1, 2, 3, 4, 6, 8):
    p = f("htlscan", "htl_k3_T0.3" if L == 4 else f"htl_k3_L{L}")
    dl, el, dc, ec, n = loss(p, VAC)
    rows.append((L, dl, el, dc, ec))
    print(f"  HTL K=3 L={L}  lead {dl:6.2f}+-{el:4.2f}  cone {dc:6.2f}+-{ec:4.2f}  n={n}")
write("lscan_L_htl.dat", ["L", "dlead", "elead", "dcone", "econe"], rows)


def ratio_rows(pairs, key):
    out = []
    for x, a, b in pairs:
        dl, el, dc, ec, _ = loss(a, VAC)
        tl, etl, tc, etc, _ = loss(b, VAC)
        out.append((x, dl / tl, dl / tl * math.hypot(el / dl, etl / tl),
                    dc / tc, dc / tc * math.hypot(ec / dc, etc / tc)))
    return out


# HTL / tune, and EMD / HTL at matched size (same kernel, same seed)
write("ratio_htl_k3.dat", ["T", "rlead", "erlead", "rcone", "ercone"],
      ratio_rows([(T, f("htlscan", f"htl_k3_T{T:g}"), f("tscan", f"tune_T{T:g}")) for T in (0.2, 0.3, 0.4, 0.5)], "T"))
print("\n# EMD / HTL in MATTER + MyEloss [lead, cone]")
for name, pairs in (
        ("emd_over_htl_l9k3", [(T, f("lscan", f"emdv_l9_T{T:g}"), f("htlscan", f"htl_k3_T{T:g}")) for T in (0.2, 0.3, 0.4, 0.5)]),
        ("emd_over_htl_l1k1", [(T, f("tscan", f"emdv_T{T:g}"), f("htlscan", f"htl_k1_T{T:g}")) for T in (0.2, 0.3, 0.4, 0.5)]),
        ("emd_over_htl_L", [(L, f("lscan", "emdv_l9_T0.3") if L == 4 else f("diag", f"L{L}_emdv"),
                             f("htlscan", "htl_k3_T0.3" if L == 4 else f"htl_k3_L{L}")) for L in (1, 2, 3, 4, 6, 8)])):
        rows = ratio_rows(pairs, "x")
        write(f"{name}.dat", ["x", "rlead", "erlead", "rcone", "ercone"], rows)
        for r in rows:
            print(f"  {name:18s} x={r[0]:<4g} lead {r[1]:5.2f}+-{r[2]:4.2f}  cone {r[3]:5.2f}+-{r[4]:4.2f}")

# ---------------------------------------------------------------- 5. qhat/T^3 curves
def emd(T):
    th = math.tanh((1000 * T - 133.6481) / 50.3029)
    return 0.5 * 3.98564 * (1 - th) + 0.5 * (1 + th) * (8.10975 + 4.10263 * math.exp(-(1000 * T - 133.6481) / 254.864))


def htl(T, E, a=0.3, nf=3):
    s = max(2 * E * T, 1.)
    l2 = math.exp(-12 * math.pi / ((33 - 2 * nf) * a))
    ar = 12 * math.pi / ((33 - 2 * nf) * math.log(s / l2)) if s >= 1 else a
    mD2 = a * 4 * math.pi * T * T * (6 + nf) / 6
    return max(0., (3 * 50.4864 / math.pi) * ar * a * math.log(s / mD2))


write("qhat_curves.dat", ["T", "emd1", "emd2", "emd9", "conf9", "htl20", "htl50", "htl100"],
      [(T, emd(T), emd(T) * 2 ** .5, emd(T) * 3, emd(0.3) * 3, htl(T, 20), htl(T, 50), htl(T, 100))
       for T in [0.15 + 0.005 * i for i in range(71)]])
print("\nwrote", OUT)
