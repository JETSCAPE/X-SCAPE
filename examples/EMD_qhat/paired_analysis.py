#!/usr/bin/env python3
"""Paired EMD-vs-HTL comparison.

The EMD and HTL runs of a given seed share the hydro event and the Pythia
hard scatterings (Pythia is seeded independently of the energy loss), so
their fluctuations are strongly correlated. Here each observable is computed
per seed for both models, d = ln(EMD/HTL) is formed seed by seed, and its
mean and standard error over seeds give the difference with the correlation
properly accounted for. The 30-50% vs 0-10% comparison combines the two
independent paired means (different seeds).
"""
import math
import os

import numpy as np

from spec_tools import read_spec

HERE = os.path.dirname(os.path.abspath(__file__))
SPEC = os.path.join(HERE, "spec")
SAMPLES = {"rhic0010": ("rhic", range(56, 96), "rs", "sphenix"),
           "rhic3050": ("rhic_c3050", range(131, 171), "rs", "sphenix"),
           "lhc0010": ("lhc", range(16, 46), "ls", "fulljet"),
           "lhc3050": ("lhc_c3050", range(116, 151), "ls", "fulljet")}
RANGES = {"rs": [(15, 30), (30, 60)], "ls": [(40, 100), (100, 316)]}


def yield_in(d, obs, lo, hi):
    e = d[obs]
    m = (e[:-1] >= lo - 1e-9) & (e[1:] <= hi + 1e-9)
    return np.sum(d[obs + "_h"][m]) / d["sumw"]


def mean_of(d, sel, obs, k):
    e = d[f"{sel}_{obs}{k}"]
    h = d[f"{sel}_{obs}{k}_h"]
    return np.sum(h * 0.5 * (e[:-1] + e[1:])) / np.sum(h)


def observables(d, sel, jet):
    out = {}
    for lo, hi in RANGES[sel]:
        out[f"jet yield {lo}-{hi}"] = yield_in(d, jet, lo, hi)
    for R in ("R06",):
        lo, hi = RANGES[sel][0][0], RANGES[sel][-1][1]
        out["yield R=0.6 / R=0.2"] = yield_in(d, f"{sel}_spec{R}", lo, hi) / yield_in(d, f"{sel}_specR02", lo, hi)
    nk = len(d[f"{sel}_njet"]) - 1
    for k in range(nk):
        lab = f"{d[f'{sel}_njet'][k]:.0f}-{d[f'{sel}_njet'][k + 1]:.0f}"
        out[f"<mass> {lab}"] = mean_of(d, sel, "mass", k)
        out[f"<girth> {lab}"] = mean_of(d, sel, "ang1", k)
        out[f"<lambda2> {lab}"] = mean_of(d, sel, "ang2", k)
    return out


res = {}
print(f"{'sample':9s} {'observable':24s} {'EMD/HTL - 1':>14s} {'signif.':>8s}  n")
for skey, (pre, seeds, sel, jet) in SAMPLES.items():
    diffs = {}
    for s in seeds:
        fe = os.path.join(SPEC, f"{pre}_emdv_l256_h{s}.spec")
        fh = os.path.join(SPEC, f"{pre}_htlv_k16_h{s}.spec")
        if not (os.path.exists(fe) and os.path.exists(fh)):
            continue
        oe = observables(read_spec(fe), sel, jet)
        oh = observables(read_spec(fh), sel, jet)
        for k in oe:
            if oe[k] > 0 and oh[k] > 0 and np.isfinite(oe[k]) and np.isfinite(oh[k]):
                diffs.setdefault(k, []).append(math.log(oe[k] / oh[k]))
    res[skey] = {}
    for k, v in diffs.items():
        v = np.array(v)
        m, se = v.mean(), v.std(ddof=1) / math.sqrt(len(v))
        res[skey][k] = (m, se)
        flag = "  <--" if abs(m) > 3 * se else ""
        print(f"{skey:9s} {k:24s} {100 * m:+7.2f} +- {100 * se:4.2f}% {m / se:+7.1f}σ  {len(v)}{flag}")

print("\n# centrality dependence: [ln(EMD/HTL)](30-50%) - [ln(EMD/HTL)](0-10%)")
for sys_ in ("rhic", "lhc"):
    a, b = res[sys_ + "3050"], res[sys_ + "0010"]
    for k in a:
        if k.startswith("jet yield") and k in b:
            m = a[k][0] - b[k][0]
            se = math.hypot(a[k][1], b[k][1])
            print(f"{sys_:5s} {k:24s} {100 * m:+7.2f} +- {100 * se:4.2f}%  {m / se:+5.1f}σ")
