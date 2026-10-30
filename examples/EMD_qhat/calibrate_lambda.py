#!/usr/bin/env python3
"""chi^2 calibration of lambda_t for the default chain with the EMD qhat
(MATTER Type 9 + LBT qhat_model=1, model 'lbtemd') against LHC data.

Grid: lambda_t = 1.0, 1.5, 2.0, 3.0, 4.5, 6.0; the JETSCAPE tune as a reference.
Samples (paired hydro seeds): Pb-Pb 5.02 TeV 0-10% (seeds 16-45, 331-360) and
30-50% (116-150). Data: ATLAS R=0.4 jet R_AA and CMS charged-hadron R_AA at
both centralities.

The model R_AA (and its statistical error) is interpolated in ln pT to each
data point inside the model's range. Two chi^2 definitions:
  uncorr: data stat and sys added in quadrature, point by point (conservative)
  corr:   data sys fully correlated within a dataset, one nuisance parameter
          s per dataset (shift d + s * sigma_sys), profiled analytically, +s^2
chi^2(lambda_t) is fitted with a parabola in ln lambda_t through the lowest
grid point and its two neighbours; the minimum gives the best value and
Delta chi^2 = 1 its 1-sigma range. A minimum at the edge of the grid is
reported as a bound.
Writes lamcal.txt and fig_lamcal.png.
"""
import csv
import glob
import math
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from spec_tools import combine

HERE = os.path.dirname(os.path.abspath(__file__))
SPEC = os.path.join(HERE, "spec")
DATA = os.path.join(HERE, "data")
LAMS = [(1.0, "l1"), (1.5, "l1p5"), (2.0, "l2"), (3.0, "l3"), (4.5, "l4p5"), (6.0, "l6")]
SAMPLES = {"0010": ("lhc", list(range(16, 46)) + list(range(331, 361))),
           "3050": ("lhc_c3050", list(range(116, 151)))}
DATASETS = {  # (observable, label, file, columns x, y, stat, sys+, sys-, pT range)
    ("0010", "atlas"): ("fulljet", "ATLAS jets 0-10%", "ATLAS_jetRAA_R04_0-10_5020_ins1673184_t19.csv", (0, 3, 6, 4, 5), (100, 320)),
    ("3050", "atlas"): ("fulljet", "ATLAS jets 30-50%", "ATLAS_jetRAA_R04_30-50_5020_ins1673184_t22t23avg.csv", (0, 3, 6, 4, 5), (100, 320)),
    ("0010", "cms"): ("hadron", "CMS hadrons 0-10%", "CMS_hRAA_0-10_5020_ins1496050_t14.csv", (0, 3, 4, 6, 7), (10, 100)),
    ("3050", "cms"): ("hadron", "CMS hadrons 30-50%", "CMS_hRAA_30-50_5020_ins1496050_t11.csv", (0, 3, 4, 6, 7), (10, 100)),
}
COMBOS = {"ATLAS 0-10%": [("0010", "atlas")], "ATLAS 30-50%": [("3050", "atlas")],
          "ATLAS both": [("0010", "atlas"), ("3050", "atlas")],
          "CMS both": [("0010", "cms"), ("3050", "cms")],
          "ATLAS + CMS, both (main)": [("0010", "atlas"), ("3050", "atlas"), ("0010", "cms"), ("3050", "cms")]}


def hepdata(name, cols):
    rows = [r for r in csv.reader(l for l in open(os.path.join(DATA, name)) if not l.startswith("#")) if r][1:]
    a = np.array([[float(r[i]) for i in cols] for r in rows])
    a[:, 2:] = np.abs(a[:, 2:])
    return a  # x, y, stat, sys+, sys-


def raa(aa, pp, obs):
    ya, yp = aa[obs + "_h"] / aa["sumw"], pp[obs + "_h"] / pp["sumw"]
    with np.errstate(divide="ignore", invalid="ignore"):
        r = ya / yp
        e = np.abs(r) * np.sqrt(aa[obs + "_h2"] / aa[obs + "_h"] ** 2 + pp[obs + "_h2"] / pp[obs + "_h"] ** 2)
    x = 0.5 * (aa[obs][:-1] + aa[obs][1:])
    ok = np.isfinite(r) & np.isfinite(e) & (e > 0)
    return x[ok], r[ok], e[ok]


def load(prefix, tag, seeds):
    fl = [os.path.join(SPEC, f"{prefix}_{tag}_h{s}.spec") for s in seeds]
    fl = [f for f in fl if os.path.exists(f)]
    return (combine(fl) if fl else None), len(fl)


PP = combine(sorted(glob.glob(os.path.join(SPEC, "pp_s*.spec"))))
DAT = {k: hepdata(v[2], v[3]) for k, v in DATASETS.items()}


def residuals(model, key):
    """Model at the data points: (m, sm, d, sstat, ssys) arrays."""
    obs, _, _, _, (lo, hi) = DATASETS[key]
    x, r, e = raa(model, PP, obs)
    d = DAT[key]
    d = d[(d[:, 0] >= max(lo, x.min())) & (d[:, 0] <= min(hi, x.max()))]
    lx = np.log(x)
    m = np.interp(np.log(d[:, 0]), lx, r)
    sm = np.interp(np.log(d[:, 0]), lx, e)
    return m, sm, d[:, 1], d[:, 2], 0.5 * (d[:, 3] + d[:, 4])


def chi2(model_by_cent, keys):
    cu = cc = 0.
    n = 0
    for key in keys:
        m, sm, d, st, sy = residuals(model_by_cent[key[0]], key)
        cu += float(np.sum((m - d) ** 2 / (sm ** 2 + st ** 2 + sy ** 2)))
        # correlated sys: minimise sum (m - d - s sy)^2 / v + s^2 over s
        v = sm ** 2 + st ** 2
        s = float(np.sum((m - d) * sy / v) / (1 + np.sum(sy ** 2 / v)))
        cc += float(np.sum((m - d - s * sy) ** 2 / v) + s ** 2)
        n += len(d)
    return cu, cc, n


def parabola_min(lams, chis):
    x = np.log(np.array(lams)); y = np.array(chis)
    a, b, c = np.polyfit(x, y, 2)
    if a <= 0:
        return float("nan"), float("nan"), float("nan"), float("nan")
    x0 = -b / (2 * a); y0 = c - b * b / (4 * a); dx = 1 / math.sqrt(a)
    return math.exp(x0), math.exp(x0 - dx), math.exp(x0 + dx), y0


out = []
def p(s=""):
    print(s); out.append(s)


MOD = {}
for lam, tag in LAMS + [(None, "tune")]:
    MOD[tag] = {}
    for cent, (pre, seeds) in SAMPLES.items():
        MOD[tag][cent], n = load(pre, tag if tag == "tune" else f"lbtemd_{tag}", seeds)
        p(f"# {tag:5s} {cent}: {n} hydro events")
avail = [(l, t) for l, t in LAMS if all(MOD[t][c] is not None for c in SAMPLES)]

p("\n# chi^2 / n_points  (uncorrelated sys | correlated sys)")
RES = {}
for name, keys in COMBOS.items():
    p(f"\n## {name}")
    rows = []
    for lam, tag in avail:
        cu, cc, n = chi2(MOD[tag], keys)
        rows.append((lam, cu, cc, n))
        p(f"  lambda_t = {lam:3.1f}: {cu:7.1f} | {cc:7.1f}   / {n}")
    cu, cc, n = chi2(MOD["tune"], keys)
    p(f"  tune          : {cu:7.1f} | {cc:7.1f}   / {n}")
    RES[name] = rows
    if len(rows) >= 3:
        for lab, i in (("uncorr", 1), ("corr", 2)):
            k = int(np.argmin([r[i] for r in rows]))
            if k in (0, len(rows) - 1):
                p(f"  best ({lab}): at the grid edge, lambda_t {'<=' if k == 0 else '>='} {rows[k][0]:g}, chi2 = {rows[k][i]:.1f} / {rows[0][3]}")
                continue
            near = rows[k - 1:k + 2]
            b, lo, hi, cmin = parabola_min([r[0] for r in near], [r[i] for r in near])
            p(f"  best ({lab}): lambda_t = {b:.2f}  (1 sigma {lo:.2f}-{hi:.2f}), chi2_min = {cmin:.1f} / {rows[0][3]}")

open(os.path.join(HERE, "lamcal.txt"), "w").write("\n".join(out) + "\n")

# ---------------------------------------------------------------- figure
fig, ax = plt.subplots(2, 3, figsize=(17, 9))
cols = {1.0: "#F2A066", 1.5: "#E06A1F", 2.0: "#C2410C", 3.0: "#8A2E06", 4.5: "#5B1F04", 6.0: "#2E1002"}
for row, cent in enumerate(SAMPLES):
    for col, kind in enumerate(("atlas", "cms")):
        a = ax[row, col]; key = (cent, kind)
        d = DAT[key]
        a.errorbar(d[:, 0], d[:, 1], yerr=np.hypot(d[:, 2], 0.5 * (d[:, 3] + d[:, 4])), fmt="k*", ms=7, label=DATASETS[key][1])
        for lam, tag in avail:
            x, r, e = raa(MOD[tag][cent], PP, DATASETS[key][0])
            a.errorbar(x, r, yerr=e, fmt="s-", ms=3, lw=1, color=cols[lam], label=f"EMD, $\\lambda_t$={lam:g}")
        x, r, e = raa(MOD["tune"][cent], PP, DATASETS[key][0])
        a.errorbar(x, r, yerr=e, fmt="o--", ms=3, lw=1, color="#1F5FA8", label="tune")
        a.set_xscale("log"); a.set_ylim(0, 1.1); a.set_xlim(*(DATASETS[key][4][0] * 0.8, DATASETS[key][4][1] * 1.1))
        a.set_xlabel(r"$p_T$ [GeV]"); a.set_ylabel(r"$R_{AA}$"); a.legend(fontsize=7)
for i, (name, lab) in enumerate((("ATLAS + CMS, both (main)", "ATLAS + CMS"), ("ATLAS both", "ATLAS only"))):
    a = ax[i, 2]; rows = RES[name]
    if rows:
        L = [r[0] for r in rows]
        a.plot(L, [r[1] for r in rows], "o-", color="#C2410C", label="uncorrelated sys")
        a.plot(L, [r[2] for r in rows], "s--", color="#555B66", label="correlated sys")
        a.set_xscale("log"); a.set_xticks(L); a.set_xticklabels([f"{l:g}" for l in L])
    a.set_xlabel(r"$\lambda_t$"); a.set_ylabel(r"$\chi^2$"); a.set_title(f"{lab}, 0-10% and 30-50%"); a.legend(fontsize=8)
fig.tight_layout(); fig.savefig(os.path.join(HERE, "fig_lamcal.png"), dpi=110)
print("wrote lamcal.txt, fig_lamcal.png")
