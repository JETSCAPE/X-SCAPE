#!/usr/bin/env python3
"""E2C (two-point energy correlator) of R = 0.4 jets in full hydro, for the
EMD qhat in MATTER+LBT around the best-fit lambda_t, against the JETSCAPE tune.

The chi^2 fit to ATLAS jets (calibrate_lambda.py) gives lambda_t = 2.4
(2.0-4.0). Runs exist at lambda_t = 2 and 3 (LHC, both centralities, same
seeds) and 2.6 (LHC and RHIC 0-10%); lambda_t = 2.4 is interpolated bin by bin
between 2 and 3, linearly in ln lambda_t. RHIC has lambda_t = 2 (and 2.6 at
0-10%) only.

E2C per jet: sum_{i<j} pT_i pT_j / pT_jet^2 per unit ln R_L (charged
constituents, pT > 1 GeV, holes subtracted), from jetana. Writes
fig_e2c_bestfit.png (AA/pp per jet, all jet-pT bins, four samples),
fig_e2c_bestfit_dist.png (the distributions in one bin per sample) and
e2c_bestfit.txt (EMD(2.4)/tune per R_L bin).
"""
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
LHC010 = list(range(16, 46)) + list(range(331, 361))
SAMPLES = {  # key: (prefix, pp pattern, title, selection, jet-pT labels, {model: seeds})
    "lhc0010": ("lhc", "pp_s*", "Pb–Pb 5.02 TeV 0–10%", "ls", ["40-60", "60-80", "80-120", "120-200"],
                {"tune": LHC010, "lbtemd_l2": LHC010, "lbtemd_l3": LHC010, "lbtemd_l2p6": list(range(16, 46))}),
    "lhc3050": ("lhc_c3050", "pp_s*", "Pb–Pb 5.02 TeV 30–50%", "ls", ["40-60", "60-80", "80-120", "120-200"],
                {"tune": list(range(116, 151)), "lbtemd_l2": list(range(116, 151)), "lbtemd_l3": list(range(116, 151))}),
    "rhic0010": ("rhic", "rhic_pp_s*", "Au–Au 200 GeV 0–10%", "rs", ["15-20", "20-30", "30-60"],
                 {"tune": list(range(56, 96)) + list(range(171, 331)),
                  "lbtemd_l2": list(range(56, 96)) + list(range(171, 331)), "lbtemd_l2p6": list(range(56, 96))}),
    "rhic3050": ("rhic_c3050", "rhic_pp_s*", "Au–Au 200 GeV 30–50%", "rs", ["15-20", "20-30", "30-60"],
                 {"tune": list(range(131, 171)), "lbtemd_l2": list(range(131, 171))}),
}
STYLE = {"tune": ("JETSCAPE tune", "#1F5FA8", "o", "--"),
         "lbtemd_l2": (r"EMD, $\lambda_t=2$", "#F2A066", "s", ":"),
         "lbtemd_l2p4": (r"EMD, $\lambda_t=2.4$ (best fit, interp.)", "#C2410C", "D", "-"),
         "lbtemd_l2p6": (r"EMD, $\lambda_t=2.6$", "#E06A1F", "^", ":"),
         "lbtemd_l3": (r"EMD, $\lambda_t=3$", "#8A2E06", "v", ":")}
W24 = math.log(2.4 / 2.0) / math.log(3.0 / 2.0)


def per_jet(d, sel, k):
    e = d[f"{sel}_eec{k}"]
    nj = d[f"{sel}_njet_h"][k]
    w = np.diff(np.log(e))
    return e, d[f"{sel}_eec{k}_h"] / nj / w, np.sqrt(d[f"{sel}_eec{k}_h2"]) / nj / w


def ratio(a, ea, b, eb):
    with np.errstate(divide="ignore", invalid="ignore"):
        r = a / b
        return r, np.abs(r) * np.hypot(ea / a, eb / b)


out = []
D = {}
for skey, (pre, pppat, title, sel, jb, models) in SAMPLES.items():
    pp = combine(sorted(glob.glob(os.path.join(SPEC, pppat + ".spec"))))
    D[skey] = {"pp": (pp, 0)}
    for m, seeds in models.items():
        fl = [os.path.join(SPEC, f"{pre}_{m}_h{s}.spec") for s in seeds]
        fl = [f for f in fl if os.path.exists(f)]
        D[skey][m] = (combine(fl), len(fl))
    out.append(f"# {skey}: " + ", ".join(f"{m} {n}" for m, (_, n) in D[skey].items() if m != "pp") + " hydro events")


def curves(skey, k):
    """{model: (x, AA/pp, err)} for jet-pT bin k, with lambda_t = 2.4 interpolated."""
    pre, pppat, title, sel, jb, models = SAMPLES[skey]
    e, yp, ep = per_jet(D[skey]["pp"][0], sel, k)
    x = np.sqrt(e[:-1] * e[1:])
    res = {}
    for m in models:
        _, y, ey = per_jet(D[skey][m][0], sel, k)
        res[m] = ratio(y, ey, yp, ep)
    if "lbtemd_l2" in res and "lbtemd_l3" in res:
        (r2, e2), (r3, e3) = res["lbtemd_l2"], res["lbtemd_l3"]
        res["lbtemd_l2p4"] = (r2 + W24 * (r3 - r2), np.hypot((1 - W24) * e2, W24 * e3))
    ok = yp > 0.02 * np.nanmax(yp)
    return x, ok, res, (e, yp, ep)


# ---------------------------------------------------------------- AA/pp grid
order = ["tune", "lbtemd_l2", "lbtemd_l2p6", "lbtemd_l3", "lbtemd_l2p4"]
fig, ax = plt.subplots(4, 4, figsize=(19, 15), squeeze=False)
for row, skey in enumerate(SAMPLES):
    pre, pppat, title, sel, jb, models = SAMPLES[skey]
    for col in range(4):
        a = ax[row, col]
        if col >= len(jb):
            a.axis("off"); continue
        x, ok, res, _ = curves(skey, col)
        for j, m in enumerate([m for m in order if m in res]):
            lab, c, mk, ls = STYLE[m]
            r, er = res[m]
            q = ok & np.isfinite(r)
            xs = x[q] * (1 + 0.012 * (j - 2))
            a.errorbar(xs, r[q], yerr=er[q], fmt=mk, ls=ls, ms=4 if m != "lbtemd_l2p4" else 5,
                       lw=2.2 if m == "lbtemd_l2p4" else 1.0, color=c, label=lab if (row == 0 and col == 0) or (row == 2 and col == 0) else None)
        a.axhline(1, color="gray", lw=0.7); a.set_xscale("log"); a.set_ylim(0.4, 2.2)
        a.set_title(f"{title}, jets {jb[col]} GeV", fontsize=10)
        a.set_xlabel(r"$R_L$"); a.set_ylabel("E2C  AA/pp (per jet)" if col == 0 else "")
ax[0, 0].legend(fontsize=8, loc="upper left"); ax[2, 0].legend(fontsize=8, loc="upper left")
fig.suptitle(r"E2C of $R=0.4$ jets, full hydro: EMD $\hat q$ in MATTER+LBT around the best-fit $\lambda_t=2.4$ vs the JETSCAPE tune", fontsize=13)
fig.tight_layout(rect=(0, 0, 1, 0.97)); fig.savefig(os.path.join(HERE, "fig_e2c_bestfit.png"), dpi=110)

# ---------------------------------------------------------------- distributions, one bin per sample
fig, ax = plt.subplots(2, 4, figsize=(19, 8), gridspec_kw={"height_ratios": [2, 1]}, sharex="col")
show = {"lhc0010": 2, "lhc3050": 2, "rhic0010": 1, "rhic3050": 1}
for col, skey in enumerate(SAMPLES):
    pre, pppat, title, sel, jb, models = SAMPLES[skey]
    k = show[skey]
    x, ok, res, (e, yp, ep) = curves(skey, k)
    a, b = ax[0, col], ax[1, col]
    a.errorbar(x[ok], yp[ok], yerr=ep[ok], fmt="k*", ms=6, label="pp")
    for m in [m for m in order if m in res]:
        lab, c, mk, ls = STYLE[m]
        r, er = res[m]
        q = ok & np.isfinite(r)
        a.plot(x[q], (r * yp)[q], ls=ls if m != "lbtemd_l2p4" else "-", color=c, lw=2.2 if m == "lbtemd_l2p4" else 1.2, label=lab)
        if m != "tune":
            rt, ert = res["tune"]
            rr, err = ratio(r, er, rt, ert)
            qq = q & np.isfinite(rr)
            b.errorbar(x[qq], rr[qq], yerr=err[qq], fmt=mk, ls=ls, ms=4, color=c, lw=2 if m == "lbtemd_l2p4" else 1)
    a.set_xscale("log"); a.set_title(f"{title}, jets {jb[k]} GeV", fontsize=10)
    a.set_ylabel(r"$(1/N_{jet})\, d\Sigma_{E2C}/d\ln R_L$" if col == 0 else "")
    b.axhline(1, color="gray", lw=0.7); b.set_ylim(0.7, 1.3); b.set_xlabel(r"$R_L$")
    b.set_ylabel("EMD / tune" if col == 0 else "")
    if col == 0:
        a.legend(fontsize=8)
fig.tight_layout(); fig.savefig(os.path.join(HERE, "fig_e2c_bestfit_dist.png"), dpi=110)

# ---------------------------------------------------------------- table: EMD(2.4) / tune
out.append("\n# E2C ratio EMD(lambda_t = 2.4, interpolated) / tune, per jet, in R_L ranges (LHC)")
for skey in ("lhc0010", "lhc3050"):
    pre, pppat, title, sel, jb, models = SAMPLES[skey]
    for k, lab in enumerate(jb):
        x, ok, res, (e, yp, ep) = curves(skey, k)
        r, er = res["lbtemd_l2p4"]; rt, ert = res["tune"]
        rr, err = ratio(r, er, rt, ert)
        row = []
        for lo, hi in ((0.0, 0.03), (0.03, 0.1), (0.1, 0.4)):
            q = ok & np.isfinite(rr) & (x >= lo) & (x < hi) & (err > 0)
            if q.any():
                w = 1 / err[q] ** 2
                row.append(f"R_L {lo:.2f}-{hi:.2f}: {np.sum(rr[q] * w) / np.sum(w):.3f}+-{1 / math.sqrt(np.sum(w)):.3f}")
        out.append(f"  {skey} {lab:>8s}  " + "   ".join(row))
open(os.path.join(HERE, "e2c_bestfit.txt"), "w").write("\n".join(out) + "\n")
print("\n".join(out))
print("wrote fig_e2c_bestfit.png, fig_e2c_bestfit_dist.png, e2c_bestfit.txt")
