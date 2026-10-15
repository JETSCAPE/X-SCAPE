#!/usr/bin/env python3
"""EMD qhat at the best-fit lambda_t = 2.4 vs the default JETSCAPE MATTER+LBT,
full hydro, from the reruns with the extended jetana (tasks_x.txt).

Models (paired: same hydro events and hard scatterings):
  tune          MATTER (Type 5) + LBT (tables)
  lbtemd_l2p4   MATTER (Type 9) + LBT with the EMD radiative qhat, lambda_t = 2.4
Samples: LHC 0-10% (60 hydro events), 30-50% (35); RHIC 0-10% (80), 30-50% (40);
pp: LHC 16 x 25k, RHIC 8 x 25k events.

Writes
  x_compare.txt          paired ln(EMD/tune) per observable, mean +- s.e.
  fig_x_raa_e2c.png      jet R_AA and E2C (AA/pp per jet), four samples
  fig_x_newobs.png       x_J, D(z), rho(r), charged-jet mass at the LHC
The new observables (x_J, D(z), rho(r), charged-jet mass) use LHC-style cuts
(see jetana.cc) and are shown for the LHC only.
"""
import glob
import math
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from spec_tools import combine, read_spec

HERE = os.path.dirname(os.path.abspath(__file__))
SPEC = os.path.join(HERE, "spec")
SAMPLES = {  # key: (prefix, seeds, pp files, title, selection, jet spectrum, jet bins, lhc?)
    "lhc0010": ("lhcx", list(range(16, 46)) + list(range(331, 361)), "ppx_s*", "Pb–Pb 5.02 TeV 0–10%", "ls", "fulljet",
                ["40-60", "60-80", "80-120", "120-200"], True),
    "lhc3050": ("lhcx_c3050", list(range(116, 151)), "ppx_s*", "Pb–Pb 5.02 TeV 30–50%", "ls", "fulljet",
                ["40-60", "60-80", "80-120", "120-200"], True),
    "rhic0010": ("rhicx", list(range(56, 96)) + list(range(171, 211)), "rhic_ppx_s*", "Au–Au 200 GeV 0–10%", "rs", "sphenix",
                 ["15-20", "20-30", "30-60"], False),
    "rhic3050": ("rhicx_c3050", list(range(131, 171)), "rhic_ppx_s*", "Au–Au 200 GeV 30–50%", "rs", "sphenix",
                 ["15-20", "20-30", "30-60"], False),
}
MODELS = {"tune": ("default MATTER+LBT", "#1F5FA8", "o", "--"),
          "lbtemd_l2p4": (r"EMD $\hat q$, $\lambda_t=2.4$", "#C2410C", "s", "-")}
YR = {"rs": [(15, 30), (30, 60)], "ls": [(40, 100), (100, 316)]}
XJB = ["100-126", "126-158", "158-200", "200-251", "251-316"]
FFB = ["126-158", "158-200", "200-251", "251-316"]
CHB = ["60-80", "80-100", "100-120"]


# ---------------------------------------------------------------- helpers
def yield_in(d, obs, lo, hi):
    e = d[obs]
    m = (e[:-1] >= lo - 1e-9) & (e[1:] <= hi + 1e-9)
    return np.sum(d[obs + "_h"][m]) / d["sumw"]


def mean_of(d, key):
    e, h = d[key], d[key + "_h"]
    return np.sum(h * 0.5 * (e[:-1] + e[1:])) / np.sum(h)


def frac(d, key, lo, hi, centers=False):
    e, h = d[key], d[key + "_h"]
    x = np.sqrt(e[:-1] * e[1:]) if centers else 0.5 * (e[:-1] + e[1:])
    m = (x >= lo) & (x < hi)
    return np.sum(h[m]) / np.sum(h)


def observables(d, sel, jet, lhc):
    out = {}
    for lo, hi in YR[sel]:
        out[f"jet yield {lo}-{hi}"] = yield_in(d, jet, lo, hi)
    hlo = 10 if lhc else 5
    out[f"hadron yield {hlo}-{4 * hlo}"] = yield_in(d, "hadron", hlo, 4 * hlo)
    lo, hi = YR[sel][0][0], YR[sel][-1][1]
    out["R=0.6 / R=0.2 yield"] = yield_in(d, f"{sel}_specR06", lo, hi) / yield_in(d, f"{sel}_specR02", lo, hi)
    nk = len(d[f"{sel}_njet"]) - 1
    for k in range(nk):
        lab = f"{d[f'{sel}_njet'][k]:.0f}-{d[f'{sel}_njet'][k + 1]:.0f}"
        out[f"<mass> {lab}"] = mean_of(d, f"{sel}_mass{k}")
        out[f"<girth> {lab}"] = mean_of(d, f"{sel}_ang1{k}")
        out[f"<lambda2> {lab}"] = mean_of(d, f"{sel}_ang2{k}")
        out[f"E2C frac R_L<0.1 {lab}"] = frac(d, f"{sel}_eec{k}", 0., 0.1, centers=True)
    if lhc:
        for k, lab in enumerate(XJB):
            out[f"<x_J> {lab}"] = mean_of(d, f"xj{k}")
            out[f"pair fraction {lab}"] = d["xj_pair_h"][k] / d["xj_lead_h"][k]
        for k, lab in enumerate(FFB):
            nj = d["ff_njet_h"][k]
            e, h = d[f"ff_z{k}"], d[f"ff_z{k}_h"]
            z = np.sqrt(e[:-1] * e[1:])
            for zlo, zhi in ((0.01, 0.04), (0.04, 0.3), (0.3, 1.0)):
                out[f"D(z) {zlo}-{zhi}, jets {lab}"] = np.sum(h[(z >= zlo) & (z < zhi)]) / nj
        for rlo, rhi in ((0., 0.1), (0.1, 0.3), (0.3, 1.0)):
            out[f"rho(r) frac {rlo}-{rhi}"] = frac(d, "js_rho", rlo, rhi)
        for k, lab in enumerate(CHB):
            out[f"<ch. mass> {lab}"] = mean_of(d, f"chm_mass{k}")
    return out


def per_jet(d, nkey, k, key, logx=False):
    e = d[key]
    nj = d[nkey + "_h"][k]
    w = np.diff(np.log(e)) if logx else np.diff(e)
    return e, d[key + "_h"] / nj / w, np.sqrt(d[key + "_h2"]) / nj / w


def ratio(a, ea, b, eb):
    with np.errstate(divide="ignore", invalid="ignore"):
        r = a / b
        return r, np.abs(r) * np.hypot(ea / a, eb / b)


def raa(aa, pp, obs):
    ya, yp = aa[obs + "_h"] / aa["sumw"], pp[obs + "_h"] / pp["sumw"]
    with np.errstate(divide="ignore", invalid="ignore"):
        r = ya / yp
        e = np.abs(r) * np.sqrt(aa[obs + "_h2"] / aa[obs + "_h"] ** 2 + pp[obs + "_h2"] / pp[obs + "_h"] ** 2)
    return aa[obs], r, e


# ---------------------------------------------------------------- 1. paired differences
lines = []
P = lambda s="": (print(s), lines.append(s))
P("# ln(EMD lambda_t=2.4 / default MATTER+LBT), paired per hydro event: mean +- s.e. (significance)")
SIG = []
for skey, (pre, seeds, ppat, title, sel, jet, jb, lhc) in SAMPLES.items():
    diffs = {}
    for s in seeds:
        fe = os.path.join(SPEC, f"{pre}_lbtemd_l2p4_h{s}.spec")
        ft = os.path.join(SPEC, f"{pre}_tune_h{s}.spec")
        if not (os.path.exists(fe) and os.path.exists(ft)):
            continue
        oe, ot = observables(read_spec(fe), sel, jet, lhc), observables(read_spec(ft), sel, jet, lhc)
        for k in oe:
            if oe[k] > 0 and ot[k] > 0 and np.isfinite(oe[k]) and np.isfinite(ot[k]):
                diffs.setdefault(k, []).append(math.log(oe[k] / ot[k]))
    P(f"\n## {title} ({len(seeds)} hydro events)")
    for k, v in diffs.items():
        v = np.array(v)
        if len(v) < 5:
            continue
        mu, se = v.mean(), v.std(ddof=1) / math.sqrt(len(v))
        SIG.append((skey, k, mu, se))
        P(f"  {k:34s} {100 * mu:+7.2f} +- {100 * se:5.2f}%  {mu / se:+5.1f}σ  n={len(v)}{'  <--' if abs(mu) > 3 * se else ''}")
z = np.array([m / s for _, _, m, s in SIG if s > 0])
P(f"\n# {len(z)} paired tests: |z|>2: {np.sum(np.abs(z) > 2)} (expected ~{0.0455 * len(z):.1f}), |z|>3: {np.sum(np.abs(z) > 3)} (expected ~{0.0027 * len(z):.2f})")

# ---------------------------------------------------------------- load combined
D = {}
for skey, (pre, seeds, ppat, *_r) in SAMPLES.items():
    D[skey] = {"pp": combine(sorted(glob.glob(os.path.join(SPEC, ppat + ".spec"))))}
    for m in MODELS:
        fl = [os.path.join(SPEC, f"{pre}_{m}_h{s}.spec") for s in seeds]
        D[skey][m] = combine([f for f in fl if os.path.exists(f)])

# ---------------------------------------------------------------- 2. R_AA and E2C
fig, ax = plt.subplots(4, 4, figsize=(19, 14), gridspec_kw={"height_ratios": [2, 1, 2, 1]})
for i, skey in enumerate(SAMPLES):
    pre, seeds, ppat, title, sel, jet, jb, lhc = SAMPLES[skey]
    r0, c0 = (i // 2) * 2, (i % 2) * 2
    pp = D[skey]["pp"]
    # jet R_AA
    a, b = ax[r0, c0], ax[r0 + 1, c0]
    res = {}
    for m, (lab, c, mk, ls) in MODELS.items():
        e, r, er = raa(D[skey][m], pp, jet)
        x = 0.5 * (e[:-1] + e[1:])
        res[m] = (r, er)
        a.errorbar(x, r, yerr=er, fmt=mk, ls=ls, ms=4, color=c, label=lab)
    q, eq = ratio(*res["lbtemd_l2p4"], *res["tune"])
    b.errorbar(x, q, yerr=eq, fmt="s", color="#C2410C", ms=4)
    a.set_xscale("log"); b.set_xscale("log"); a.set_ylim(0, 1.2); b.set_ylim(0.8, 1.2); b.axhline(1, color="gray", lw=0.7)
    a.set_title(f"{title}: R=0.4 jet $R_{{AA}}$", fontsize=10); a.set_ylabel(r"$R_{AA}$"); b.set_ylabel("EMD / default")
    b.set_xlabel(r"jet $p_T$ [GeV]"); a.legend(fontsize=8)
    # E2C AA/pp, one jet bin
    k = 2 if lhc else 1
    a, b = ax[r0, c0 + 1], ax[r0 + 1, c0 + 1]
    e, yp, ep = per_jet(pp, f"{sel}_njet", k, f"{sel}_eec{k}", logx=True)
    x = np.sqrt(e[:-1] * e[1:]); ok = yp > 0.02 * np.nanmax(yp)
    res = {}
    for m, (lab, c, mk, ls) in MODELS.items():
        _, y, ey = per_jet(D[skey][m], f"{sel}_njet", k, f"{sel}_eec{k}", logx=True)
        r, er = ratio(y, ey, yp, ep); res[m] = (r, er)
        a.errorbar(x[ok], r[ok], yerr=er[ok], fmt=mk, ls=ls, ms=4, color=c)
    q, eq = ratio(*res["lbtemd_l2p4"], *res["tune"])
    b.errorbar(x[ok], q[ok], yerr=eq[ok], fmt="s", color="#C2410C", ms=4)
    for z_ in (a, b):
        z_.set_xscale("log"); z_.axhline(1, color="gray", lw=0.7)
    a.set_ylim(0.4, 2.2); b.set_ylim(0.8, 1.2)
    a.set_title(f"{title}: E2C AA/pp, jets {jb[k]} GeV", fontsize=10); b.set_xlabel(r"$R_L$")
fig.suptitle(r"Full hydro: EMD $\hat q$ at $\lambda_t=2.4$ vs default JETSCAPE MATTER+LBT (paired hydro events)", fontsize=13)
fig.tight_layout(rect=(0, 0, 1, 0.97)); fig.savefig(os.path.join(HERE, "fig_x_raa_e2c.png"), dpi=105)

# ---------------------------------------------------------------- 3. new observables (LHC)
fig, ax = plt.subplots(4, 4, figsize=(19, 14), gridspec_kw={"height_ratios": [2, 1, 2, 1]})
for i, skey in enumerate(("lhc0010", "lhc3050")):
    pre, seeds, ppat, title, sel, jet, jb, lhc = SAMPLES[skey]
    pp = D[skey]["pp"]
    r0 = 2 * i
    panels = []
    # (a) x_J, per-pair normalised, pT1 126-158 GeV
    k = 1
    e = pp[f"xj{k}"]; x = 0.5 * (e[:-1] + e[1:]); w = np.diff(e)
    def xjdist(d):
        n = d["xj_pair_h"][k]
        return d[f"xj{k}_h"] / n / w, np.sqrt(d[f"xj{k}_h2"]) / n / w
    panels.append(("x_J", x, xjdist, f"$x_J$ (per pair), $p_{{T,1}}$ {XJB[k]} GeV", False, (0, 4)))
    # (b) D(z) per jet, jets 126-158
    ef = pp["ff_z0"]; xf = np.sqrt(ef[:-1] * ef[1:])
    def ffdist(d):
        _, y, ey = per_jet(d, "ff_njet", 0, "ff_z0", logx=False)
        return y, ey
    panels.append(("D(z)", xf, ffdist, f"$D(z)$ per jet, jets {FFB[0]} GeV", True, None))
    # (c) rho(r), normalised to unit integral over r < 1
    er_ = pp["js_rho"]; xr = 0.5 * (er_[:-1] + er_[1:]); wr = np.diff(er_)
    def rho(d):
        h = d["js_rho_h"]; s = np.sum(h)
        return h / s / wr, np.sqrt(d["js_rho_h2"]) / s / wr
    panels.append(("rho", xr, rho, r"$\rho(r)$, jets $p_T>120$ GeV, $0.3<|\eta|<1.6$", False, None))
    # (d) charged-jet mass per jet, ch. jets 60-80
    em = pp["chm_mass0"]; xm = 0.5 * (em[:-1] + em[1:])
    def chm(d):
        _, y, ey = per_jet(d, "chm_njet", 0, "chm_mass0")
        return y, ey
    panels.append(("m", xm, chm, f"charged-jet mass per jet, $p_{{T,ch}}$ {CHB[0]} GeV", False, None))
    for c_, (name, xx, fn, lab, logx, _yl) in enumerate(panels):
        a, b = ax[r0, c_], ax[r0 + 1, c_]
        yp, ep = fn(pp)
        ok = np.isfinite(yp) & (yp > 0.01 * np.nanmax(yp))
        a.errorbar(xx[ok], yp[ok], yerr=ep[ok], fmt="k*", ms=5, label="pp")
        res = {}
        for m, (ml, c, mk, ls) in MODELS.items():
            y, ey = fn(D[skey][m]); res[m] = (y, ey)
            a.errorbar(xx[ok], y[ok], yerr=ey[ok], fmt=mk, ls=ls, ms=3.5, color=c, label=ml)
        q, eq = ratio(*res["lbtemd_l2p4"], *res["tune"])
        b.errorbar(xx[ok], q[ok], yerr=eq[ok], fmt="s", ms=4, color="#C2410C")
        b.axhline(1, color="gray", lw=0.7); b.set_ylim(0.8, 1.2); b.set_ylabel("EMD / default" if c_ == 0 else "")
        if logx:
            a.set_xscale("log"); b.set_xscale("log")
        if name in ("D(z)", "rho"):
            a.set_yscale("log")
        a.set_title(f"{title}: {lab}", fontsize=9)
        b.set_xlabel({"x_J": "$x_J$", "D(z)": "$z$", "rho": "$r$", "m": "$m_{ch}$ [GeV]"}[name])
        if c_ == 0:
            a.legend(fontsize=8)
fig.suptitle(r"Full hydro, LHC: EMD $\hat q$ at $\lambda_t=2.4$ vs default MATTER+LBT — dijet $x_J$, fragmentation, jet shape, charged-jet mass", fontsize=13)
fig.tight_layout(rect=(0, 0, 1, 0.97)); fig.savefig(os.path.join(HERE, "fig_x_newobs.png"), dpi=105)

open(os.path.join(HERE, "x_compare.txt"), "w").write("\n".join(lines) + "\n")
print("wrote x_compare.txt, fig_x_raa_e2c.png, fig_x_newobs.png")
