#!/usr/bin/env python3
"""Analysis of the ACCRE production (after the shower-graph fix).

Models (all normalised to the ATLAS 0-10% jet R_AA at 5.02 TeV):
  tune      MATTER+LBT (JETSCAPE)
  htlv_k16  MATTER(Type 5)+MyEloss, HTL qhat, K = 16.4
  emdv_l256 MATTER(Type 9)+MyEloss, EMD qhat, lambda_t = 256
Samples: RHIC Au-Au 200 GeV and LHC Pb-Pb 5.02 TeV, 0-10% and 30-50%.

Writes fig_data_0010.png, fig_data_3050.png, fig_rhic.png, fig_centrality.png,
fig_sub_<sample>.png, fig_radius.png and prints all tables, including the
EMD-vs-HTL significance of every observable.
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

MODELS = {"tune": ("MATTER+LBT (JETSCAPE)", "#1f5fa8", "o"),
          "htlv_k16": (r"MyEloss, HTL $\hat q$ (K=16.4)", "#15803d", "D"),
          "emdv_l256": (r"MyEloss, EMD $\hat q$ ($\lambda_t$=256)", "#c2410c", "s")}
SAMPLES = {  # key: (prefix, seeds, pp pattern, label, substructure selection, jet-pT labels)
    "rhic0010": ("rhic", range(56, 96), "rhic_pp_s*", "Au-Au 200 GeV 0-10%", "rs", ["15-20", "20-30", "30-60"]),
    "rhic3050": ("rhic_c3050", range(131, 171), "rhic_pp_s*", "Au-Au 200 GeV 30-50%", "rs", ["15-20", "20-30", "30-60"]),
    "lhc0010": ("lhc", range(16, 46), "pp_s*", "Pb-Pb 5.02 TeV 0-10%", "ls", ["40-60", "60-80", "80-120", "120-200"]),
    "lhc3050": ("lhc_c3050", range(116, 151), "pp_s*", "Pb-Pb 5.02 TeV 30-50%", "ls", ["40-60", "60-80", "80-120", "120-200"]),
}
SIG = []  # (sample, observable, bin, chi2, n)


def load(prefix, seeds):
    fl = [os.path.join(SPEC, f"{prefix}_h{s}.spec") for s in seeds]
    fl = [f for f in fl if os.path.exists(f)]
    return combine(fl), len(fl)


def raa(aa, pp, obs):
    ya, yp = aa[obs + "_h"] / aa["sumw"], pp[obs + "_h"] / pp["sumw"]
    with np.errstate(divide="ignore", invalid="ignore"):
        r = ya / yp
        e = np.abs(r) * np.sqrt(aa[obs + "_h2"] / aa[obs + "_h"] ** 2 + pp[obs + "_h2"] / pp[obs + "_h"] ** 2)
    return aa[obs], r, e


def ratio(a, ea, b, eb):
    with np.errstate(divide="ignore", invalid="ignore"):
        r = a / b
        return r, np.abs(r) * np.sqrt((ea / a) ** 2 + (eb / b) ** 2)


def chi2_diff(a, ea, b, eb):
    m = np.isfinite(a) & np.isfinite(b) & (ea > 0) & (eb > 0)
    return float(np.sum((a - b)[m] ** 2 / (ea[m] ** 2 + eb[m] ** 2))), int(m.sum())


def pvalue(c, n):
    if n <= 0:
        return float("nan")
    z = ((c / n) ** (1 / 3) - (1 - 2 / (9 * n))) / math.sqrt(2 / (9 * n))
    return 0.5 * math.erfc(z / math.sqrt(2))


def hepdata(name, xcol, ycol, stat, sysp, sysm):
    rows = [r for r in csv.reader(l for l in open(os.path.join(DATA, name)) if not l.startswith("#")) if r][1:]
    a = np.array([[float(r[i]) for i in (xcol, ycol, stat, sysp, sysm)] for r in rows])
    a[:, 2:] = np.abs(a[:, 2:])
    return a


# ---------------------------------------------------------------- load
D = {}
for skey, (pre, seeds, pppat, *_rest) in SAMPLES.items():
    pp = combine(sorted(glob.glob(os.path.join(SPEC, pppat + ".spec"))))
    D[skey] = {"pp": pp}
    for m in MODELS:
        D[skey][m] = load(f"{pre}_{m}", seeds)
print("hydro events per model:", {k: {m: D[k][m][1] for m in MODELS} for k in D})

# ---------------------------------------------------------------- (1) LHC data
DATASETS = {
    "0010": [("hadron", "CMS charged hadrons, |η|<1", hepdata("CMS_hRAA_0-10_5020_ins1496050_t14.csv", 0, 3, 4, 6, 7), (10, 160)),
             ("fulljet", "ATLAS R=0.4 jets, |y|<2.8", hepdata("ATLAS_jetRAA_R04_0-10_5020_ins1673184_t19.csv", 0, 3, 6, 4, 5), (30, 320)),
             ("chjet", "ALICE R=0.4 charged jets, |η|<0.5", hepdata("ALICE_chjetRAA_R04_0-10_5020_ins2637686_Figure_4b_top_ML.csv", 0, 1, 2, 4, 5), (30, 140))],
    "3050": [("hadron", "CMS charged hadrons, |η|<1", hepdata("CMS_hRAA_30-50_5020_ins1496050_t11.csv", 0, 3, 4, 6, 7), (10, 160)),
             ("fulljet", "ATLAS R=0.4 jets, |y|<2.8", hepdata("ATLAS_jetRAA_R04_30-50_5020_ins1673184_t22t23avg.csv", 0, 3, 6, 4, 5), (30, 320)),
             ("chjet", "ALICE R=0.4 charged jets, |η|<0.5", hepdata("ALICE_chjetRAA_R04_30-50_5020_ins2637686_Figure_4b_middle_ML.csv", 0, 1, 2, 4, 5), (30, 140))],
}
print("\n# chi2/n against LHC data (data stat+sys, model stat)")
for cent, sets in DATASETS.items():
    S = D["lhc" + cent]
    fig, ax = plt.subplots(2, 3, figsize=(15, 8), gridspec_kw={"height_ratios": [2.2, 1]}, sharex="col")
    for col, (obs, title, d, xr) in enumerate(sets):
        a, b = ax[0, col], ax[1, col]
        sel = (d[:, 0] >= xr[0]) & (d[:, 0] <= xr[1]); d = d[sel]
        a.fill_between(d[:, 0], d[:, 1] - d[:, 4], d[:, 1] + d[:, 3], color="0.85", step="mid")
        a.errorbar(d[:, 0], d[:, 1], yerr=d[:, 2], fmt="k*", ms=6, label="data")
        for i, (m, (lab, c, mk)) in enumerate(MODELS.items()):
            e, r, er = raa(S[m][0], S["pp"], obs)
            x = 0.5 * (e[:-1] + e[1:]) * (1 + 0.02 * (i - 1))
            ok = np.isfinite(r) & (x >= xr[0]) & (x <= xr[1])
            a.errorbar(x[ok], r[ok], yerr=er[ok], fmt=mk, ms=4, color=c, label=f"{lab} [{S[m][1]}]")
            dat = np.interp(x, d[:, 0], d[:, 1], left=np.nan, right=np.nan)
            dse = np.interp(x, d[:, 0], np.hypot(d[:, 2], 0.5 * (d[:, 3] + d[:, 4])), left=np.nan, right=np.nan)
            q = ok & np.isfinite(dat)
            b.errorbar(x[q], (r / dat)[q], yerr=(er / dat)[q], fmt=mk, ms=4, color=c)
            if q.any():
                chi = np.mean(((r - dat) / np.hypot(er, dse))[q] ** 2)
                print(f"  {cent} {obs:8s} {m:10s} n={q.sum():2d} chi2/n = {chi:5.2f}")
        a.set_title(title); a.set_xscale("log"); a.set_ylim(0, 1.2); a.set_ylabel(r"$R_{AA}$")
        b.axhspan(0.9, 1.1, color="0.92"); b.axhline(1, color="gray", lw=0.8); b.set_ylim(0.4, 1.8)
        b.set_xlabel(r"$p_T$ [GeV]"); b.set_ylabel("model/data")
    ax[0, 0].legend(fontsize=7)
    fig.suptitle(f"Pb-Pb 5.02 TeV {cent[:2]}-{cent[2:]}%")
    fig.tight_layout(); fig.savefig(os.path.join(HERE, f"fig_data_{cent}.png"), dpi=120)

# ---------------------------------------------------------------- (2) RHIC, (3) centrality
fig, ax = plt.subplots(1, 3, figsize=(16, 4.6))
for col, skey in enumerate(("rhic0010", "rhic3050")):
    for i, (m, (lab, c, mk)) in enumerate(MODELS.items()):
        e, r, er = raa(D[skey][m][0], D[skey]["pp"], "sphenix")
        x = 0.5 * (e[:-1] + e[1:]) * (1 + 0.015 * (i - 1))
        ax[col].errorbar(x, r, yerr=er, fmt=mk, color=c, ms=4, label=f"{lab} [{D[skey][m][1]}]")
    ax[col].set_title(f"{SAMPLES[skey][3]}: R=0.4 jets, |η|<0.7")
    ax[col].set_xscale("log"); ax[col].set_ylim(0, 1.1); ax[col].set_xlabel(r"jet $p_T$ [GeV]"); ax[col].set_ylabel(r"$R_{AA}$")
ax[0].legend(fontsize=7)
for skey, mk2 in (("rhic0010", "s"), ("rhic3050", "o"), ("lhc0010", "^"), ("lhc3050", "v")):
    obs = "sphenix" if skey.startswith("rhic") else "fulljet"
    e, ra, ea = raa(D[skey]["emdv_l256"][0], D[skey]["pp"], obs)
    _, rh, eh = raa(D[skey]["htlv_k16"][0], D[skey]["pp"], obs)
    q, eq = ratio(ra, ea, rh, eh)
    ax[2].errorbar(0.5 * (e[:-1] + e[1:]), q, yerr=eq, fmt=mk2, ms=4, label=SAMPLES[skey][3])
    c2, n = chi2_diff(ra, ea, rh, eh); SIG.append((skey, "jet R_AA", "all", c2, n))
ax[2].axhline(1, color="gray", lw=0.8); ax[2].set_xscale("log"); ax[2].set_ylim(0.6, 1.5)
ax[2].set_xlabel(r"jet $p_T$ [GeV]"); ax[2].set_ylabel(r"$R_{AA}^{EMD}/R_{AA}^{HTL}$"); ax[2].legend(fontsize=7)
ax[2].set_title("EMD vs HTL (same MyEloss chain)")
fig.tight_layout(); fig.savefig(os.path.join(HERE, "fig_rhic.png"), dpi=120)

print("\n# centrality dependence R_AA(30-50%) / R_AA(0-10%)")
fig, ax = plt.subplots(1, 3, figsize=(16, 4.6))
panels = [("rhic", "sphenix", "RHIC: R=0.4 jets, |η|<0.7"), ("lhc", "fulljet", "LHC: R=0.4 jets, |y|<2.8"),
          ("lhc", "hadron", "LHC: charged hadrons, |η|<1")]
for col, (sys_, obs, title) in enumerate(panels):
    cr = {}
    for i, (m, (lab, c, mk)) in enumerate(MODELS.items()):
        e, r3, e3 = raa(D[sys_ + "3050"][m][0], D[sys_ + "3050"]["pp"], obs)
        _, r0, e0 = raa(D[sys_ + "0010"][m][0], D[sys_ + "0010"]["pp"], obs)
        q, eq = ratio(r3, e3, r0, e0)
        cr[m] = (q, eq)
        x = 0.5 * (e[:-1] + e[1:]) * (1 + 0.015 * (i - 1))
        ax[col].errorbar(x, q, yerr=eq, fmt=mk, color=c, ms=4, label=lab)
        ok = np.isfinite(q) & (eq > 0)
        w = 1 / eq[ok] ** 2
        print(f"  {sys_:4s} {obs:8s} {m:10s} weighted mean = {np.sum(q[ok] * w) / np.sum(w):.3f} +- {1 / math.sqrt(np.sum(w)):.3f}")
    c2, n = chi2_diff(*cr["emdv_l256"], *cr["htlv_k16"])
    SIG.append((sys_, f"R_AA(30-50)/R_AA(0-10) {obs}", "all", c2, n))
    ax[col].axhline(1, color="gray", lw=0.8); ax[col].set_xscale("log"); ax[col].set_ylim(0.5, 2.5)
    ax[col].set_title(title + f"   EMD-HTL $\\chi^2$/n={c2:.1f}/{n}", fontsize=9)
    ax[col].set_xlabel(r"$p_T$ [GeV]"); ax[col].set_ylabel(r"$R_{AA}(30$-$50\%)/R_{AA}(0$-$10\%)$")
ax[0].legend(fontsize=7)
fig.tight_layout(); fig.savefig(os.path.join(HERE, "fig_centrality.png"), dpi=120)

# ---------------------------------------------------------------- (4) substructure
OBS = [("eec", "EEC", True), ("mass", "jet mass", False), ("ang05", "λ0.5", False), ("ang1", "girth λ1", False),
       ("ang2", "λ2", False)]


def per_jet(d, sel, obs, k):
    edges = d[f"{sel}_{obs}{k}"]
    nj = d[f"{sel}_njet_h"][k]
    w = np.diff(np.log(edges)) if obs == "eec" else np.diff(edges)
    return edges, d[f"{sel}_{obs}{k}_h"] / nj / w, np.sqrt(d[f"{sel}_{obs}{k}_h2"]) / nj / w


print("\n# mean jet mass and girth (R=0.4, charged constituents)")
for skey, (pre, seeds, pppat, title, sel, jbins) in SAMPLES.items():
    pp = D[skey]["pp"]
    nk = len(jbins)
    fig, ax = plt.subplots(len(OBS), nk, figsize=(4.2 * nk, 2.8 * len(OBS)), squeeze=False)
    for i, (obs, oname, logx) in enumerate(OBS):
        for k in range(nk):
            a = ax[i, k]
            edges, yp, ep = per_jet(pp, sel, obs, k)
            x = np.sqrt(edges[:-1] * edges[1:]) if logx else 0.5 * (edges[:-1] + edges[1:])
            res = {}
            for m, (lab, c, mk) in MODELS.items():
                _, y, e = per_jet(D[skey][m][0], sel, obs, k)
                r, er = ratio(y, e, yp, ep)
                res[m] = (r, er)
                ok = np.isfinite(r) & (yp > 0.02 * np.nanmax(yp))
                a.errorbar(x[ok], r[ok], yerr=er[ok], fmt=mk, ms=3, color=c, label=lab)
            c2, n = chi2_diff(*res["emdv_l256"], *res["htlv_k16"])
            SIG.append((skey, oname, jbins[k], c2, n))
            a.axhline(1, color="gray", lw=0.7); a.set_ylim(0, 2.5)
            if logx:
                a.set_xscale("log")
            a.set_title(f"{oname}, {jbins[k]} GeV   EMD-HTL $\\chi^2$/n={c2:.1f}/{n}", fontsize=8)
            if k == 0:
                a.set_ylabel("AA/pp (per jet)", fontsize=8)
    ax[0, 0].legend(fontsize=7)
    fig.suptitle(f"{title}: R=0.4 jet substructure"); fig.tight_layout()
    fig.savefig(os.path.join(HERE, f"fig_sub_{skey}.png"), dpi=100)
    for obs in ("mass", "ang1"):
        for k in range(nk):
            line = f"  {skey} {obs:5s} {jbins[k]:>7s}"
            for key, d in [("pp", pp)] + [(m, D[skey][m][0]) for m in MODELS]:
                ed = d[f"{sel}_{obs}{k}"]; h = d[f"{sel}_{obs}{k}_h"]
                line += f"  {key}={np.sum(h * 0.5 * (ed[:-1] + ed[1:])) / np.sum(h):.3f}"
            print(line)

# ---------------------------------------------------------------- (5) radius dependence
print("\n# R_AA(0.6)/R_AA(0.2), weighted mean over jet pT")
fig, ax = plt.subplots(2, 4, figsize=(19, 8))
for col, (skey, (pre, seeds, pppat, title, sel, jbins)) in enumerate(SAMPLES.items()):
    pp = D[skey]["pp"]
    rr = {}
    for m, (lab, c, mk) in MODELS.items():
        d = D[skey][m][0]
        vals = {}
        for R in ("R02", "R04", "R06"):
            vals[R] = raa(d, pp, f"{sel}_spec{R}")
        e = vals["R02"][0]; x = 0.5 * (e[:-1] + e[1:])
        for R, ls in (("R02", ":"), ("R04", "-"), ("R06", "--")):
            ax[0, col].errorbar(x, vals[R][1], yerr=vals[R][2], fmt=mk, ls=ls, ms=3, color=c,
                                label=f"{lab} R={int(R[1:]) / 10:.1f}" if col == 0 else None)
        q, eq = ratio(vals["R06"][1], vals["R06"][2], vals["R02"][1], vals["R02"][2])
        rr[m] = (q, eq)
        ax[1, col].errorbar(x, q, yerr=eq, fmt=mk, ms=4, color=c)
        ok = np.isfinite(q) & (eq > 0); w = 1 / eq[ok] ** 2
        print(f"  {skey} {m:10s} {np.sum(q[ok] * w) / np.sum(w):.3f} +- {1 / math.sqrt(np.sum(w)):.3f}")
    c2, n = chi2_diff(*rr["emdv_l256"], *rr["htlv_k16"])
    SIG.append((skey, "R_AA(0.6)/R_AA(0.2)", "all", c2, n))
    ax[0, col].set_title(title); ax[1, col].set_title(f"R_AA(0.6)/R_AA(0.2)   EMD-HTL χ²/n={c2:.1f}/{n}", fontsize=9)
    for a in ax[:, col]:
        a.set_xscale("log"); a.axhline(1, color="gray", lw=0.7); a.set_xlabel(r"jet $p_T$ [GeV]")
    ax[0, col].set_ylim(0, 1.2); ax[1, col].set_ylim(0.4, 2.0)
ax[0, 0].legend(fontsize=6, ncol=2)
fig.tight_layout(); fig.savefig(os.path.join(HERE, "fig_radius.png"), dpi=110)

# ---------------------------------------------------------------- summary
print("\n# EMD vs HTL: all observables")
print(f"{'sample':9s} {'observable':34s} {'bin':>8s} {'chi2/n':>11s} {'p':>7s}")
for skey, oname, jb, c2, n in SIG:
    p = pvalue(c2, n)
    print(f"{skey:9s} {oname:34s} {jb:>8s} {c2:6.1f}/{n:<4d} {p:7.3f}{'  <-- p<0.01' if p < 0.01 else ''}")
ps = [pvalue(c, n) for _, _, _, c, n in SIG if n > 0]
print(f"\n{len(ps)} comparisons; min p = {min(ps):.4f}; expected number with p<0.05 by chance: {0.05 * len(ps):.1f}, observed: {sum(p < 0.05 for p in ps)}")
