#!/usr/bin/env python3
"""Option 1 ("basics"): stock JETSCAPE MATTER+LBT with only the radiative
qhat(T) swapped for the EMD holographic one, compared with the default tune.

Models (RHIC and LHC 0-10%, same hydro seeds as the tune production):
  tune          MATTER (Type 5) + LBT (tables)               - default
  lbtemd_l2     MATTER (Type 9) + LBT (EMD radiative qhat), lambda_t = 2
                (matches the tune's 0-10% LHC jet suppression; RHIC has
                80 hydro events, seeds 56-95 and 171-210)
  lbtemd_l2p6   MATTER (Type 9) + LBT (EMD radiative qhat), lambda_t = 2.6
                (same 0-10% LHC jet suppression as the tune)
  lbtemd_l8     same with lambda_t = 8 (fits the ATLAS 0-10% jet R_AA)

1. Paired per-seed differences ln(EMD/tune) of yields and substructure means,
   with the standard error over seeds (the runs of a seed share the hydro event
   and the Pythia hard scatterings).
2. Figure option1.png: jet R_AA with data, R_AA(0.6)/R_AA(0.2), and the EEC
   and girth AA/pp per jet for one jet-pT bin, RHIC and LHC.
"""
import csv
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
DATA = os.path.join(HERE, "data")

MODELS = {"tune": ("MATTER+LBT default", "#1f5fa8", "o"),
          "lbtemd_l2": (r"MATTER+LBT, EMD $\hat q$, $\lambda_t$=2", "#15803d", "D"),
          "lbtemd_l2p6": (r"MATTER+LBT, EMD $\hat q$, $\lambda_t$=2.6", "#c2410c", "s"),
          "lbtemd_l8": (r"MATTER+LBT, EMD $\hat q$, $\lambda_t$=8", "#7e22ce", "^")}
SAMPLES = {  # key: (prefix, seeds, pp pattern, title, selection, jet spectrum, jet-pT bins, bin shown)
    "rhic0010": ("rhic", list(range(56, 96)) + list(range(171, 331)), "rhic_pp_s*", "Au-Au 200 GeV 0-10%", "rs", "sphenix",
                 ["15-20", "20-30", "30-60"], 1),
    "rhic3050": ("rhic_c3050", range(131, 171), "rhic_pp_s*", "Au-Au 200 GeV 30-50%", "rs", "sphenix",
                 ["15-20", "20-30", "30-60"], 1),
    "lhc0010": ("lhc", list(range(16, 46)) + list(range(331, 361)), "pp_s*", "Pb-Pb 5.02 TeV 0-10%", "ls", "fulljet",
                ["40-60", "60-80", "80-120", "120-200"], 2),
    "lhc3050": ("lhc_c3050", range(116, 151), "pp_s*", "Pb-Pb 5.02 TeV 30-50%", "ls", "fulljet",
                ["40-60", "60-80", "80-120", "120-200"], 2),
}
RANGES = {"rs": [(15, 30), (30, 60)], "ls": [(40, 100), (100, 316)]}


# ---------------------------------------------------------------- helpers
def yield_in(d, obs, lo, hi):
    e = d[obs]
    m = (e[:-1] >= lo - 1e-9) & (e[1:] <= hi + 1e-9)
    return np.sum(d[obs + "_h"][m]) / d["sumw"]


def mean_of(d, sel, obs, k):
    e = d[f"{sel}_{obs}{k}"]
    h = d[f"{sel}_{obs}{k}_h"]
    return np.sum(h * 0.5 * (e[:-1] + e[1:])) / np.sum(h)


def eec_narrow_frac(d, sel, k, cut=0.1):
    e = d[f"{sel}_eec{k}"]
    h = d[f"{sel}_eec{k}_h"]
    return np.sum(h[e[1:] <= cut + 1e-9]) / np.sum(h)


def observables(d, sel, jet):
    out = {}
    for lo, hi in RANGES[sel]:
        out[f"jet yield {lo}-{hi}"] = yield_in(d, jet, lo, hi)
    hlo = 10 if sel == "ls" else 5
    out[f"hadron yield {hlo}-{4 * hlo}"] = yield_in(d, "hadron", hlo, 4 * hlo)
    lo, hi = RANGES[sel][0][0], RANGES[sel][-1][1]
    out["yield R=0.6 / R=0.2"] = yield_in(d, f"{sel}_specR06", lo, hi) / yield_in(d, f"{sel}_specR02", lo, hi)
    nk = len(d[f"{sel}_njet"]) - 1
    for k in range(nk):
        lab = f"{d[f'{sel}_njet'][k]:.0f}-{d[f'{sel}_njet'][k + 1]:.0f}"
        out[f"<mass> {lab}"] = mean_of(d, sel, "mass", k)
        out[f"<girth> {lab}"] = mean_of(d, sel, "ang1", k)
        out[f"<lambda2> {lab}"] = mean_of(d, sel, "ang2", k)
        out[f"EEC(R_L<0.1) frac {lab}"] = eec_narrow_frac(d, sel, k)
    return out


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


def per_jet(d, sel, obs, k):
    edges = d[f"{sel}_{obs}{k}"]
    nj = d[f"{sel}_njet_h"][k]
    w = np.diff(np.log(edges)) if obs == "eec" else np.diff(edges)
    return edges, d[f"{sel}_{obs}{k}_h"] / nj / w, np.sqrt(d[f"{sel}_{obs}{k}_h2"]) / nj / w


def hepdata(name, xcol, ycol, stat, sysp, sysm):
    rows = [r for r in csv.reader(l for l in open(os.path.join(DATA, name)) if not l.startswith("#")) if r][1:]
    a = np.array([[float(r[i]) for i in (xcol, ycol, stat, sysp, sysm)] for r in rows])
    a[:, 2:] = np.abs(a[:, 2:])
    return a


# ---------------------------------------------------------------- 1. paired
PAIRED = {}
print(f"{'sample':9s} {'model':12s} {'observable':26s} {'EMD/tune - 1':>15s} {'signif.':>8s}  n")
for skey, (pre, seeds, _pp, _t, sel, jet, _jb, _k) in SAMPLES.items():
    for m in list(MODELS)[1:]:
        diffs = {}
        for s in seeds:
            fe = os.path.join(SPEC, f"{pre}_{m}_h{s}.spec")
            ft = os.path.join(SPEC, f"{pre}_tune_h{s}.spec")
            if not (os.path.exists(fe) and os.path.exists(ft)):
                continue
            oe = observables(read_spec(fe), sel, jet)
            ot = observables(read_spec(ft), sel, jet)
            for k in oe:
                if oe[k] > 0 and ot[k] > 0 and np.isfinite(oe[k]) and np.isfinite(ot[k]):
                    diffs.setdefault(k, []).append(math.log(oe[k] / ot[k]))
        for k, v in diffs.items():
            v = np.array(v)
            mu, se = v.mean(), v.std(ddof=1) / math.sqrt(len(v))
            PAIRED[(m, skey, k)] = (mu, se)
            flag = "  <--" if abs(mu) > 3 * se else ""
            print(f"{skey:9s} {m:12s} {k:26s} {100 * mu:+7.2f} +- {100 * se:4.2f}% {mu / se:+7.1f}σ  {len(v)}{flag}")
        print()

# ------------------------------------------------- 1b. RHIC vs LHC double ratio
# [ln(EMD/tune)]_RHIC - [ln(EMD/tune)]_LHC: the T-shape of qhat at fixed LHC
# normalisation (independent seeds, so the errors add in quadrature)
JETOBS = {"rs": ("jet yield 15-30", "jet yield 30-60", "hadron yield 5-20"),
          "ls": ("jet yield 40-100", "jet yield 100-316", "hadron yield 10-40")}
print("# double ratios: [ln(EMD/tune)](A) - [ln(EMD/tune)](B)")
for m, (lab, _c, _mk) in list(MODELS.items())[1:]:
    for sa, sb in (("rhic0010", "lhc0010"), ("rhic3050", "lhc0010"), ("rhic3050", "rhic0010"),
                   ("lhc3050", "lhc0010"), ("rhic3050", "lhc3050")):
        oa, ob = JETOBS[SAMPLES[sa][4]], JETOBS[SAMPLES[sb][4]]
        for ra, rb in ((oa[0], ob[0]), (oa[0], ob[1]), (oa[1], ob[0]), (oa[2], ob[2])):
            if (m, sa, ra) in PAIRED and (m, sb, rb) in PAIRED:
                (a, ea), (b, eb) = PAIRED[(m, sa, ra)], PAIRED[(m, sb, rb)]
                d, ed = a - b, math.hypot(ea, eb)
                print(f"  {m:12s} {sa} {ra:18s} vs {sb} {rb:18s} {100 * d:+7.2f} +- {100 * ed:4.2f}%  {d / ed:+5.1f}σ")
print()

# ---------------------------------------------------------------- 2. figure
fig, ax = plt.subplots(len(SAMPLES), 4, figsize=(19, 4.25 * len(SAMPLES)))
for row, (skey, (pre, seeds, pppat, title, sel, jet, jbins, kshow)) in enumerate(SAMPLES.items()):
    pp = combine(sorted(glob.glob(os.path.join(SPEC, pppat + ".spec"))))
    D = {}
    for m in MODELS:
        fl = [os.path.join(SPEC, f"{pre}_{m}_h{s}.spec") for s in seeds]
        fl = [f for f in fl if os.path.exists(f)]
        print(f"{skey} {m}: {len(fl)} hydro events")
        if fl:
            D[m] = combine(fl)
    # (a) jet R_AA
    a = ax[row, 0]
    if skey.startswith("lhc"):
        d = hepdata("ATLAS_jetRAA_R04_0-10_5020_ins1673184_t19.csv" if skey == "lhc0010" else
                    "ATLAS_jetRAA_R04_30-50_5020_ins1673184_t22t23avg.csv", 0, 3, 6, 4, 5)
        a.errorbar(d[:, 0], d[:, 1], yerr=np.hypot(d[:, 2], 0.5 * (d[:, 3] + d[:, 4])), fmt="k*", ms=6,
                   label="ATLAS R=0.4, |y|<2.8")
    for m, (lab, c, mk) in [(k, v) for k, v in MODELS.items() if k in D]:
        e, r, er = raa(D[m], pp, jet)
        a.errorbar(0.5 * (e[:-1] + e[1:]), r, yerr=er, fmt=mk, ms=4, color=c, label=lab)
    a.set_xscale("log"); a.set_ylim(0, 1.2); a.axhline(1, color="gray", lw=0.7)
    a.set_title(f"{title}: R=0.4 jet $R_{{AA}}$"); a.set_xlabel(r"jet $p_T$ [GeV]"); a.legend(fontsize=7)
    # (b) R_AA(0.6)/R_AA(0.2)
    a = ax[row, 1]
    for m, (lab, c, mk) in [(k, v) for k, v in MODELS.items() if k in D]:
        e, r6, e6 = raa(D[m], pp, f"{sel}_specR06")
        _, r2, e2 = raa(D[m], pp, f"{sel}_specR02")
        q, eq = ratio(r6, e6, r2, e2)
        a.errorbar(0.5 * (e[:-1] + e[1:]), q, yerr=eq, fmt=mk, ms=4, color=c)
    a.set_xscale("log"); a.set_ylim(0.4, 2.0); a.axhline(1, color="gray", lw=0.7)
    a.set_title(r"$R_{AA}(R=0.6)/R_{AA}(R=0.2)$"); a.set_xlabel(r"jet $p_T$ [GeV]")
    # (c, d) EEC and girth, AA/pp per jet
    for col, (obs, name, logx) in ((2, ("eec", "EEC", True)), (3, ("ang1", "girth", False))):
        a = ax[row, col]
        edges, yp, ep = per_jet(pp, sel, obs, kshow)
        x = np.sqrt(edges[:-1] * edges[1:]) if logx else 0.5 * (edges[:-1] + edges[1:])
        for m, (lab, c, mk) in [(k, v) for k, v in MODELS.items() if k in D]:
            _, y, e = per_jet(D[m], sel, obs, kshow)
            r, er = ratio(y, e, yp, ep)
            ok = np.isfinite(r) & (yp > 0.02 * np.nanmax(yp))
            a.errorbar(x[ok], r[ok], yerr=er[ok], fmt=mk, ms=3, color=c)
        if logx:
            a.set_xscale("log"); a.set_xlabel(r"$R_L$")
        else:
            a.set_xlabel("girth")
        a.set_ylim(0, 2.5); a.axhline(1, color="gray", lw=0.7)
        a.set_title(f"{name} AA/pp per jet, {jbins[kshow]} GeV")
fig.tight_layout()
fig.savefig(os.path.join(HERE, "option1.png"), dpi=110)
print("wrote option1.png")
