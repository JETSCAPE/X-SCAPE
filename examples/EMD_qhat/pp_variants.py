#!/usr/bin/env python3
"""Why do the pp jets carry too much soft, wide-angle energy?

pp 5.02 TeV, PYTHIA 8 hard process + MATTER vacuum shower (Q0 = 1 GeV), with
  ppx             colorless hadronization, PYTHIA MPI on   (production baseline)
  ppv_nompi       colorless, MPI off
  ppv_colored     colored hadronization, MPI on
  ppv_colored_nompi  colored, MPI off
  ppv_hybrid      hybrid hadronization (recombination + strings), MPI on
The colored samples are made with the fixed shared build (2026-10-26). Before
the fix, MATTER left max_color uninitialized after photon radiation, colour
tags overflowed, and ColoredHadronization (which ignored a failed
pythia.next()) wrote unhadronized partons for every later event; that earlier
sample needed 2 x 64 short runs and jetana_data's skipping of parton events.
(PythiaGun switches MPI on and passes all status-62 partons, MPI included, to
MATTER; "PartonLevel:MPI = off" in LinesToRead overrides it.)  Compared with
the pp data used in step 2, with the measurements' UE treatment applied to the
model (data_compare.py).

Writes pp_variants.txt and fig_pp_variants.png.

    python3 pp_variants.py shower

compares instead the parton shower (colorless hadronization, MPI on):
  ppx             MATTER vacuum shower, Q0 = 1 GeV, vir_factor = 0.25 (production)
  ppm_pyfsr       PYTHIA's own final-state shower (FSR_on = 1); MATTER passes the
                  partons through (vir_factor = 2e-6 keeps their virtuality below QS)
  ppm_q02         MATTER stops at Q0 = 2 GeV
  ppm_vir0p1      MATTER starting virtuality vir_factor = 0.1
and writes pp_shower.txt and fig_pp_shower.png.

    python3 pp_variants.py pythia

compares the JETSCAPE production pp with standalone PYTHIA 8 (pythia_pp.cc:
default Monash tune, full event, the same pTHat bias and decays), with MPI on
and off, and both with the pion mass for charged particles as in ALICE's
track-based jets (spec_data_pimass/, jetana_data with JETANA_PIONMASS=1);
writes pp_pythia.txt and fig_pp_pythia.png.  In the other modes the ALICE
rows always use the pion mass.
"""
import glob
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from data_compare import HERE, SPEC, SPEC_PI, chi2, draw_data, hep_table, per_jet, shape, shape_sub, unit
from spec_tools import combine

VARIANTS = {"ppx": ("colorless, MPI on (production)", "#1F5FA8", "-"),
            "ppv_nompi": ("colorless, MPI off", "#1F5FA8", ":"),
            "ppv_colored": ("colored, MPI on", "#C2410C", "-"),
            "ppv_colored_nompi": ("colored, MPI off", "#C2410C", ":"),
            "ppv_hybrid": ("hybrid, MPI on", "#15803D", "-")}
FILES = {"ppx": ["ppx_s%d" % i for i in range(1, 9)],
         "ppv_nompi": ["ppv_nompi_s%d" % i for i in range(1, 9)],
         "ppv_colored": ["ppv_colored_s%d" % i for i in range(1, 9)],
         "ppv_colored_nompi": ["ppv_colored_nompi_s%d" % i for i in range(1, 9)],
         "ppv_hybrid": ["ppv_hybrid_s%d" % i for i in range(1, 9)]}


def counts(paths):
    """Events used and skipped (events with unhadronized partons)."""
    used = skipped = 0
    for f in paths:
        t = open(f).readline().split()
        used += int(t[3])
        skipped += int(t[5]) if len(t) > 5 else 0
    return used, skipped


def step(ax, e, y, v):
    lab, col, ls = VARIANTS[v]
    ax.step(e, np.append(y, y[-1]), where="post", color=col, ls=ls, lw=1.5, label=lab)


OUT = ("fig_pp_variants.png", "pp_variants.txt")
SHOWER_VARIANTS = {"ppx": ("MATTER, Q0 = 1, vir_factor 0.25 (production)", "#1F5FA8", "-"),
                   "ppm_pyfsr": ("PYTHIA shower (FSR_on = 1)", "#C2410C", "-"),
                   "ppm_q02": ("MATTER, Q0 = 2 GeV", "#15803D", "-"),
                   "ppm_vir0p1": ("MATTER, vir_factor 0.1", "#7C3AED", "--")}
for _v in SHOWER_VARIANTS:
    FILES.setdefault(_v, ["%s_s%d" % (_v, i) for i in range(1, 9)])
PYTHIA_VARIANTS = {"ppx": ("JETSCAPE pp (production)", "#1F5FA8", "-"),
                   "pyth": ("PYTHIA 8 standalone, MPI on", "#C2410C", "-"),
                   "pyth_nompi": ("PYTHIA 8 standalone, MPI off", "#C2410C", ":"),
                   "ppx_pimass": ("JETSCAPE pp, pion mass", "#1F5FA8", "--"),
                   "pyth_pimass": ("PYTHIA 8, MPI on, pion mass", "#C2410C", "--")}
FILES.update({"pyth": ["pyth_s%d" % i for i in range(1, 9)],
              "pyth_nompi": ["pyth_nompi_s%d" % i for i in range(1, 9)],
              "ppx_pimass": ["../spec_data_pimass/ppx_s%d" % i for i in range(1, 9)],
              "pyth_pimass": ["../spec_data_pimass/pyth_s%d" % i for i in range(1, 9)]})


ALICE_PION_MASS = True  # ALICE rows from spec_data_pimass, as ALICE assigns tracks the pion mass


def main():
    D, DA = {}, {}
    for v in VARIANTS:
        fl = [os.path.join(SPEC, f"{n}.spec") for n in FILES[v]]
        fl = [f for f in fl if os.path.exists(f)]
        if fl:
            D[v] = combine(fl)
            D[v]["nfiles"] = len(fl)
            D[v]["nused"], D[v]["nskip"] = counts(fl)
            fa = [os.path.join(SPEC_PI, os.path.basename(f)) for f in fl] if ALICE_PION_MASS else fl
            DA[v] = combine(fa)
    V = list(D)
    rows = []  # (observable, {variant: (chi2, n)})

    t11 = hep_table("ins1658057", "Table 11")
    bpp = unit([b for b in t11 if "P P" in b["quals"].get("RE", "")][0])
    fig, axs = plt.subplots(3, 4, figsize=(14, 10))
    # CMS rho, UE subtracted and not
    for col, (sub, lab) in enumerate([(True, "CMS $\\rho$ pp, UE subtracted (sideband)"), (False, "CMS $\\rho$ pp, no UE subtraction")]):
        ax = axs[0, col]
        res = {}
        for v in V:
            e, y, ey = (shape_sub(D[v], "cms_rhoall", "cms_rhosb") if sub else shape(D[v], "cms_rhoall"))
            res[v] = chi2(e, y, ey, bpp)
            step(ax, e, y, v)
        rows.append((lab.replace("$\\rho$", "rho"), res))
        draw_data(ax, bpp, "CMS")
        ax.set_yscale("log")
        ax.set_xlabel(r"$\Delta r$")
        ax.set_title(lab, fontsize=9)
    axs[0, 0].legend(fontsize=7, frameon=False)
    # ATLAS D(z) pp, two pT bins, UE subtracted
    for col, k in enumerate([0, 2]):
        dz = hep_table("ins1673177", f"Table {4 * k}")
        b = [x for x in dz if x["quals"].get("Centrality") == "pp"][0]
        ax = axs[0, 2 + col]
        res = {}
        for v in V:
            e, y, ey = per_jet(D[v], f"atl_ffz{k}", "atl_ffnjet", k, f"atl_ffzbkg{k}")
            res[v] = chi2(e, y, ey, b)
            yd, _ = chi2_ratio(e, y, b)
            step(ax, e, yd, v)
        rows.append((f"ATLAS D(z) pp {['126-158', '158-200', '200-251'][k]}", res))
        ax.axhline(1, color="0.6", lw=0.6)
        x = 0.5 * (b["lo"] + b["hi"])
        ax.fill_between(x, 1 - np.hypot(b["stat"], b["sys"]) / b["y"], 1 + np.hypot(b["stat"], b["sys"]) / b["y"],
                        color="0.8", step="mid", label="ATLAS unc.")
        ax.set_xscale("log")
        ax.set_ylim(0.6, 1.6)
        ax.set_xlabel("$z$")
        ax.set_title(f"D(z) pp model/ATLAS, {['126-158', '158-200', '200-251'][k]} GeV", fontsize=9)
    # ALICE pp mass ungroomed / groomed
    for g in range(2):
        for k in range(3):
            b = unit(hep_table("ins2845788", f"Table {33 + 2 * k + g}")[0])
            ax = axs[1 + g, k]
            res = {}
            for v in V:
                e, y, ey = shape(DA[v], (f"ali_mg{k}" if g else f"ali_m{k}"))
                res[v] = chi2(e, y, ey, b)
                step(ax, e, y, v)
            lab = f"ALICE pp {'groomed ' if g else ''}mass {['40-60', '60-80', '80-100'][k]}"
            rows.append((lab, res))
            draw_data(ax, b, "ALICE")
            ax.set_title(lab, fontsize=9)
            ax.set_xlabel("$m$ (GeV)")
    # ATLAS x_J pp, all pT1 bins summed in chi2
    xj = hep_table("ins2075431", "pp self norm")
    res = {v: [0., 0] for v in V}
    for k in range(11):
        b = unit(xj[k])
        for v in V:
            e, y, ey = shape(D[v], f"atl_xj{k}")
            c, n = chi2(e, y, ey, b)
            res[v][0] += c
            res[v][1] += n
    rows.append(("ATLAS xJ pp (11 pT1 bins)", {v: tuple(res[v]) for v in V}))
    # mean jet-mass shift table
    ax = axs[1, 3]
    means = {}
    for v in V:
        mm = means[v] = []
        for k in range(3):
            e, y, _ = shape(DA[v], f"ali_m{k}")
            mm.append(np.sum(0.5 * (e[:-1] + e[1:]) * y * np.diff(e)))
        ax.plot([50, 70, 90], mm, marker="o", color=VARIANTS[v][1], ls=VARIANTS[v][2], label=VARIANTS[v][0])
    md = []
    for k in range(3):
        b = unit(hep_table("ins2845788", f"Table {33 + 2 * k}")[0])
        md.append(np.sum(0.5 * (b["lo"] + b["hi"]) * b["y"] * (b["hi"] - b["lo"])))
    ax.plot([50, 70, 90], md, "ko", label="ALICE")
    ax.set_xlabel(r"$p_T^{ch\,jet}$ (GeV)")
    ax.set_ylabel(r"$\langle m \rangle$ (GeV)")
    ax.set_title("ALICE pp mean jet mass", fontsize=9)
    ax.legend(fontsize=6, frameon=False)
    axs[2, 3].axis("off")
    fig.tight_layout()
    fig.savefig(os.path.join(HERE, OUT[0]), dpi=140)
    plt.close(fig)

    out = ["pp variants vs data: chi2 / n_bins",
           "files (events used / skipped for unhadronized partons): "
           + ", ".join(f"{v} {D[v]['nfiles']} ({D[v]['nused']}/{D[v]['nskip']})" for v in V), "",
           f"  {'observable':<42s}" + "".join(f"{v:>20s}" for v in V)]
    for lab, res in rows:
        out.append(f"  {lab:<42s}" + "".join(f"{res[v][0]:12.1f} /{res[v][1]:<4d}  " for v in V))
    out += ["", "ALICE pp <m> (GeV) at 40-60, 60-80, 80-100 (HEPData bins):",
            f"  {'ALICE':<20s}" + " ".join(f"{x:6.2f}" for x in md)]
    out += [f"  {v:<20s}" + " ".join(f"{x:6.2f}" for x in means[v])
            + "   (model - data: " + " ".join(f"{x - y:+.2f}" for x, y in zip(means[v], md)) + ")" for v in V]
    open(os.path.join(HERE, OUT[1]), "w").write("\n".join(out) + "\n")
    print("\n".join(out))


def chi2_ratio(e, y, b):
    """Model / data on the data bins (NaN where the bins do not match)."""
    r = np.full(len(e) - 1, np.nan)
    for i, (lo, hi) in enumerate(zip(b["lo"], b["hi"])):
        j = np.where(np.abs(e[:-1] - lo) < 1e-3 * hi)[0]
        if len(j):
            r[j[0]] = y[j[0]] / b["y"][i]
    return r, None


if __name__ == "__main__":
    import sys
    if sys.argv[1:] == ["shower"]:
        VARIANTS.clear()
        VARIANTS.update(SHOWER_VARIANTS)
        OUT = ("fig_pp_shower.png", "pp_shower.txt")
    elif sys.argv[1:] == ["pythia"]:
        ALICE_PION_MASS = False  # this mode shows both conventions explicitly
        VARIANTS.clear()
        VARIANTS.update(PYTHIA_VARIANTS)
        OUT = ("fig_pp_pythia.png", "pp_pythia.txt")
    main()
