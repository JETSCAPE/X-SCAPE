#!/usr/bin/env python3
"""Step 2: default MATTER+LBT and EMD qhat (lambda_t = 2.4) vs LHC jet data.

Models from the lambda_t = 2.4 reruns (tasks_x.txt), analysed from the kept
event files with jetana_data (spec_data/), binned as the HEPData records:
  ATLAS x_J        ins2075431  self-normalised, pp and 0-10%
  ATLAS D(z)       ins1673177  pp D(z); R_D(z) 0-10% (and 30-40% vs our 30-50%)
  CMS rho(dr)      ins1658057  pp, 0-10%, 30-50%; PbPb/pp
  ALICE m, lambda  ins2845788  R = 0.2 charged jets; pp mass; Pb-Pb mass and
                               angularities (0-10% assumed, central data)
Shape observables are normalised as the data (unit area over the data bins,
or per jet for D(z)).  The underlying event is subtracted from the model as in
the measurements: CMS rho with the 1.5 < |deta| < 2.5 sideband (pp and PbPb),
ATLAS D(z) with perpendicular cones; ALICE subtracts nothing in pp, and the
Pb-Pb background is the thermal medium, absent from the model events.
Groomed ALICE observables use Soft Drop with the holes (signed reclustering).
chi2 = (m - d)^T C^-1 (m - d) with C = data stat + sys + model stat.  For ALICE
the signed "trkeff" and "generator" sources are fully correlated across bins
(as stated in the record) and the rest uncorrelated; for the other records
all uncertainties are taken as uncorrelated.  The ALICE totals are also given
with all uncertainties uncorrelated (diag).

Writes data_compare.txt and fig_data_{xj,ff,cms,alice_mass,alice_ang}.png.
"""
import glob
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from spec_tools import combine

HERE = os.path.dirname(os.path.abspath(__file__))
SPEC = os.path.join(HERE, "spec_data")
HEP = next((p for p in (os.path.join(HERE, "data", "hepdata"), os.path.join(HERE, "data", "hepdata_raw"))
            if os.path.isdir(p)), os.path.join(HERE, "data", "hepdata"))
SAMPLES = {"lhc0010": ("lhcx", list(range(16, 46)) + list(range(331, 361))),
           "lhc3050": ("lhcx_c3050", list(range(116, 151)))}
MODELS = {"tune": ("default MATTER+LBT", "#1F5FA8", "--"),
          "lbtemd_l2p4": (r"EMD $\hat q$, $\lambda_t=2.4$", "#C2410C", "-")}
XJPT = [100, 112, 126, 141, 158, 178, 200, 224, 251, 282, 316, 398]
FFPT = ["126-158", "158-200", "200-251", "251-316", "316-398"]
ALIPT = ["40-60", "60-80", "80-100", "100-150"]
ALPHA = [("1", "1"), ("15", "1.5"), ("2", "2"), ("3", "3")]


# ---------------------------------------------------------------- HEPData
def hep_table(ins, name):
    """Blocks of a HEPData csv table: dicts with quals, head, lo, hi, y, stat, sys."""
    for f in glob.glob(os.path.join(HEP, ins, "**", "*.csv"), recursive=True):
        lines = open(f, errors="ignore").read().splitlines()
        if not any(l.startswith("#: name:") and l[8:].strip() == name for l in lines[:5]):
            continue
        blocks, quals, cur = [], {}, None
        for l in lines:
            if l.startswith("#:"):
                if ",,," in l:
                    k, v = l[2:].split(",,,", 1)
                    quals[k.strip()] = v.strip()
                continue
            if not l.strip():
                cur = None
                continue
            c = l.split(",")
            try:
                v = [float(x) for x in c]
            except ValueError:  # column header (may contain quoted commas)
                cur = {"quals": dict(quals), "head": l, "rows": []}
                blocks.append(cur)
                hdr = []
                for part in __import__("csv").reader([l]):
                    hdr = part
                cur["cols"] = hdr
                continue
            cur["rows"].append(v)
        out = []
        for b in blocks:
            if not b["rows"]:  # description lines continued outside "#:"
                continue
            r = np.array(b["rows"])
            stat2 = np.zeros(len(r))
            sys2 = np.zeros(len(r))
            corr = []  # signed, bin-to-bin correlated sources
            for i in range(4, r.shape[1], 2):
                e2 = np.maximum(np.abs(r[:, i]), np.abs(r[:, i + 1])) ** 2
                name = b["cols"][i].lower()
                if "stat" in name:
                    stat2 += e2
                else:
                    sys2 += e2
                    if "trkeff" in name or "generator" in name:
                        corr.append(r[:, i])
            out.append({"quals": b["quals"], "head": b["head"], "lo": r[:, 1], "hi": r[:, 2],
                        "y": r[:, 3], "stat": np.sqrt(stat2), "sys": np.sqrt(sys2),
                        "corr": corr})
        return out
    raise KeyError(f"{ins} {name}")


def unit(b):
    """Renormalise a data block to unit area (shape comparison)."""
    a = np.sum(b["y"] * (b["hi"] - b["lo"]))
    return {**b, "y": b["y"] / a, "stat": b["stat"] / a, "sys": b["sys"] / a,
            "corr": [c / a for c in b["corr"]]}


# ---------------------------------------------------------------- model
def load():
    D = {"pp": {}}
    for m in MODELS:
        D["pp"][m] = combine(sorted(glob.glob(os.path.join(SPEC, "ppx_s*.spec"))))
    for s, (pre, seeds) in SAMPLES.items():
        D[s] = {}
        for m in MODELS:
            fl = [os.path.join(SPEC, f"{pre}_{m}_h{i}.spec") for i in seeds]
            fl = [f for f in fl if os.path.exists(f)]
            D[s][m] = combine(fl)
            D[s][m]["nfiles"] = len(fl)
    return D


def shape(d, key):
    """Unit-area distribution with errors."""
    e, h, h2 = d[key], d[key + "_h"], d[key + "_h2"]
    w = np.diff(e)
    a = np.sum(h)
    return e, h / a / w, np.sqrt(h2) / a / w


def per_jet(d, key, nkey, k, bkg=None):
    """Per-jet density; bkg: underlying-event histogram to subtract."""
    e, h, h2 = d[key], d[key + "_h"], d[key + "_h2"]
    if bkg:
        h, h2 = h - d[bkg + "_h"], h2 + d[bkg + "_h2"]
    n = d[nkey + "_h"][k]
    w = np.diff(e)
    return e, h / n / w, np.sqrt(h2) / n / w


def shape_sub(d, key, bkg):
    """Unit-area distribution after subtracting the background histogram."""
    e = d[key]
    h, h2 = d[key + "_h"] - d[bkg + "_h"], d[key + "_h2"] + d[bkg + "_h2"]
    w = np.diff(e)
    a = np.sum(h)
    return e, h / a / w, np.sqrt(h2) / a / w


def ratio(a, ea, b, eb):
    with np.errstate(divide="ignore", invalid="ignore"):
        r = a / b
        return r, np.abs(r) * np.sqrt((ea / a) ** 2 + (eb / b) ** 2)


def chi2(e, m, em, b, mask=None, diag=False):
    """Model bins matched to the data bins by their edges (NaN where none)."""
    mm, ee = np.full(len(b["y"]), np.nan), np.full(len(b["y"]), np.nan)
    for i, (lo, hi) in enumerate(zip(b["lo"], b["hi"])):
        j = np.where((np.abs(e[:-1] - lo) < 1e-3 * max(hi, 1e-3)) & (np.abs(e[1:] - hi) < 1e-3 * max(hi, 1e-3)))[0]
        if len(j):
            mm[i], ee[i] = m[j[0]], em[j[0]]
    m, em = mm, ee
    s2 = b["stat"] ** 2 + b["sys"] ** 2 + em ** 2
    ok = np.isfinite(m) & (s2 > 0) & np.isfinite(em)
    if mask is not None:
        ok &= mask
    if diag or not b.get("corr"):
        return float(np.sum((m[ok] - b["y"][ok]) ** 2 / s2[ok])), int(ok.sum())
    cv = [c[ok] for c in b["corr"]]
    C = np.diag(s2[ok] - sum(c ** 2 for c in cv)) + sum(np.outer(c, c) for c in cv)
    r = m[ok] - b["y"][ok]
    return float(r @ np.linalg.solve(C, r)), int(ok.sum())


# ---------------------------------------------------------------- plotting
def draw_data(ax, b, label="data"):
    x = 0.5 * (b["lo"] + b["hi"])
    tot = np.hypot(b["stat"], b["sys"])
    ax.errorbar(x, b["y"], yerr=tot, xerr=0.5 * (b["hi"] - b["lo"]), fmt="o", ms=3.5,
                color="k", elinewidth=0.8, capsize=0, label=label, zorder=5)


def draw_model(ax, e, y, ey, m, label=True):
    lab, col, ls = MODELS[m]
    yy = np.append(y, y[-1])
    ax.step(e, yy, where="post", color=col, ls=ls, lw=1.6, label=lab if label else None)
    xc = 0.5 * (e[:-1] + e[1:])
    ax.fill_between(xc, y - ey, y + ey, color=col, alpha=0.18, step="mid", lw=0)


# ---------------------------------------------------------------- main
def main():
    D = load()
    report = []
    tot = {}  # group -> model -> [chi2, n]

    def add(group, table, res, resd=None):
        if resd is not None:
            add(group + " (diag)", table + " diag", resd)
        line = f"  {table:<44s}"
        for m in MODELS:
            c, n = res[m]
            line += f"  {m}: {c:7.1f}/{n:<3d}"
            tot.setdefault(group, {}).setdefault(m, [0., 0])
            tot[group][m][0] += c
            tot[group][m][1] += n
        report.append(line)

    # ---------------- ATLAS x_J
    report.append("ATLAS dijet x_J (self-normalised), chi2/n per pT1 bin")
    xj_aa = hep_table("ins2075431", "PbPb self norm cent0")
    xj_pp = hep_table("ins2075431", "pp self norm")
    show = [0, 2, 4, 6, 8]
    fig, axs = plt.subplots(2, len(show), figsize=(3.2 * len(show), 6.4), sharey="row")
    for row, (lab, blocks, d) in enumerate([("pp", xj_pp, D["pp"]), ("0-10%", xj_aa, D["lhc0010"])]):
        for k in range(len(XJPT) - 1):
            b = unit(blocks[k])
            res = {}
            for m in MODELS:
                e, y, ey = shape(d[m], f"atl_xj{k}")
                res[m] = chi2(e, y, ey, b)
                if k in show:
                    ax = axs[row, show.index(k)]
                    draw_model(ax, e, y, ey, m, label=(row == 0 and k == 0))
            add(f"ATLAS xJ {lab}", f"{lab} pT1 {XJPT[k]}-{XJPT[k + 1]}", res)
            if k in show:
                ax = axs[row, show.index(k)]
                draw_data(ax, b, "ATLAS" if (row == 0 and k == 0) else None)
                ax.set_xlim(0.3, 1.0)
                ax.set_title(f"{lab}, $p_{{T,1}}$ {XJPT[k]}–{XJPT[k + 1]} GeV", fontsize=9)
                if row == 1:
                    ax.set_xlabel("$x_J$")
        axs[row, 0].set_ylabel(r"$(1/N_{\rm pair})\,dN/dx_J$")
    axs[0, 0].legend(fontsize=7, frameon=False)
    fig.tight_layout()
    fig.savefig(os.path.join(HERE, "fig_data_xj.png"), dpi=150)
    plt.close(fig)

    # ---------------- ATLAS D(z)
    report.append("ATLAS D(z): pp D(z), R_D(z) = PbPb/pp (UE subtracted, perpendicular cones)")
    resn = {}
    fig, axs = plt.subplots(3, len(FFPT), figsize=(3.2 * len(FFPT), 9), sharex=True)
    for k, ptl in enumerate(FFPT):
        dz = hep_table("ins1673177", f"Table {4 * k}")
        rd = hep_table("ins1673177", f"Table {4 * k + 2}")
        bpp = [b for b in dz if b["quals"].get("Centrality") == "pp"][0]
        r00 = [b for b in rd if b["quals"].get("Centrality") == "0-10%"][0]
        r34 = [b for b in rd if b["quals"].get("Centrality") == "30-40%"][0]
        for row, (b, lab) in enumerate([(bpp, "pp D(z)"), (r00, "R_D 0-10%"), (r34, "R_D 30-40% (model 30-50%)")]):
            ax = axs[row, k]
            res = {}
            for m in MODELS:
                e, ypp, epp = per_jet(D["pp"][m], f"atl_ffz{k}", "atl_ffnjet", k, f"atl_ffzbkg{k}")
                if row == 0:
                    y, ey = ypp, epp
                    _, yn, en = per_jet(D["pp"][m], f"atl_ffz{k}", "atl_ffnjet", k)
                    resn = resn if m != list(MODELS)[0] else {}
                    resn[m] = chi2(e, yn, en, b)
                    if m == list(MODELS)[-1]:
                        add("ATLAS D(z) pp D(z), no UE sub", f"pp D(z) no UE sub pT {ptl}", resn)
                    ax.step(e, np.append(yn, yn[-1]), where="post", color="0.5", ls=":", lw=1.2,
                            label="no UE subtraction" if (k == 0 and m == list(MODELS)[0]) else None)
                else:
                    s = "lhc0010" if row == 1 else "lhc3050"
                    _, ya, ea = per_jet(D[s][m], f"atl_ffz{k}", "atl_ffnjet", k, f"atl_ffzbkg{k}")
                    y, ey = ratio(ya, ea, ypp, epp)
                res[m] = chi2(e, y, ey, b)
                draw_model(ax, e, y, ey, m, label=(row == 0 and k == 0))
            add(f"ATLAS D(z) {lab.split(' (')[0]}", f"{lab} pT {ptl}", res)
            draw_data(ax, b, "ATLAS" if (row == 0 and k == 0) else None)
            ax.set_xscale("log")
            if row == 0:
                ax.set_yscale("log")
                ax.set_title(f"$p_T^{{jet}}$ {ptl} GeV", fontsize=9)
            else:
                ax.axhline(1, color="0.6", lw=0.6)
                ax.set_ylim(0.4, 2.2)
            if row == 2:
                ax.set_xlabel("$z$")
            if k == 0:
                ax.set_ylabel(["$D(z)$ pp", "$R_{D(z)}$ 0–10%", "$R_{D(z)}$ 30–40% / 30–50%"][row])
    axs[0, 0].legend(fontsize=7, frameon=False)
    fig.tight_layout()
    fig.savefig(os.path.join(HERE, "fig_data_ff.png"), dpi=150)
    plt.close(fig)

    # ---------------- CMS rho(dr)
    report.append("CMS jet shape rho(dr) (unit area over 0 < dr < 1) and PbPb/pp, "
                  "UE subtracted with the 1.5 < |deta| < 2.5 sideband")
    t11 = hep_table("ins1658057", "Table 11")
    t12 = hep_table("ins1658057", "Table 12")

    def cblock(t, c):
        return [b for b in t if b["quals"].get("CENTRALITY") == c or
                (c == "pp" and "P P" in b["quals"].get("RE", ""))][0]
    def rho(d):
        return shape_sub(d, "cms_rhoall", "cms_rhosb")
    fig, axs = plt.subplots(2, 4, figsize=(13.3, 6.4))
    # pp before and after the underlying-event subtraction
    ax = axs[0, 3]
    b = unit(cblock(t11, "pp"))
    for m in list(MODELS)[:1]:
        e, y, ey = shape(D["pp"][m], "cms_rhoall")
        ax.step(e, np.append(y, y[-1]), where="post", color="0.5", ls=":", lw=1.4, label="pp, no UE subtraction")
        e, y, ey = shape_sub(D["pp"][m], "cms3_rho", "cms3_rhobkg")
        ax.step(e, np.append(y, y[-1]), where="post", color="0.5", ls="--", lw=1.0, label="pp, eta-reflection")
        e, y, ey = rho(D["pp"][m])
        draw_model(ax, e, y, ey, m, label=False)
        add("CMS rho pp, no UE sub", "rho pp no UE sub", {mm: chi2(e, *shape(D["pp"][mm], "cms_rhoall")[1:], b) for mm in MODELS})
        add("CMS rho pp, eta-reflection", "rho pp eta-reflection",
            {mm: chi2(e, *shape_sub(D["pp"][mm], "cms3_rho", "cms3_rhobkg")[1:], b) for mm in MODELS})
    draw_data(ax, b, "CMS pp")
    ax.set_yscale("log")
    ax.legend(fontsize=7, frameon=False)
    ax.set_title("pp: effect of the UE subtraction", fontsize=9)
    axs[1, 3].axis("off")
    for col, (c, s) in enumerate([("pp", "pp"), ("0-10%", "lhc0010"), ("30-50%", "lhc3050")]):
        b = unit(cblock(t11, c))
        ax = axs[0, col]
        res = {}
        for m in MODELS:
            e, y, ey = rho(D[s][m])
            res[m] = chi2(e, y, ey, b)
            draw_model(ax, e, y, ey, m, label=(col == 0))
        add("CMS rho", f"rho {c}", res)
        draw_data(ax, b, "CMS" if col == 0 else None)
        ax.set_yscale("log")
        ax.set_title(f"$\\rho(\\Delta r)$ {c}", fontsize=9)
        if col == 0:
            continue
        b = cblock(t12, c)
        ax = axs[1, col]
        res = {}
        for m in MODELS:
            e, ya, ea = rho(D[s][m])
            _, yp, ep = rho(D["pp"][m])
            y, ey = ratio(ya, ea, yp, ep)
            res[m] = chi2(e, y, ey, b)
            draw_model(ax, e, y, ey, m, label=False)
        add("CMS rho ratio", f"rho PbPb/pp {c}", res)
        draw_data(ax, b, None)
        ax.axhline(1, color="0.6", lw=0.6)
        ax.set_ylim(0.6, 2.4)
        ax.set_xlabel(r"$\Delta r$")
        ax.set_title(f"PbPb/pp {c}", fontsize=9)
    # diagnostic: 0-10% PbPb/pp without hole subtraction (particles only)
    ax = axs[1, 0]
    for m in MODELS:
        if "cms_rhopos" not in D["lhc0010"][m]:
            break
        e, ya, ea = shape(D["lhc0010"][m], "cms_rhopos")
        _, yp, ep = shape(D["pp"][m], "cms_rhopos")
        y, ey = ratio(ya, ea, yp, ep)
        draw_model(ax, e, y, ey, m, label=False)
    draw_data(ax, cblock(t12, "0-10%"), None)
    ax.axhline(1, color="0.6", lw=0.6)
    ax.set_ylim(0.6, 2.4)
    ax.set_xlabel(r"$\Delta r$")
    ax.set_title("PbPb/pp 0-10%, no hole subtraction", fontsize=9)
    axs[0, 0].legend(fontsize=7, frameon=False)
    fig.tight_layout()
    fig.savefig(os.path.join(HERE, "fig_data_cms.png"), dpi=150)
    plt.close(fig)

    # ---------------- ALICE mass
    report.append("ALICE R = 0.2 charged-jet mass (unit area over data bins)")
    fig, axs = plt.subplots(3, 4, figsize=(13, 9))
    rows = [("pp m", "pp", lambda k: 33 + 2 * k if k < 3 else None, "ali_m"),
            ("PbPb m", "lhc0010", lambda k: 39 + 2 * k, "ali_m"),
            ("PbPb m_g", "lhc0010", lambda k: 40 + 2 * k, "ali_mg")]
    for row, (lab, s, tab, key) in enumerate(rows):
        for k, ptl in enumerate(ALIPT):
            ax = axs[row, k]
            if tab(k) is None:
                ax.axis("off")
                continue
            b = unit(hep_table("ins2845788", f"Table {tab(k)}")[0])
            res, resd = {}, {}
            for m in MODELS:
                e, y, ey = shape(D[s][m], f"{key}{k}")
                res[m] = chi2(e, y, ey, b)
                resd[m] = chi2(e, y, ey, b, diag=True)
                draw_model(ax, e, y, ey, m, label=(row == 0 and k == 0))
            add(f"ALICE {lab}", f"{lab} pT {ptl}", res, resd)
            draw_data(ax, b, "ALICE" if (row == 0 and k == 0) else None)
            ax.set_title(f"{lab}, $p_T^{{ch}}$ {ptl} GeV", fontsize=9)
            if row == 2:
                ax.set_xlabel("$m$ (GeV)")
    axs[0, 0].legend(fontsize=7, frameon=False)
    fig.tight_layout()
    fig.savefig(os.path.join(HERE, "fig_data_alice_mass.png"), dpi=150)
    plt.close(fig)

    # ---------------- ALICE angularities (Pb-Pb)
    report.append("ALICE R = 0.2 angularities, Pb-Pb (unit area over data bins)")
    fig, axs = plt.subplots(4, 8, figsize=(22, 10.5))
    for k, ptl in enumerate(ALIPT):
        for g in range(2):
            for a, (akey, alab) in enumerate(ALPHA):
                t = 1 + 8 * k + 4 * g + a
                b = unit(hep_table("ins2845788", f"Table {t}")[0])
                key = ("ali_lg" if g else "ali_l") + f"{akey}_{k}"
                ax = axs[k, 4 * g + a]
                res, resd = {}, {}
                for m in MODELS:
                    e, y, ey = shape(D["lhc0010"][m], key)
                    res[m] = chi2(e, y, ey, b)
                    resd[m] = chi2(e, y, ey, b, diag=True)
                    draw_model(ax, e, y, ey, m, label=(k == 0 and g == 0 and a == 0))
                gl = "groomed " if g else ""
                add(f"ALICE {gl}lambda", f"{gl}lambda_{alab} pT {ptl}", res, resd)
                draw_data(ax, b, "ALICE" if (k == 0 and g == 0 and a == 0) else None)
                ax.set_title(f"{gl}$\\lambda_{{{alab}}}$, {ptl} GeV", fontsize=8)
    axs[0, 0].legend(fontsize=7, frameon=False)
    fig.tight_layout()
    fig.savefig(os.path.join(HERE, "fig_data_alice_ang.png"), dpi=130)
    plt.close(fig)

    # ---------------- summary
    out = ["Step 2: model vs data, chi2 / n_bins (ALICE with correlated sources; see header)",
           "files: " + ", ".join(f"{s} {m} {D[s][m]['nfiles']}" for s in SAMPLES for m in MODELS),
           "", "Totals per observable group:"]
    out.append(f"  {'group':<28s}" + "".join(f"{MODELS[m][0][:22]:>26s}" for m in MODELS) + "   delta chi2 (EMD - default)")
    for g, v in tot.items():
        c = [v[m] for m in MODELS]
        out.append(f"  {g:<28s}" + "".join(f"{x[0]:14.1f} / {x[1]:<3d} ({x[0] / max(x[1], 1):4.2f})" for x in c)
                   + f"   {c[1][0] - c[0][0]:+8.1f}")
    out += ["", "Per table:"] + report
    open(os.path.join(HERE, "data_compare.txt"), "w").write("\n".join(out) + "\n")
    print("\n".join(out))


if __name__ == "__main__":
    main()
