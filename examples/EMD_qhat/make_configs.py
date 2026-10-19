#!/usr/bin/env python3
"""Write every JETSCAPE configuration of the EMD qhat study, by name.

    python3 make_configs.py [group ...]     (default: all groups)

Groups and what they reproduce (task lists go to tasks/<group>.txt, XMLs to
configs/<name>.xml):

  bricks   parton-level bricks of brick_figures.pdf (output in bricks/<set>/)
  pp       pp baselines: LHC 5.02 TeV (16 x 25k events), RHIC 200 GeV (8 x 25k)
  routeA   MATTER+MyEloss with EMD (lambda_t=256) and HTL (K=16.4) qhat, and the
           JETSCAPE tune; RHIC and LHC, 0-10% and 30-50%
  routeB   MATTER Type 9 + LBT with the EMD radiative qhat (lambda_t = 2, 2.6, 8)
           and the extra tune runs it is paired with
  lamcal   lambda_t = 1, 1.5, 3 at the LHC, 0-10% and 30-50% (chi^2 calibration)
  bestfit  the best fit lambda_t = 2.4 against the default MATTER+LBT, RHIC and LHC,
           0-10% and 30-50%, plus pp; run with KEEP_EVENTS=1 to keep the events
           (these are the inputs of x_compare.py: x_J, D(z), rho(r), jet mass, ...,
           and, through jetana_data, of data_compare.py)
  ppvar    pp hadronization study (pp_variants.py): the LHC pp of bestfit, seeds 1-8,
           with PYTHIA MPI off, colored or hybrid hadronization, plus 2 x 64 short
           colored runs (see group_ppvar); run with KEEP_EVENTS=1

Hydro tasks write their final-state hadrons to the named pipe
/tmp/fifo/<name>, which run_task.sh reads with jetana; bricks write
<work>/bricks/<set>/<name>_final_state_partons.dat. Each hydro event is
reused for 2500 hard events; tasks of different models with the same seed use
the same hydro event and the same Pythia hard scatterings ("paired").
"""
import os
import sys

import mkcfg

HERE = os.path.dirname(os.path.abspath(__file__))

LHC010 = list(range(16, 46))          # route A and the first route-B runs
LHC010_EXTRA = list(range(331, 361))  # added for route B / lambda_t calibration
LHC3050 = list(range(116, 151))
RHIC010 = list(range(56, 96))
RHIC010_EXTRA = list(range(171, 331))
RHIC3050 = list(range(131, 171))


def hydro(system, prefix, model, tag, seeds, cent, **kw):
    """One task per hydro seed: 2500 jet events on one reused hydro event."""
    clow, chigh = cent
    out = {}
    for s in seeds:
        name = f"{prefix}_{tag}_h{s}"
        out[name] = mkcfg.collision2(system, model, f"/tmp/fifo/{name}", 2500, s, nreuse=2500,
                                     clow=clow, chigh=chigh, **kw)
    return out


def group_pp():
    out = {}
    # LHC pp was made with the first-generation template (mkcfg.collision: no
    # explicit maxT or eCMforHadronization), kept here so it reproduces exactly
    for s in range(1, 17):
        out[f"pp_s{s}"] = mkcfg.collision("vac", f"/tmp/fifo/pp_s{s}", (10, 500), 25000, 100 + s)
    for s in range(1, 9):
        out[f"rhic_pp_s{s}"] = mkcfg.collision2("RHIC", "vac", f"/tmp/fifo/rhic_pp_s{s}", 25000, 300 + s)
    return out


SAMPLES_A = (("LHC", "lhc", LHC010, (0, 10)), ("LHC", "lhc_c3050", LHC3050, (30, 50)),
             ("RHIC", "rhic", RHIC010, (0, 10)), ("RHIC", "rhic_c3050", RHIC3050, (30, 50)))


def group_routeA():
    out = {}
    for system, prefix, seeds, cent in SAMPLES_A:
        out.update(hydro(system, prefix, "tune", "tune", seeds, cent))
        out.update(hydro(system, prefix, "emdv", "emdv_l256", seeds, cent, lam=256.0))
        out.update(hydro(system, prefix, "htlv", "htlv_k16", seeds, cent, K=16.4))
    return out


def group_routeB():
    out = {}
    # lambda_t = 2: all four samples, with the extra statistics
    for system, prefix, seeds, cent in (("LHC", "lhc", LHC010 + LHC010_EXTRA, (0, 10)),
                                        ("LHC", "lhc_c3050", LHC3050, (30, 50)),
                                        ("RHIC", "rhic", RHIC010 + RHIC010_EXTRA, (0, 10)),
                                        ("RHIC", "rhic_c3050", RHIC3050, (30, 50))):
        out.update(hydro(system, prefix, "lbtemd", "lbtemd_l2", seeds, cent, lam=2.0))
    # lambda_t = 2.6 and 8 (0-10% only)
    for lam, tag in ((2.6, "lbtemd_l2p6"), (8.0, "lbtemd_l8")):
        out.update(hydro("LHC", "lhc", "lbtemd", tag, LHC010, (0, 10), lam=lam))
        out.update(hydro("RHIC", "rhic", "lbtemd", tag, RHIC010, (0, 10), lam=lam))
    # tune runs on the extra seeds (route A has the others)
    out.update(hydro("LHC", "lhc", "tune", "tune", LHC010_EXTRA, (0, 10)))
    out.update(hydro("RHIC", "rhic", "tune", "tune", RHIC010_EXTRA, (0, 10)))
    return out


def group_lamcal():
    out = {}
    for lam, tag in ((1.0, "lbtemd_l1"), (1.5, "lbtemd_l1p5"), (3.0, "lbtemd_l3")):
        out.update(hydro("LHC", "lhc", "lbtemd", tag, LHC010 + LHC010_EXTRA, (0, 10), lam=lam))
        out.update(hydro("LHC", "lhc_c3050", "lbtemd", tag, LHC3050, (30, 50), lam=lam))
    return out


def group_bricks(work):
    """Name -> XML; names are '<set>/<run>'. 100 GeV parton at eta = 0; seed 7
    (seed 5 for the path-length scan)."""
    out = {}

    def b(sset, name, model, nev=1000, seed=7, **kw):
        path = os.path.join(work, "bricks", sset, name)
        out[f"{sset}/{name}"] = mkcfg.brick(model, path, nev=nev, seed=seed, **kw)

    Ts = (0.18, 0.2, 0.25, 0.3, 0.4, 0.5)
    b("tscan", "vac", "vac", T=0.3)
    for T in Ts:                                    # route A at lambda_t = 1 (Figs 2-4)
        b("tscan", f"tune_T{T:g}", "tune", T=T, L=4.0)
        b("tscan", f"emd_T{T:g}", "emd", T=T, L=4.0)
        b("tscan", f"emdv_T{T:g}", "emdv", T=T, L=4.0)
    for T in (0.2, 0.3, 0.4, 0.5):                  # lambda_t scan, conformal control, HTL
        for lam in (4, 9, 16):
            b("lscan", f"emdv_l{lam}_T{T:g}", "emdv", T=T, L=4.0, lam=float(lam))
        b("lscan", f"emdc_l9_T{T:g}", "emdv", T=T, L=4.0, lam=9.0, ctref=0.3)
        for K in (1, 3):
            b("htlscan", f"htl_k{K}_T{T:g}", "htlv", T=T, L=4.0, K=K)
    for L in (1, 2, 3, 6, 8):                       # path length at T = 0.3 GeV (Fig 5)
        b("diag", f"L{L}_tune", "tune", nev=2000, seed=5, T=0.3, L=float(L))
        b("diag", f"L{L}_emdv", "emdv", nev=2000, seed=5, T=0.3, L=float(L), lam=9.0)
        b("htlscan", f"htl_k3_L{L}", "htlv", nev=2000, seed=5, T=0.3, L=float(L), K=3)
    for pid in (1, 21):                             # route B (Fig 7)
        b("opt1", f"vac_pid{pid}", "vac", T=0.3, pid=pid)
        b("opt1", f"tune_pid{pid}", "tune", T=0.3, L=4.0, pid=pid)
        for lam in (1, 2, 4, 8):
            b("opt1", f"lbtemd_l{lam}_pid{pid}", "lbtemd", T=0.3, L=4.0, pid=pid, lam=float(lam))
    return out


def group_bestfit():
    """EMD at the best-fit lambda_t = 2.4 vs the default (tune), new task names
    (prefix x) so they sit next to the earlier runs; default and pp configs are
    identical to the production ones of the same seed."""
    out = {}
    for s in range(1, 17):
        n = f"ppx_s{s}"
        out[n] = mkcfg.collision("vac", f"/tmp/fifo/{n}", (10, 500), 25000, 100 + s)
    models = (("tune", "tune", {}), ("lbtemd", "lbtemd_l2p4", {"lam": 2.4}))
    for system, prefix, seeds, cent in (("LHC", "lhcx", LHC010 + LHC010_EXTRA, (0, 10)),
                                        ("LHC", "lhcx_c3050", LHC3050, (30, 50)),
                                        ("RHIC", "rhicx", RHIC010 + RHIC010_EXTRA[:40], (0, 10)),
                                        ("RHIC", "rhicx_c3050", RHIC3050, (30, 50))):
        for s in seeds:
            for model, tag, kw in models:
                n = f"{prefix}_{tag}_h{s}"
                out[n] = mkcfg.collision2(system, model, f"/tmp/fifo/{n}", 2500, s, nreuse=2500,
                                          clow=cent[0], chigh=cent[1], **kw)
    for s in range(1, 9):
        n = f"rhic_ppx_s{s}"
        out[n] = mkcfg.collision2("RHIC", "vac", f"/tmp/fifo/{n}", 25000, 300 + s)
    return out


def group_ppvar():
    """LHC pp with the hadronization and MPI varied (pp_variants.py). PythiaGun
    switches MPI on in the code; the LinesToRead line overrides it."""
    out = {}
    for s in range(1, 9):
        base = mkcfg.collision("vac", "/tmp/fifo/NAME", (10, 500), 25000, 100 + s)
        for v in ("nompi", "colored", "colored_nompi", "hybrid"):
            n = f"ppv_{v}_s{s}"
            x = base.replace("/tmp/fifo/NAME", f"/tmp/fifo/{n}")
            had = "colored" if v.startswith("colored") else "hybrid" if v == "hybrid" else "colorless"
            x = x.replace("<name>colorless</name>", f"<name>{had}</name>")
            if v.endswith("nompi"):
                x = x.replace("PhaseSpace:bias2SelectionRef = 10.0\n",
                              "PhaseSpace:bias2SelectionRef = 10.0\n      PartonLevel:MPI = off\n")
            out[n] = x
    # ColoredHadronization never recovers after its first PYTHIA junction failure
    # in a run (it then writes unhadronized partons, which jetana_data skips),
    # so the colored samples are mostly made of many short runs
    for s in range(1, 65):
        for v, nev, seed0 in (("coloredr", 600, 2000), ("colorednompir", 300, 3000)):
            n = f"ppv_{v}_s{s}"
            x = mkcfg.collision("vac", f"/tmp/fifo/{n}", (10, 500), nev, seed0 + s)
            x = x.replace("<name>colorless</name>", "<name>colored</name>")
            if v == "colorednompir":
                x = x.replace("PhaseSpace:bias2SelectionRef = 10.0\n",
                              "PhaseSpace:bias2SelectionRef = 10.0\n      PartonLevel:MPI = off\n")
            out[n] = x
    return out


def main(groups, work=HERE):
    gens = {"bricks": lambda: group_bricks(work), "pp": group_pp, "routeA": group_routeA,
            "routeB": group_routeB, "lamcal": group_lamcal, "bestfit": group_bestfit,
            "ppvar": group_ppvar}
    os.makedirs(os.path.join(work, "configs"), exist_ok=True)
    os.makedirs(os.path.join(work, "tasks"), exist_ok=True)
    os.makedirs(os.path.join(work, "logs"), exist_ok=True)   # SLURM output of slurm/array.sbatch
    for g in groups or gens:
        cfgs = gens[g]()
        for name, xml in cfgs.items():
            path = os.path.join(work, "configs", name + ".xml")
            os.makedirs(os.path.dirname(path), exist_ok=True)
            os.makedirs(os.path.join(work, "bricks", os.path.dirname(name)), exist_ok=True) if g == "bricks" else None
            open(path, "w").write(xml)
        open(os.path.join(work, "tasks", f"{g}.txt"), "w").write("\n".join(cfgs) + "\n")
        print(f"{g:7s} {len(cfgs):5d} configurations -> tasks/{g}.txt")


if __name__ == "__main__":
    main(sys.argv[1:])
