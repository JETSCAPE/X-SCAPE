#!/usr/bin/env python3
"""Validate FastRootBulkWriter output against the legacy RootBulkWriter.

Run the three example configs first (from build_gpu/):
    ./runJetscape ../config/BulkFastTest/OO_one_event.xml         -> OO_test.root        (legacy)
    ./runJetscape ../config/BulkFastTest/OO_one_event_fast.xml    -> OO_test_fast.root   (native)
    ./runJetscape ../config/BulkFastTest/OO_one_event_fastgrid.xml-> OO_test_fastgrid.root (grid)

Then, from build_gpu/, with a Python that has numpy+uproot (e.g. conda `fno_env`):
    python ../config/BulkFastTest/validate.py

The system python3 on this machine has no numpy; use the conda env.
"""
import os
import sys
import numpy as np
import uproot


def load(fn):
    f = uproot.open(fn)
    t = f["t"]
    res = t["user_res"].array(library="np")
    ntau = t["ntau_freezeout"].array(library="np")
    P = lambda k: (f[k].member("fVal") if k in f else None)
    meta = {k: P(k) for k in ["nx", "ny", "neta", "nFeatures",
                              "tau_min", "dtau", "tau_stride"]}
    meta["grid_mode"] = f["grid_mode"].member("fTitle") if "grid_mode" in f else None
    meta["file_bytes"] = os.path.getsize(fn)
    return res, ntau, meta


def describe(tag, fn):
    res, ntau, m = load(fn)
    L = len(res[0])
    per = m["nx"] * m["ny"] * m["neta"] * m["nFeatures"]
    print(f"[{tag}] {fn}")
    print(f"    file_bytes = {m['file_bytes']:,}")
    print(f"    grid_mode  = {m['grid_mode']}")
    print(f"    nx={m['nx']} ny={m['ny']} neta={m['neta']} nF={m['nFeatures']} "
          f"ntau={int(ntau[0])} tau_min={m['tau_min']:.4f} dtau={m['dtau']:.5f} "
          f"tau_stride={m['tau_stride']}")
    print(f"    len(user_res[0]) = {L:,}  (= ntau*nx*ny*neta*nF ? {L == int(ntau[0])*per})")
    return res, ntau, m


def main():
    here = os.getcwd()
    need = ["OO_test.root", "OO_test_fastgrid.root", "OO_test_fast.root"]
    missing = [f for f in need if not os.path.exists(f)]
    if missing:
        print("Missing ROOT files (run the configs first, from build_gpu/):", missing)
        sys.exit(1)

    rL, nL, mL = describe("legacy", "OO_test.root")
    rG, nG, mG = describe("grid",   "OO_test_fastgrid.root")
    rN, nN, mN = describe("native", "OO_test_fast.root")

    print("\n--- grid mode vs legacy (should match closely) ---")
    a = np.asarray(rL[0], np.float64)
    b = np.asarray(rG[0], np.float64)
    if a.size != b.size:
        print(f"    SHAPE MISMATCH: legacy {a.size} vs grid {b.size}")
    else:
        d = np.abs(a - b)
        # relative diff on the energy channel (feature 0), where values are large
        eL = a[0::4]
        eG = b[0::4]
        rel = np.abs(eL - eG) / np.maximum(np.abs(eL), 1e-6)
        print(f"    elements          = {a.size:,}")
        print(f"    max |abs diff|    = {d.max():.3e}")
        print(f"    mean |abs diff|   = {d.mean():.3e}")
        print(f"    energy max value  = {eL.max():.5g}  (max |abs diff| / energy max = {d.max()/eL.max():.2e})")
        print(f"    energy max |rel|  = {rel.max():.3e}")
        print(f"    bit-for-bit equal = {np.array_equal(a, b)}")

    print("\n--- native mode sanity (full MUSIC grid) ---")
    E = np.asarray(rN[0], np.float64).reshape(
        int(nN[0]), mN["nx"], mN["ny"], mN["neta"], mN["nFeatures"])[..., 0]
    print(f"    energy min/median/max = {E.min():.4g} / {np.median(E):.4g} / {E.max():.5g}")
    print(f"    legacy energy max     = {np.asarray(rL[0], np.float64)[0::4].max():.5g}")
    print(f"    native nan={bool(np.isnan(E).any())} inf={bool(np.isinf(E).any())}")


if __name__ == "__main__":
    main()
