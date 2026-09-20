#!/usr/bin/env python
"""Phase 2/3 of FNO4d/PLAN_fv_vs_music.md: run ONE MUSIC leg on the shared fast_data IC.

Builds the task list in Python so the initial condition can be injected from an HDF5 file
without any new C++:

    FvMusicInitialState (this file)  ->  NullPreDynamics  ->  MUSIC

`NullPreDynamics` copies `ini->GetEntropyDensityDistribution()` straight into `e_` as ENERGY
density with u = (1,0,0,0), pi = 0, Pi = 0 (NullPreDynamics.cc:36-56), which is exactly the
state the FV solver starts from.  MUSIC's Initial_profile 42 reader then indexes that vector
as `idx = (ny*neta)*ix + neta*iy + ieta` (init.cpp:1239) -- plain C order on (nx, ny, neta),
so `e0[event].ravel()` is the right thing to hand over.

FastRootBulkWriter is NOT in the task list by default: it calls
`clear_hydro_info_from_memory()` at the end of its Exec (FastRootBulkWriter.cc:271), which
would destroy the native store before we can read it.  We read the store directly with
`get_native_evolution_numpy()` inside `per_event_loop`, which yields after the event has run
but before its memory is released.  Pass --root to add the writer as well (then the numpy
export is taken first).

Usage
-----
    cd $XSCAPE/build_gpu
    MUSIC_FORCE_CPU=1 python ../config/FVvsMUSIC/run_music_leg.py \
        --leg ideal_conformal \
        --ic  ~/FNO4d/workflow_fastdata/out/ic_fvmusic_AuAu_b0.h5 \
        --out out/music_ideal_conformal.npz

    # the Phase 2 gate -- asymmetric probe, two tau steps, exact IC round-trip check
    MUSIC_FORCE_CPU=1 python ../config/FVvsMUSIC/run_music_leg.py \
        --leg ideal_conformal --probe \
        --ic  ~/FNO4d/workflow_fastdata/out/ic_fvmusic_probe.h5 \
        --out out/music_probe.npz
"""

from __future__ import annotations

import argparse
import os
import shutil
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
XSCAPE = os.path.abspath(os.path.join(HERE, "..", ".."))
LEGS = ("ideal_conformal", "ideal_eos91", "is_eos91")


def _import_jetscape():
    """Import the PyJetscape package, adding its in-tree location to sys.path."""
    for p in (os.path.join(XSCAPE, "external_packages", "js-contrib", "contribs",
                           "PyJetscape", "python"),
              os.path.join(XSCAPE, "build_gpu", "external_packages", "js-contrib",
                           "contribs", "PyJetscape", "python")):
        if os.path.isdir(p) and p not in sys.path:
            sys.path.insert(0, p)
    import jetscape                                   # noqa: E402
    from jetscape.run_jetscape import per_event_loop  # noqa: E402
    return jetscape, per_event_loop


def make_ic_class(jetscape):
    """Build the Python InitialState subclass against the imported bindings.

    Only `Exec` is overridden.  `Init` is deliberately left alone: the trampoline binds
    `Init` to `&JetScapeModuleBase::Init` (bind_framework.cc:162), which dispatches
    virtually, so a Python `Init` that called `super().Init()` would re-enter itself.
    Leaving it to C++ also means `InitialState::Init()` reads <IS><grid_*> from the XML
    (InitialState.cc:24-40) as usual, which we then cross-check against the IC file.
    `InitialState::ExecuteTask()` is an empty stub, so overriding `Exec` without chaining
    to the base loses nothing.
    """

    class FvMusicInitialState(jetscape.InitialState):
        def __init__(self, ic_path, event=0, verbose=True):
            super().__init__()
            self.SetId("FvMusicInitialState")
            self.ic_path = ic_path
            self.event = int(event)
            self.verbose = verbose
            self.e0 = None          # (nx, ny, neta) float64, GeV/fm^3
            self.meta = None

        def _load(self):
            import h5py
            with h5py.File(self.ic_path, "r") as f:
                fmt = f.attrs.get("format", b"")
                fmt = fmt.decode() if isinstance(fmt, bytes) else fmt
                if fmt != "mc_glauber_tilted/initial_state":
                    raise ValueError(f"{self.ic_path}: unexpected format {fmt!r}")
                self.meta = {k: (v.decode() if isinstance(v, bytes) else v)
                             for k, v in f.attrs.items()}
                e = f["e0"][self.event]
            # float64 and C-contiguous: set_entropy_density_from_numpy forcecasts, but
            # doing it here makes the array we assert against identical to the one sent.
            self.e0 = np.ascontiguousarray(e, dtype=np.float64)

        def Exec(self):
            if self.e0 is None:
                self._load()
            m = self.meta
            nx, ny, neta = (int(m["nx"]), int(m["ny"]), int(m["neta"]))
            dx, dy, deta = (float(m["dx"]), float(m["dy"]), float(m["deta"]))
            if self.e0.shape != (nx, ny, neta):
                raise ValueError(f"IC array {self.e0.shape} != attrs ({nx},{ny},{neta})")
            if nx != ny:
                raise ValueError(
                    f"Initial_profile 42 recovers nx = sqrt(size/neta) and sets ny = nx "
                    f"(init.cpp:135-137); the grid must be square, got {nx}x{ny}")

            # MUSIC's axis is x = i*dx - nx*dx/2, and GetXSize() = ceil(2*grid_max/step),
            # so grid_max = n*d/2 makes the framework, MUSIC and this array agree.
            self.SetRanges(nx * dx / 2.0, ny * dy / 2.0, neta * deta / 2.0)
            self.SetSteps(dx, dy, deta)

            got = (self.GetXSize(), self.GetYSize(), self.GetZSize())
            if got != (nx, ny, neta):
                raise ValueError(f"InitialState grid {got} != IC grid {(nx, ny, neta)}")

            self.set_entropy_density_from_numpy(self.e0)   # energy density, GeV/fm^3

            if self.verbose:
                ijk = np.unravel_index(int(np.argmax(self.e0)), self.e0.shape)
                flat = (ny * neta) * ijk[0] + neta * ijk[1] + ijk[2]
                print(f"[FvMusicInitialState] injected {self.e0.shape} from "
                      f"{os.path.basename(self.ic_path)} event {self.event}", flush=True)
                print(f"[FvMusicInitialState]   tau0 = {float(m['tau0'])} fm/c, "
                      f"d = ({dx}, {dy}, {deta}) fm, sizes = {got}", flush=True)
                print(f"[FvMusicInitialState]   e_max = {self.e0.max():.6g} GeV/fm^3 "
                      f"at index {ijk}, MUSIC flat idx {flat}", flush=True)

        def Clear(self):
            pass

    return FvMusicInitialState


def xml_int(user_xml, path):
    import xml.etree.ElementTree as ET
    node = ET.parse(user_xml).getroot().find(path)
    if node is None:
        raise ValueError(f"{user_xml}: no <{path}>")
    return int(node.text)


def check_taus_matches_ic(user_xml, tau0_ic):
    """<Preequilibrium><taus>, not <Initial_time_tau_0>, is MUSIC's tau0 on this path."""
    import xml.etree.ElementTree as ET
    node = ET.parse(user_xml).getroot().find("./Preequilibrium/taus")
    if node is None:
        raise ValueError(f"{user_xml}: no <Preequilibrium><taus>")
    taus = float(node.text)
    if abs(taus - tau0_ic) > 1e-12:
        raise ValueError(
            f"<Preequilibrium><taus> = {taus} but the IC file says tau0 = {tau0_ic}.\n"
            f"On the Initial_profile 42 path MUSIC takes tau0 from <taus> "
            f"(MusicWrapper.cc:391 -> PreequilibriumDynamics.cc:58), so the two codes "
            f"would start at different times with no warning.")
    return taus


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--leg", choices=LEGS, required=True)
    ap.add_argument("--ic", required=True, help="the shared IC HDF5 from fvmusic_make_ic.py")
    ap.add_argument("--out", required=True, help="output .npz")
    ap.add_argument("--event", type=int, default=0, help="event index inside the IC file")
    ap.add_argument("--workdir", default=os.path.join(XSCAPE, "build_gpu"),
                    help="run here; MUSIC resolves EOS/ and music_input relative to it")
    ap.add_argument("--main-xml", default=os.path.join(XSCAPE, "config", "jetscape_main.xml"))
    ap.add_argument("--user-xml", default=None, help="override the per-leg XML")
    ap.add_argument("--tau-stride", type=int, default=1)
    ap.add_argument("--probe", action="store_true",
                    help="Phase 2 gate: cap the run at a couple of tau steps and assert "
                         "that MUSIC's frame 0 reproduces the injected IC exactly")
    ap.add_argument("--tau-max", type=float, default=None,
                    help="Total_evolution_time_tau, i.e. the DURATION of the evolution from "
                         "tau0 [fm/c].  Always set this: MUSIC only self-terminates when "
                         "`frozen == 1 && tau > source_tau_max` (evolve.cpp:396), which does "
                         "not trigger in a source-free run, so the default 30 fm would be "
                         "evolved in full -- hundreds of vacuum frames that the FV solver "
                         "cannot follow stably anyway.  Default: 0.7 with --probe, else 6.0 "
                         "for the conformal leg and 12.0 for the EOS-91 legs.")
    ap.add_argument("--probe-rtol", type=float, default=2e-6,
                    help="tolerance on the frame-0 round trip (float32 store)")
    ap.add_argument("--root", action="store_true",
                    help="also add FastRootBulkWriter to the task list")
    a = ap.parse_args(argv)

    user_xml = a.user_xml or os.path.join(HERE, f"fv_vs_music_{a.leg}.xml")
    ic_path = os.path.abspath(os.path.expanduser(a.ic))
    out_path = os.path.abspath(os.path.expanduser(a.out))
    workdir = os.path.abspath(os.path.expanduser(a.workdir))
    for p in (user_xml, a.main_xml, ic_path):
        if not os.path.exists(p):
            raise SystemExit(f"not found: {p}")

    import h5py
    with h5py.File(ic_path, "r") as f:
        tau0_ic = float(f.attrs["tau0"])
    check_taus_matches_ic(user_xml, tau0_ic)

    # MpiMusic rewrites EOS_to_use / Include_Bulk_Visc in the MUSIC input file in place
    # (MusicWrapper.cc:124,243), so run off a copy and leave the template pristine.
    # One copy per leg: MpiMusic rewrites the file in place, so two legs sharing a
    # filename would clobber each other's EOS_to_use if run concurrently.  The <MUSIC_input_file>
    # tag in each leg XML names the matching per-leg copy.
    src = os.path.join(HERE, "music_input_fv")
    dst = os.path.join(workdir, f"music_input_fv_{a.leg}")
    shutil.copyfile(src, dst)
    tau_max = a.tau_max
    if tau_max is None:
        tau_max = 0.7 if a.probe else (6.0 if a.leg == "ideal_conformal" else 12.0)
    import re as _re
    with open(dst) as fh:
        txt = fh.read()
    txt = _re.sub(r"^Total_evolution_time_tau .*$",
                  f"Total_evolution_time_tau {tau_max}   # set by run_music_leg.py",
                  txt, flags=_re.M)
    # EOS_to_use MUST be written into the file, not left to <MUSIC><EOS>.
    # MUSIC builds its EoS in the constructor initialiser list, `eos(DATA.whichEOS)`
    # (music.cpp:24-26), from the value read out of THIS file.  The wrapper's later
    # set_parameter("EOS", v) only assigns DATA.whichEOS (read_in_parameters.cpp:1012) and
    # never rebuilds that object, and update_music_input_parameter() then rewrites the file
    # for iSS's benefit -- so the file ends up LOOKING right while the run used the old EoS.
    eos_id = xml_int(user_xml, "./Hydro/MUSIC/EOS")
    txt = _re.sub(r"^EOS_to_use .*$",
                  f"EOS_to_use {eos_id}   # set by run_music_leg.py from <MUSIC><EOS>",
                  txt, flags=_re.M)
    with open(dst, "w") as fh:
        fh.write(txt)
    print(f"MUSIC input : {dst}  (copy of {src})")
    print(f"tau window  : {tau0_ic} .. {tau0_ic + tau_max} fm/c "
          f"(Total_evolution_time_tau = {tau_max})")
    print(f"EOS_to_use  : {eos_id}  (written into the input file, see the comment above)")
    print(f"user XML    : {user_xml}")
    print(f"IC          : {ic_path}  (tau0 = {tau0_ic})")
    print(f"workdir     : {workdir}")
    print(f"MUSIC_FORCE_CPU = {os.environ.get('MUSIC_FORCE_CPU', '<unset>')}")

    jetscape, per_event_loop = _import_jetscape()
    ICClass = make_ic_class(jetscape)

    os.chdir(workdir)
    ini = ICClass(ic_path, event=a.event)
    modules = [ini, jetscape.create_module("NullPreDynamics"),
               jetscape.create_module("MUSIC")]
    if a.root:
        modules.append(jetscape.create_module("FastRootBulkWriter"))

    saved = None
    for js in per_event_loop(a.main_xml, user_xml, modules=modules, n_events=1):
        hydro = jetscape.JetScapeSignalManager.Instance().GetHydroPointer()
        if hydro is None:
            raise SystemExit("no hydro module registered")
        g = hydro.get_bulk_info()
        arr = hydro.get_native_evolution_numpy(a.tau_stride)   # (ntau,nx,ny,neta,4)
        saved = dict(
            arr=arr,
            tau_min=np.float64(g.tau_min), dtau=np.float64(g.dtau * a.tau_stride),
            ntau=np.int64(arr.shape[0]),
            x_min=np.float64(g.x_min), dx=np.float64(g.dx),
            y_min=np.float64(g.y_min), dy=np.float64(g.dy),
            eta_min=np.float64(g.eta_min), deta=np.float64(g.deta),
            nx=np.int64(g.nx), ny=np.int64(g.ny), neta=np.int64(g.neta),
            tau_stride=np.int64(a.tau_stride),
            leg=np.str_(a.leg), ic=np.str_(ic_path), tau0_ic=np.float64(tau0_ic),
        )
        break

    if saved is None:
        raise SystemExit("the event loop yielded nothing")

    arr = saved["arr"]
    print(f"\nnative store: {arr.shape} {arr.dtype}  "
          f"tau = {saved['tau_min']} + k*{saved['dtau']}, ntau = {saved['ntau']}")
    print(f"  grid: x_min={saved['x_min']} dx={saved['dx']}  "
          f"eta_min={saved['eta_min']} deta={saved['deta']}")
    print(f"  e: min={arr[..., 0].min():.4g} max={arr[..., 0].max():.4g} GeV/fm^3; "
          f"finite = {bool(np.isfinite(arr).all())}")

    rc = 0
    if a.probe:
        rc = _probe_gate(arr, saved, ini, a.probe_rtol)

    os.makedirs(os.path.dirname(out_path) or ".", exist_ok=True)
    np.savez_compressed(out_path, **saved)
    print(f"\nwrote {out_path}")
    return rc


def _probe_gate(arr, saved, ini, rtol):
    """Assert MUSIC's frame 0 reproduces the injected IC, cell for cell."""
    print("\n--- Phase 2 gate: IC round trip ---")
    e_in = ini.e0
    e_out = arr[0, ..., 0].astype(np.float64)
    ok = True

    if e_out.shape != e_in.shape:
        print(f"  FAIL shape: MUSIC {e_out.shape} vs IC {e_in.shape}")
        return 1

    i_in = np.unravel_index(int(np.argmax(e_in)), e_in.shape)
    i_out = np.unravel_index(int(np.argmax(e_out)), e_out.shape)
    print(f"  arg-max cell : IC {i_in}  MUSIC {i_out}   "
          f"{'OK' if i_in == i_out else 'FAIL -- transposed or shifted axis'}")
    ok &= (i_in == i_out)

    # The native store is float32, so compare at float32 resolution.
    denom = np.maximum(np.abs(e_in), 1e-30)
    rel = np.abs(e_out - e_in) / denom
    # The vacuum floor is clamped by MUSIC to eps >= 1e-16 fm^-4, so judge on the fluid.
    fluid = e_in > 1e-3 * e_in.max()
    print(f"  peak         : IC {e_in.max():.8g}  MUSIC {e_out.max():.8g}")
    print(f"  max rel diff : all cells {rel.max():.3e}, "
          f"fluid cells (e > 1e-3 e_max, n = {fluid.sum()}) {rel[fluid].max():.3e}")
    ok &= bool(rel[fluid].max() <= rtol)

    # An eta offset would show as an asymmetry between the eta axis of the two arrays.
    eta_in = (np.arange(e_in.shape[2]) - 0.5 * (e_in.shape[2] - 1)) * float(saved["deta"])
    eta_mus = float(saved["eta_min"]) + np.arange(e_in.shape[2]) * float(saved["deta"])
    print(f"  eta axis     : fast_data [{eta_in[0]:+.5f}, {eta_in[-1]:+.5f}]  "
          f"MUSIC [{eta_mus[0]:+.5f}, {eta_mus[-1]:+.5f}]  "
          f"(offset {eta_mus[0] - eta_in[0]:+.5f} fm -- expected -deta/2 = "
          f"{-float(saved['deta']) / 2:+.5f}; see mismatch 2 in the plan)")

    print(f"  => {'PASS' if ok else 'FAIL'} (tolerance {rtol:.1e})")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
