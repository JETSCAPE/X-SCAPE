// JEWEL 2.6.0 C++ port -- interface for use as a final-state shower in JETSCAPE.
//
// This header is plain C++17 (JETSCAPE's standard); the JEWEL library itself is
// built with a newer standard. JEWEL keeps its global state (PYTHIA 6 common
// blocks, lookup tables), so there is one JEWEL instance per process and the
// functions below must not be called from several threads at once.
//
// Use:
//   1. jewel::api::initialise(settings, cell_query) once. This reads the JEWEL
//      parameter file (same labels as JEWEL 2.6.0) and computes or reads the
//      lookup tables (splitint.dat, pdfs.dat, xsecs.dat).
//   2. jewel::api::next_event(tau0) at the start of every JETSCAPE event.
//   3. jewel::api::shower(...) for every hard parton. It runs JEWEL's complete
//      medium-modified shower (vacuum splittings, elastic and inelastic
//      scatterings, LPM interference) for that parton in the medium described by
//      cell_query, and returns the final partons, the medium recoils and the
//      thermal momenta of the scattering centres before the scattering ("holes").
//
// The medium replaces medium-simple.f: temperature and flow velocity come from
// cell_query (e.g. a JETSCAPE hydro), scattering centres are thermal quarks and
// gluons sampled in the fluid rest frame and boosted with the cell velocity, as
// in LBT. Debye mass, density and the T_c cut follow medium-simple.f.
#pragma once

#include <functional>
#include <string>
#include <vector>

namespace jewel::api {

// Local medium at a lab-frame space-time point (t, x, y, z) [fm].
struct MediumCell {
  double temperature = 0.0;  // GeV; <= 0 means no medium
  double vx = 0.0, vy = 0.0, vz = 0.0;  // fluid velocity in the lab frame
};
using CellQuery = std::function<MediumCell(double t, double x, double y, double z)>;

struct MediumSettings {
  // Upper bound of the temperature anywhere in the medium [GeV]. The tables and
  // the veto algorithm for the scattering rate need it; cells above it are
  // counted (n_above_tmax) and JEWEL logs "weight > 1".
  double tmax = 0.6;
  double tc = 0.17;          // no scattering centres below T_c [GeV] (TC)
  int nf = 3;                // thermal light flavours (NF)
  double mdfactor = 0.45;    // minimum Debye mass [GeV] (MDFACTOR)
  double mdscalefac = 1.0;   // m_D = 3 MDSCALEFAC T (MDSCALEFAC)
  double lifetime = 20.0;    // no medium interactions after this lab time [fm]
  double tau0 = 0.6;         // medium starts at proper time tau0 [fm] (overridden by next_event)
  // before tau0: 0 = no medium (as in LBT), 1 = T(tau0) * tau / tau0 at the same
  // (x, y, eta_s) (the linear rise of medium-simple.f)
  int pre_tau0 = 0;
  // bound on the flux factor gamma (1 - v.beta) of the scattering rate (enters
  // the veto algorithm only)
  double fluxmax = 6.0;
  double etamax = 3.1;       // |eta_s| beyond which there is no medium
  bool vacuum = false;       // no medium at all (cell_query is not called)
};

struct Settings {
  std::string param_file;            // JEWEL parameter file ("" = JEWEL defaults)
  std::string run_param_file = "jewel_params.dat";  // effective parameter file written here
  std::vector<std::string> extra_params;  // "LABEL value" lines appended (override the file)
  long seed = 0;                     // PYR seed, 0 <= seed < 900000000
  MediumSettings medium;
};

// One parton out of JEWEL. status follows LBT: 0 shower parton, 1 medium recoil,
// -1 thermal momentum of a scattering centre before the scattering (to be
// subtracted).
struct OutParton {
  int id = 0;
  int status = 0;
  double p[4] = {0, 0, 0, 0};  // E, px, py, pz [GeV]
  double x[4] = {0, 0, 0, 0};  // t, x, y, z of production [fm]
  double mass = 0.0;           // invariant mass [GeV]
  int col = 0, acol = 0;       // JEWEL colour indices (0 = none)
};

struct Stats {
  long n_showers = 0, n_discarded = 0, n_above_tmax = 0, n_scatterings = 0;
  long n_above_fluxmax = 0;   // rate queries with gamma (1 - v.beta) > fluxmax
  double max_flux = 0.0;      // largest gamma (1 - v.beta) seen
  long n_holes_dropped = 0;   // stored scattering centres whose scattering was undone
};

void initialise(const Settings& settings, CellQuery cell_query);
void next_event(double tau0);
// Shower one parton with flavour id (PDG), four-momentum p = (E, px, py, pz) and
// production point x = (t, x, y, z). qmax is the maximum virtuality [GeV] (JEWEL
// uses p_T of the hard scattering). Returns false if JEWEL discarded the shower
// (out is then empty).
bool shower(int id, const double p[4], const double x[4], double qmax, std::vector<OutParton>& out);
const Stats& stats();

}  // namespace jewel::api
