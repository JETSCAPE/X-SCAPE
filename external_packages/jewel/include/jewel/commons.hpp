// JEWEL 2.6.0 C++ port -- global state of the JEWEL core.
//
// Every Fortran COMMON block of jewel-2.6.0.f is one namespace
// jewel::cb::<blockname> holding its members as inline variables (static
// storage, zero-initialised like Fortran COMMON storage in .bss). A ported
// routine pulls in exactly the blocks the Fortran routine declared:
//     using namespace jewel::cb::param;
//     using namespace jewel::py::pyjets;
// Member names are the Fortran names in lower case (PYTHIA blocks keep the
// conventional upper-case names, see pythia6.hpp). Arrays keep their Fortran
// shape and 1-based (or declared) index bounds.
//
// None of these blocks is shared with the Fortran PYTHIA library except
// those declared in pythia6.hpp (PYJETS, PYDAT1, PYDAT3, PYSUBS, PYPARS,
// PYDATR, PDFAS, WEXPO), so plain C++ storage is used here.
#pragma once

#include <string>

#include "jewel/farray.hpp"
#include "jewel/pythia6.hpp"

namespace jewel::cb {

constexpr int kLines = f77::kMaxLines;  // 23000

// C--pdfset
namespace pdf { inline int pdfset = 0; }
// C--number of protons
namespace np { inline int nproton = 0, mass = 0; }
// C--Parameter common block
namespace param {
inline double q0 = 0, lps = 0, lqcd = 0, ltime = 0, scalefacm = 0, recsoftcut = 0, rechardcut = 0;
inline int angord = 0, nf = 0, kinmode = 0, recmode = 0;
inline bool aofirstsplit = false, scatrecoil = false, allhad = false, compress = false, mpifsr = false;
}  // namespace param
// C--coherence parameters
namespace cohparams {
inline double cohlengthfac = 0, pcohrej = 0;
inline bool cohscat = false;
}  // namespace cohparams
// C--splitting integral
namespace splitint {
inline F2<double, 1000, 1000> splitiggv, splitiqqv, splitiqgv;
inline F1<double, 1000> qval, zmval;
inline double qmax = 0, zmmin = 0;
inline int npoint = 0;
}  // namespace splitint
// C--pdf common block
namespace pdfs { inline F2<double, 2, 1000> qinqx, ginqx, qingx, gingx; }
// C--cross section common block
namespace xsecs { inline F2<double, 1001, 101> intq1, intq2, intg1, intg2; }
// C--Sudakov common block
namespace insuda { inline F2<double, 1000, 2> sudaqq, sudaqg, sudagg, sudagc; }
// C--exponential integral for negative arguments
namespace expint {
inline F2<double, 3, 1000> eixs;
inline double valmax = 0;
inline int nval = 0;
}  // namespace expint
// C--discard event flag
namespace disc {
inline int ndisc = 0, nstrange = 0, ngood = 0, errcount = 0;
inline double wdisc = 0;
inline bool discard = false;
}  // namespace disc
// C--factor in front of formation times
namespace ftimefac { inline double ftfac = 0; }
// C--factor in front of alphas argument
namespace alphasfac { inline double ptfac = 0; }
// C--factor to scale cross section
namespace scalfac { inline double xsecscalefac = 0; }
// C--number of extrapolations in tables
namespace extrapolations {
inline int ntotspliti = 0, noverspliti = 0, ntotpdf = 0, noverpdf = 0;
inline int ntotxsec = 0, noverxsec = 0, ntotsuda = 0, noversuda = 0;
}  // namespace extrapolations
// C--event weight
namespace weight { inline double evweight = 0, sumofweights = 0; }
// C--identifier of file for hepmc output and logfile (Fortran unit numbers;
// the C++ side maps them onto streams, see io.hpp)
namespace hepmcid { inline int hpmcfid = 0, logfid = 0; }
// C--max rapidity
namespace rapmax { inline double etamax = 0; }
// C--memory for error message from getdeltat
namespace errline { inline int errl = 0; }
// C--total number of generated events (all isospin channels)
namespace totnumevents { inline int nevdone = 0; }
// C--organisation of event record
namespace evrecord {
inline int nsim = 0, npart = 0, offset = 0, hadrotype = 0;
inline double sqrts = 0;
inline std::string collider, channel;  // CHARACTER*4
inline std::string isochannel;         // CHARACTER*2
inline bool hadro = false, shorthepmc = false, isrscat = false;
}  // namespace evrecord
// C--extra storage for scattering centres before interactions
namespace storescatcen {
inline int nscatcen = 0;
inline int maxnscatcen = 22990;  // DATA maxnscatcen/22990/
inline F1<int, kLines> scatflav;
inline F2<double, kLines, 5> scatcen;
inline bool writescatcen = false, writedummies = false, dosubtraction = false;
inline int putonshell = 0;
inline double rsub = 0;
}  // namespace storescatcen
// position (x,y,z,t) of the stored scattering centres (library mode: holes for JETSCAPE)
namespace storescatcen { inline F2<double, kLines, 4> scatpos; }
// JEWEL used as a library (JETSCAPE module) instead of the stand-alone program
namespace lib {
inline bool active = false;
// origin of a line, for lines whose mother pointer COMPRESSEVENT has cleared:
// -1 not set, 0 shower, 1 medium recoil, 2 deleted (descends from a status -1 line)
inline F1<int, kLines> origin;
// line of the recoil that belongs to stored scattering centre j (scatcen(j,.));
// -2: the recoil line was compressed away after it had evolved further
inline F1<int, kLines> scatrec;
// number of stored scattering centres before a coherent scattering (restored when it is undone)
inline int nscatbef = 0;
}  // namespace lib
// C--Pythia parameters
namespace pythiaparams {
inline double ptmin = 0, ptmax = 0;
inline bool weighted = false, mpi = false, pi0dec = false;
}  // namespace pythiaparams
// C--variables for angular ordering
namespace angor {
inline F1<double, kLines> zd, thetad;
inline F1<bool, kLines> qqbard, aoreq;
}  // namespace angor
// C--time of evolution start and medium flag per line
namespace evoltimes {
inline F1<double, kLines> starttime;
inline F1<bool, kLines> medind;
}  // namespace evoltimes
// C--production point and time of partons
namespace time { inline F2<double, kLines, 5> mv; }
// C--colour index
namespace colour {
inline F1<int, kLines> trip, anti;
inline int colmax = 0;
}  // namespace colour
// C--position of the jet production point
namespace jetpoint { inline double x0 = 0, y0 = 0; }
// C--line numbers of special partons
namespace speclines {
inline int lme1 = 0, lme2 = 0, lv = 0, llep1 = 0, llep2 = 0, lme1orig = 0, lme2orig = 0;
inline std::string type1, type2;  // CHARACTER*2
}  // namespace speclines
// C--variables for coherent scattering
namespace coherent {
inline int nstart = 0, nend = 0;
inline F2<double, 10000, 6> allqs;
inline F2<double, 10000, 10> scatcentres;
inline F1<double, 4> qsumvec;
inline double qsum2 = 0;
}  // namespace coherent
// C--information needed to revert a coherent scattering
namespace revertcohscat {
inline double decaytl = 0, decaytsib = 0;
inline int nbef = 0, lbef = 0, lsib = 0, klbef1 = 0, klsib1 = 0;
inline bool didcohscat = false;
}  // namespace revertcohscat
// C--diagnostic histograms
namespace histograms {
inline F1<double, 100> histot, histothat, histotacc, histotrej, histothatacc, histothatrej,
    histotcoh, histothatcoh, histotheta, histoz, histomass, histotau;
inline double tmin = 0, timax = 0, thatmin = 0, thatmax = 0, thetamin = 0, thetamax = 0,
              zmin = 0, zmax = 0, massmin = 0, massmax = 0, taumin = 0, taumax = 0;
inline int ntbins = 0, nthatbins = 0, ntbinsaccrej = 0, nthatbinsaccrej = 0, nthetabins = 0,
           nzbins = 0, nmassbins = 0, ntaubins = 0;
}  // namespace histograms
// C--dummy particles (thermal momenta) for subtraction output
namespace storedummies { inline F2<double, kLines, 5> dummies; }
// C--variables passed to the Sudakov integrand (DERIV)
namespace sudaint {
inline double qa = 0, za2 = 0, eb = 0, t = 0;
inline bool instate = false;
inline std::string typ;  // CHARACTER*2
}  // namespace sudaint
// C--variables passed to the splitting-function integrand
namespace intsplitf { inline double qquad = 0, fm = 0; }
// C--variables passed to the pdf integrand
namespace pdfintv { inline double xmax = 0, z = 0; }
// C--variables passed to the cross-section integrand
namespace xsecv { inline double qlow = 0, mdx = 0; }

}  // namespace jewel::cb
