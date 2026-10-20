// JEWEL 2.6.0 C++ port -- interface to the modified PYTHIA 6.4.25 (Fortran).
//
// PYTHIA 6 remains a Fortran library. JEWEL exchanges data with it through
// COMMON blocks, which are mapped here onto extern "C" structs with identical
// layout (gfortran: COMMON /NAME/ -> symbol name_). Fortran INTEGER is int,
// DOUBLE PRECISION is double.
//
// Within JEWEL code the PYTHIA common block members are used through the
// namespaces jewel::py::<block>, e.g.
//     using namespace jewel::py::pyjets;   // N, NPAD, K, P, V
//     K(i,2) = 21;  P(i,4) += 1.0;
// which mirrors the COMMON declarations of the original Fortran routines.
#pragma once

#include <cstddef>

#include "jewel/farray.hpp"

namespace jewel::f77 {

constexpr int kMaxLines = 23000;  // enlarged PYJETS record of the JEWEL PYTHIA

extern "C" {

struct PyjetsBlock {
  int n;
  int npad;
  FArray<int, Dim<1, kMaxLines>, Dim<1, 5>> k;
  FArray<double, Dim<1, kMaxLines>, Dim<1, 5>> p;
  FArray<double, Dim<1, kMaxLines>, Dim<1, 5>> v;
};
struct Pydat1Block {
  F1<int, 200> mstu;
  F1<double, 200> paru;
  F1<int, 200> mstj;
  F1<double, 200> parj;
};
struct Pydat3Block {
  F2<int, 500, 3> mdcy;
  F2<int, 8000, 2> mdme;
  F1<double, 8000> brat;
  F2<int, 8000, 5> kfdp;
};
struct PysubsBlock {
  int msel;
  int mselpd;
  F1<int, 500> msub;
  FArray<int, Dim<1, 2>, Dim<-40, 40>> kfin;
  F1<double, 200> ckin;
};
struct PyparsBlock {
  F1<int, 200> mstp;
  F1<double, 200> parp;
  F1<int, 200> msti;
  F1<double, 200> pari;
};
struct PydatrBlock {
  F1<int, 6> mrpy;
  F1<double, 100> rrpy;
};
struct PdfasBlock {  // set by the modified PYTHIA from LHAPDF, read by JEWEL
  double pdfalphas;
};
struct WexpoBlock {  // set by JEWEL, read by the customised PYEVWT
  double weightex;
};

extern PyjetsBlock pyjets_;
extern Pydat1Block pydat1_;
extern Pydat3Block pydat3_;
extern PysubsBlock pysubs_;
extern PyparsBlock pypars_;
extern PydatrBlock pydatr_;
extern PdfasBlock pdfas_;
extern WexpoBlock wexpo_;

// PYTHIA routines (gfortran ABI: everything by reference, hidden trailing
// std::size_t string lengths for CHARACTER arguments).
double pyr_(const int* idummy);
double pyp_(const int* i, const int* j);
double pymass_(const int* kf);
int pycomp_(const int* kf);
// three-jet matrix elements, massive (Norrbin & Sjostrand); used for the heavy-flavour e+e- ME correction
double pymael_(const int* ni, const double* x1, const double* x2, const double* r1, const double* r2,
               const double* alpha);
void pyinit_(const char* frame, const char* beam, const char* target, const double* win,
             std::size_t lframe, std::size_t lbeam, std::size_t ltarget);
void pyevnt_();
void pyexec_();
void pylist_(const int* mlist);
void pyedit_(const int* medit);
void pyrget_(const int* lfn, const int* move);
void pyrset_(const int* lfn, const int* move);
void pyrobo_(const int* imi, const int* ima, const double* the, const double* phi,
             const double* bex, const double* bey, const double* bez);

// Small Fortran I/O shim (fortran/fio.f90): lets C++ open/write/close
// Fortran units, needed because PYTHIA writes its messages to unit MSTU(11)
// and PYRGET/PYRSET use unformatted Fortran units.
void jewel_fopen_(const int* unit, const char* path, const int* pathlen, const int* unformatted,
                  const int* status_old, int* ierr);
void jewel_fclose_(const int* unit);
void jewel_fwrite_(const int* unit, const char* text, const int* len);
void jewel_fflush_(const int* unit);

}  // extern "C"

}  // namespace jewel::f77

// ---------------------------------------------------------------------------
// Fortran-style access to the PYTHIA common blocks.
namespace jewel::py {

namespace pyjets {
inline int& N = f77::pyjets_.n;
inline int& NPAD = f77::pyjets_.npad;
inline auto& K = f77::pyjets_.k;
inline auto& P = f77::pyjets_.p;
inline auto& V = f77::pyjets_.v;
}  // namespace pyjets
namespace pydat1 {
inline auto& MSTU = f77::pydat1_.mstu;
inline auto& PARU = f77::pydat1_.paru;
inline auto& MSTJ = f77::pydat1_.mstj;
inline auto& PARJ = f77::pydat1_.parj;
}  // namespace pydat1
namespace pydat3 {
inline auto& MDCY = f77::pydat3_.mdcy;
inline auto& MDME = f77::pydat3_.mdme;
inline auto& BRAT = f77::pydat3_.brat;
inline auto& KFDP = f77::pydat3_.kfdp;
}  // namespace pydat3
namespace pysubs {
inline int& MSEL = f77::pysubs_.msel;
inline int& MSELPD = f77::pysubs_.mselpd;
inline auto& MSUB = f77::pysubs_.msub;
inline auto& KFIN = f77::pysubs_.kfin;
inline auto& CKIN = f77::pysubs_.ckin;
}  // namespace pysubs
namespace pypars {
inline auto& MSTP = f77::pypars_.mstp;
inline auto& PARP = f77::pypars_.parp;
inline auto& MSTI = f77::pypars_.msti;
inline auto& PARI = f77::pypars_.pari;
}  // namespace pypars
namespace pydatr {
inline auto& MRPY = f77::pydatr_.mrpy;
inline auto& RRPY = f77::pydatr_.rrpy;
}  // namespace pydatr
namespace pdfas {
inline double& pdfalphas = f77::pdfas_.pdfalphas;
}
namespace wexpo {
inline double& weightex = f77::wexpo_.weightex;
}

// Thin, by-value wrappers around the PYTHIA functions.
inline double pyr(int idummy = 0) { return f77::pyr_(&idummy); }
inline double pyp(int i, int j) { return f77::pyp_(&i, &j); }
inline double pymass(int kf) { return f77::pymass_(&kf); }
inline int pycomp(int kf) { return f77::pycomp_(&kf); }
inline void pyevnt() { f77::pyevnt_(); }
inline void pyexec() { f77::pyexec_(); }
inline void pylist(int mlist) { f77::pylist_(&mlist); }
inline void pyedit(int medit) { f77::pyedit_(&medit); }
inline void pyrget(int lfn, int move) { f77::pyrget_(&lfn, &move); }
inline void pyrset(int lfn, int move) { f77::pyrset_(&lfn, &move); }
// PYROBO(IMI,IMA,THE,PHI,BEX,BEY,BEZ): rotate/boost lines imi..ima.
inline void pyrobo(int imi, int ima, double the, double phi, double bex, double bey, double bez) {
  f77::pyrobo_(&imi, &ima, &the, &phi, &bex, &bey, &bez);
}
// PYINIT: strings are passed with their exact length (Fortran blank-pads).
void pyinit(const char* frame, const char* beam, const char* target, double win);

}  // namespace jewel::py
