// JEWEL 2.6.0 C++ port -- helpers reproducing Fortran numeric semantics.
//
// The port aims at results identical to the gfortran-compiled original for
// the same random-number sequence, so the Fortran rules that affect rounding
// are reproduced explicitly here (see PORTING.md for the full list).
#pragma once

#include <cmath>
#include <cstdlib>

namespace jewel {

// Fortran default-REAL (single precision) literal, e.g. `0.4` in Fortran
// source is REAL*4 and becomes 0.4000000059604645 when stored in a DOUBLE
// PRECISION variable. Write r4(0.4f) for such literals.
constexpr double r4(float v) noexcept { return static_cast<double>(v); }

// X**N for a constant integer N, with exactly the multiplication chain that
// gfortran emits (GCC powi_table / powi_as_mults). All exponents occurring in
// JEWEL are covered; larger ones fall back to binary powering (libgcc
// __powidf2), which is what gfortran uses for non-constant exponents.
template <int Nexp>
constexpr double ipow(double x) noexcept {
  if constexpr (Nexp < 0) {
    return 1.0 / ipow<-Nexp>(x);
  } else if constexpr (Nexp == 0) {
    return 1.0;
  } else if constexpr (Nexp == 1) {
    return x;
  } else if constexpr (Nexp == 2) {
    return x * x;
  } else if constexpr (Nexp == 3) {
    const double x2 = x * x;
    return x * x2;
  } else if constexpr (Nexp == 4) {
    const double x2 = x * x;
    return x2 * x2;
  } else if constexpr (Nexp == 5) {
    const double x2 = x * x, x3 = x * x2;
    return x2 * x3;
  } else if constexpr (Nexp == 6) {
    const double x2 = x * x, x3 = x * x2;
    return x3 * x3;
  } else if constexpr (Nexp == 8) {
    const double x2 = x * x, x4 = x2 * x2;
    return x4 * x4;
  } else {
    double r = (Nexp & 1) ? x : 1.0;
    int n = Nexp;
    while (n >>= 1) {
      x = x * x;
      if (n & 1) r = r * x;
    }
    return r;
  }
}
// INTEGER**INTEGER
template <int Nexp>
constexpr int ipowi(int x) noexcept {
  static_assert(Nexp >= 0);
  int r = 1;
  for (int i = 0; i < Nexp; ++i) r *= x;
  return r;
}
// X**N with a run-time integer N (libgcc __powidf2).
inline double powi(double x, int n) noexcept {
  unsigned m = n < 0 ? -static_cast<unsigned>(n) : static_cast<unsigned>(n);
  double r = (m % 2) ? x : 1.0;
  while (m >>= 1) {
    x = x * x;
    if (m % 2) r = r * x;
  }
  return n < 0 ? 1.0 / r : r;
}

// gfortran merges SIN(X) and COS(X) of the same argument within a routine
// into one call cexp(0+iX) (GCC's cexpi, macOS has no sincos). Where the
// original is compiled that way, the port uses these (see PORTING.md);
// isolated SIN/COS stay std::sin/std::cos.
#if defined(__APPLE__)
extern "C" _Complex double jewel_libm_cexp(_Complex double) noexcept __asm__("_cexp");
inline _Complex double gf_cexpi(double x) noexcept {
  _Complex double z;
  __real__ z = 0.0;
  __imag__ z = x;
  return jewel_libm_cexp(z);
}
#else
// glibc: gfortran lowers cexpi to sincos(), so call that for the same bits.
extern "C" void jewel_libm_sincos(double, double*, double*) noexcept __asm__("sincos");
inline _Complex double gf_cexpi(double x) noexcept {
  double s, c;
  jewel_libm_sincos(x, &s, &c);
  _Complex double z;
  __real__ z = c;
  __imag__ z = s;
  return z;
}
#endif
inline double gf_cos(double x) noexcept { return __real__ gf_cexpi(x); }
inline double gf_sin(double x) noexcept { return __imag__ gf_cexpi(x); }

// SIGN(A,B)
inline double fsign(double a, double b) noexcept { return std::copysign(std::fabs(a), b); }
inline int fsign(int a, int b) noexcept { return b >= 0 ? std::abs(a) : -std::abs(a); }
// NINT(X): round half away from zero
inline int nint(double x) noexcept { return static_cast<int>(std::lround(x)); }
// INT(X): truncate towards zero
inline int fint(double x) noexcept { return static_cast<int>(x); }
// MOD for reals (A - INT(A/P)*P, same as fmod)
inline double fmod_f(double a, double p) noexcept { return std::fmod(a, p); }

}  // namespace jewel
