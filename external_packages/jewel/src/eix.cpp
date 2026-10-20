// JEWEL 2.6.0 C++ port -- exponential integral Ei(x). Translated from meix.f:1-38.
//
// C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// C++ Copyright S. Zhang and J. M. Jing                               ++
// C++ Computation of Special Functions, John Wiley & Sons, New York,  ++
// C++ 1996 and http://jin.ece.illinois.edu .                          ++
// C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include <cmath>

#include "jewel/fortran_compat.hpp"
#include "jewel/routines.hpp"

namespace jewel {

//
//       ============================================
//       Purpose: Compute exponential integral Ei(x)
//       Input :  x  --- Argument of Ei(x)
//       Output:  EI --- Ei(x) ( x > 0 )
//       ============================================
//
// IMPLICIT DOUBLE PRECISION (A-H,O-Z): X, EI, R, GA are double; K is INTEGER.
void eix(const double& x, double& ei) {
  if (x == 0.0) {
    ei = -1.0e+300;
  } else if (x <= 40.0) {
    ei = 1.0;
    double r = 1.0;
    for (int k = 1; k <= 100; ++k) {
      r = r * static_cast<double>(k) * x / ipow<2>(static_cast<double>(k) + 1.0);
      ei = ei + r;
      if (std::fabs(r / ei) <= 1.0e-15) break;  // GO TO 20
    }
    // 20
    const double ga = 0.5772156649015328;
    ei = ga + std::log(x) + x * ei;
  } else {
    ei = 1.0;
    double r = 1.0;
    for (int k = 1; k <= 20; ++k) {
      r = r * static_cast<double>(k) / x;
      ei = ei + r;
    }
    ei = std::exp(x) / x * ei;
  }
}

}  // namespace jewel
