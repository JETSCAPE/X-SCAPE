// JEWEL 2.6.0 C++ port -- Fortran-layout arrays.
//
// FArray reproduces Fortran array semantics exactly: column-major storage,
// arbitrary lower bounds (default 1) and no hidden members, so an FArray
// placed inside an extern "C" struct has the same memory layout as the
// corresponding Fortran COMMON block member.  This lets the C++ code share
// the PYTHIA 6 common blocks with the Fortran library directly.
//
// Access: a(i) / a(i,j) or the C++23 multidimensional subscript a[i,j].
// Define JEWEL_BOUNDS_CHECK to enable bounds assertions.
#pragma once

#include <array>
#include <cassert>
#include <cstddef>

namespace jewel {

template <int Lo, int Hi>
struct Dim {
  static constexpr int lo = Lo;
  static constexpr int hi = Hi;
  static constexpr std::size_t extent = static_cast<std::size_t>(Hi - Lo + 1);
  static_assert(Hi >= Lo);
};

template <typename T, typename... Dims>
class FArray {
  static_assert(sizeof...(Dims) >= 1 && sizeof...(Dims) <= 3);
  static constexpr std::size_t total = (Dims::extent * ...);
  std::array<T, total> d_{};

  static constexpr std::array<int, sizeof...(Dims)> lo_{Dims::lo...};
  static constexpr std::array<int, sizeof...(Dims)> hi_{Dims::hi...};
  static constexpr std::array<std::size_t, sizeof...(Dims)> ext_{Dims::extent...};

  // Signed offset: the original Fortran occasionally reads one element
  // outside an array (e.g. K(0,1) in PYJETS, which aliases NPAD); the same
  // address is produced here rather than an unsigned wrap-around.
  template <typename... I>
  static constexpr std::ptrdiff_t offset(I... idx) noexcept {
    static_assert(sizeof...(I) == sizeof...(Dims), "wrong number of indices");
    const std::array<int, sizeof...(Dims)> ix{static_cast<int>(idx)...};
    std::ptrdiff_t off = 0, stride = 1;
    for (std::size_t r = 0; r < sizeof...(Dims); ++r) {
#ifdef JEWEL_BOUNDS_CHECK
      assert(ix[r] >= lo_[r] && ix[r] <= hi_[r] && "FArray index out of bounds");
#endif
      off += static_cast<std::ptrdiff_t>(ix[r] - lo_[r]) * stride;
      stride *= static_cast<std::ptrdiff_t>(ext_[r]);
    }
    return off;
  }

 public:
  using value_type = T;

  template <typename... I>
  constexpr T& operator()(I... idx) noexcept { return *(d_.data() + offset(idx...)); }
  template <typename... I>
  constexpr const T& operator()(I... idx) const noexcept { return *(d_.data() + offset(idx...)); }

  template <typename... I>
  constexpr T& operator[](I... idx) noexcept { return *(d_.data() + offset(idx...)); }
  template <typename... I>
  constexpr const T& operator[](I... idx) const noexcept { return *(d_.data() + offset(idx...)); }

  constexpr T* data() noexcept { return d_.data(); }
  constexpr const T* data() const noexcept { return d_.data(); }
  static constexpr std::size_t size() noexcept { return total; }
  constexpr void fill(const T& v) { d_.fill(v); }
};

// Shorthands for the common 1-based cases: F1<T,n> ~ T x(n), F2<T,n,m> ~ T x(n,m).
template <typename T, int N1>
using F1 = FArray<T, Dim<1, N1>>;
template <typename T, int N1, int N2>
using F2 = FArray<T, Dim<1, N1>, Dim<1, N2>>;
template <typename T, int N1, int N2, int N3>
using F3 = FArray<T, Dim<1, N1>, Dim<1, N2>, Dim<1, N3>>;

}  // namespace jewel
