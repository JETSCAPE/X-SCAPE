// JEWEL 2.6.0 C++ port -- active medium model.
#include "jewel/medium.hpp"

#include <cstdio>
#include <cstdlib>

namespace jewel {

namespace {
std::unique_ptr<Medium>& active() {
  static std::unique_ptr<Medium> m;
  return m;
}
}  // namespace

Medium& medium() {
  if (!active()) {
    std::fputs("JEWEL: no medium model selected\n", stderr);
    std::abort();
  }
  return *active();
}
void set_medium(std::unique_ptr<Medium> m) { active() = std::move(m); }

}  // namespace jewel
