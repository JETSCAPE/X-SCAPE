// JEWEL 2.6.0 C++ port -- PYTHIA 6 glue.
#include "jewel/pythia6.hpp"

#include <cstring>

namespace jewel::py {

void pyinit(const char* frame, const char* beam, const char* target, double win) {
  f77::pyinit_(frame, beam, target, &win, std::strlen(frame), std::strlen(beam),
               std::strlen(target));
}

}  // namespace jewel::py
