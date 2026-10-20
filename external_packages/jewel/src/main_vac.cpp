// JEWEL 2.6.0 C++ port -- executable jewel-2.6.0-vac (medium-vac.f).
#include "jewel/medium.hpp"
#include "jewel/routines.hpp"

int main(int argc, char** argv) {
  jewel::set_medium(jewel::make_medium_vac());
  return jewel::run_jewel(argc, argv);
}
