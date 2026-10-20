// JEWEL 2.6.0 C++ port -- executable jewel-2.6.0-simple (medium-simple.f).
#include "jewel/medium.hpp"
#include "jewel/routines.hpp"

int main(int argc, char** argv) {
  jewel::set_medium(jewel::make_medium_simple());
  return jewel::run_jewel(argc, argv);
}
