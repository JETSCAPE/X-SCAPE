// JEWEL 2.6.0 C++ port -- output channels (see io.hpp).
#include "jewel/io.hpp"

#include <charconv>
#include <cstdio>

#include "jewel/commons.hpp"
#include "jewel/pythia6.hpp"

namespace jewel::io {

bool fortran_open(int unit, std::string_view path, bool unformatted, bool must_exist) {
  const int len = static_cast<int>(path.size());
  const int unf = unformatted ? 1 : 0, old = must_exist ? 1 : 0;
  int ierr = 0;
  f77::jewel_fopen_(&unit, path.data(), &len, &unf, &old, &ierr);
  return ierr == 0;
}
void fortran_close(int unit) { f77::jewel_fclose_(&unit); }
void fortran_flush(int unit) { f77::jewel_fflush_(&unit); }
void fortran_write(int unit, std::string_view line) {
  const int len = static_cast<int>(line.size());
  f77::jewel_fwrite_(&unit, line.data(), &len);
}

void log_line(std::string_view line) { fortran_write(cb::hepmcid::logfid, line); }
void log_flush() { fortran_flush(cb::hepmcid::logfid); }

std::string flist(double x) {
  // gfortran list-directed REAL(8): 17 significant digits.
  return jewel::jfmt::format("{:.17g}", x);
}

std::FILE*& hepmc() {
  static std::FILE* f = nullptr;
  return f;
}

}  // namespace jewel::io
