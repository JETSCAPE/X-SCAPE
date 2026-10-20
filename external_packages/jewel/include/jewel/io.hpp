// JEWEL 2.6.0 C++ port -- output channels.
//
// The log file stays a Fortran unit (cb::hepmcid::logfid) because PYTHIA
// writes into it via MSTU(11); JEWEL messages are routed through the small
// Fortran shim so that both appear in their original order. The HepMC file
// is written directly from C++.
#pragma once

#include <cstdio>
#include "jewel/jfmt.hpp"
#include <string>
#include <string_view>

namespace jewel::io {

// Open a Fortran unit (formatted or unformatted). Returns false on failure.
bool fortran_open(int unit, std::string_view path, bool unformatted = false, bool must_exist = false);
void fortran_close(int unit);
void fortran_flush(int unit);
void fortran_write(int unit, std::string_view line);

// One line in the JEWEL log (Fortran `write(logfid,*) ...`; list-directed
// output starts with a blank, which callers include in the format string).
void log_line(std::string_view line);
template <typename... Args>
void log(jewel::jfmt::format_string<Args...> fmt, Args&&... args) {
  log_line(jewel::jfmt::format(fmt, std::forward<Args>(args)...));
}
void log_flush();

// Format a DOUBLE PRECISION value the way gfortran list-directed output does
// (e.g. "0.11800000071525574"); used where log output is compared by eye.
std::string flist(double x);

// The HepMC event file (unit HPMCFID in the Fortran code).
std::FILE*& hepmc();

}  // namespace jewel::io
