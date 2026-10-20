// JEWEL 2.6.0 C++ port -- main program, initialisation and event generation driver
// (PROGRAM JEWEL, INIT, INITPYTHIA, GENEVENT). Translated from jewel-2.6.0.f:65-1400.
#include <algorithm>
#include <charconv>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <filesystem>
#include "jewel/jfmt.hpp"
#include <fstream>
#include "jewel/jfmt.hpp"
#include <sstream>
#include <string>
#include <string_view>
#include <system_error>
#include <vector>

#include "jewel/heavyflavour.hpp"
#include "jewel/commons.hpp"
#include "jewel/fortran_compat.hpp"
#include "jewel/io.hpp"
#include "jewel/medium.hpp"
#include "jewel/pythia6.hpp"
#include "jewel/routines.hpp"

namespace jewel {

namespace {

// Command line (replaces IARGC()/GETARG()); set by run_jewel before init().
int g_argc = 0;
char** g_argv = nullptr;

// ---------------------------------------------------------------------------
// Helpers emulating the Fortran I/O used in init().

// Fortran CHARACTER*n assignment: truncate/pad to n characters.
std::string fchar(std::string_view s, std::size_t n) {
  std::string r(s.substr(0, std::min(n, s.size())));
  r.resize(n, ' ');
  return r;
}
// TRIM()
std::string ftrim(std::string_view s) {
  const auto e = s.find_last_not_of(' ');
  return e == std::string_view::npos ? std::string() : std::string(s.substr(0, e + 1));
}
// Value as it appears in a list-directed WRITE of a LOGICAL.
const char* flog(bool b) { return b ? "T" : "F"; }

bool is_lsep(char c) { return c == ' ' || c == '\t' || c == '\r' || c == '\n' || c == ',' || c == '/'; }

// First item of a list-directed READ(value,*): returns false if no value is
// available (end of record, null value or '/'), in which case Fortran leaves
// the variable unchanged (iostat /= 0 or null value).
bool list_token(std::string_view s, std::string& tok) {
  std::size_t i = 0;
  while (i < s.size() && (s[i] == ' ' || s[i] == '\t' || s[i] == '\r' || s[i] == '\n')) ++i;
  if (i >= s.size() || s[i] == ',' || s[i] == '/') return false;
  if (s[i] == '\'' || s[i] == '"') {  // quoted character constant
    const char q = s[i++];
    tok.clear();
    while (i < s.size()) {
      if (s[i] == q) {
        if (i + 1 < s.size() && s[i + 1] == q) {
          tok += q;
          i += 2;
          continue;
        }
        return true;
      }
      tok += s[i++];
    }
    return true;
  }
  const std::size_t b = i;
  while (i < s.size() && !is_lsep(s[i])) ++i;
  tok.assign(s.substr(b, i - b));
  return true;
}

bool parse_int(const std::string& tok, int& v) {
  std::string_view t = tok;
  if (!t.empty() && t[0] == '+') t.remove_prefix(1);
  if (t.empty()) return false;
  int x = 0;
  const auto [p, ec] = std::from_chars(t.data(), t.data() + t.size(), x);
  if (ec != std::errc() || p != t.data() + t.size()) return false;
  v = x;
  return true;
}

bool parse_real(const std::string& tok, double& v) {
  std::string t = tok;
  if (t.empty()) return false;
  for (char& c : t) {
    if (c == 'x' || c == 'X') return false;  // no hex floats in Fortran
    if (c == 'd' || c == 'D' || c == 'q' || c == 'Q') c = 'e';
  }
  char* end = nullptr;
  const double x = std::strtod(t.c_str(), &end);
  if (end != t.c_str() + t.size()) return false;
  v = x;
  return true;
}

// READ(value,*,iostat=ios) <integer>
void read_list(std::string_view s, int& v) {
  std::string tok;
  if (list_token(s, tok)) parse_int(tok, v);
}
// READ(value,*,iostat=ios) <double precision>
void read_list(std::string_view s, double& v) {
  std::string tok;
  if (list_token(s, tok)) parse_real(tok, v);
}
// READ(value,*,iostat=ios) <logical>: optional '.', then T/t or F/f.
void read_list(std::string_view s, bool& v) {
  std::string tok;
  if (!list_token(s, tok) || tok.empty()) return;
  std::size_t i = 0;
  if (tok[0] == '.') i = 1;
  if (i >= tok.size()) return;
  if (tok[i] == 'T' || tok[i] == 't') v = true;
  else if (tok[i] == 'F' || tok[i] == 'f') v = false;
}
// READ(value,*,iostat=ios) <character*n>
void read_list(std::string_view s, std::string& v, std::size_t n) {
  std::string tok;
  if (list_token(s, tok)) v = ftrim(fchar(tok, n));
}
// READ(value,'(a)',iostat=ios) <character>: the whole value; trailing blanks are
// irrelevant (Fortran blank-pads, OPEN trims them), leading blanks are kept.
void read_a(std::string_view s, std::string& v) { v = ftrim(s); }

// Table files: list-directed READ(10,*) of n items. Every READ starts a new
// record; items are taken from as many records as needed, the rest of the last
// record is skipped.
std::vector<std::string> table_record(std::istream& in, std::size_t n, const std::string& file) {
  std::vector<std::string> toks;
  std::string line;
  while (toks.size() < n) {
    if (!std::getline(in, line)) {
      // NOTE(port): Fortran has no iostat here; gfortran aborts with a run-time error.
      jewel::jfmt::println(stderr, "Fortran runtime error: End of file while reading {}", file);
      std::exit(2);
    }
    std::istringstream ls(line);
    std::string t;
    while (toks.size() < n && ls >> t) {
      // list-directed separators: blanks and commas
      std::size_t b = 0;
      while (b <= t.size()) {
        const std::size_t e = t.find(',', b);
        const std::string part = t.substr(b, e == std::string::npos ? std::string::npos : e - b);
        if (!part.empty() && toks.size() < n) toks.push_back(part);
        if (e == std::string::npos) break;
        b = e + 1;
      }
    }
  }
  return toks;
}
[[noreturn]] void table_error(const std::string& file) {
  jewel::jfmt::println(stderr, "Fortran runtime error: Bad value during read of {}", file);
  std::exit(2);
}
double table_real(const std::string& tok, const std::string& file) {
  double v = 0;
  if (!parse_real(tok, v)) table_error(file);
  return v;
}
int table_int(const std::string& tok, const std::string& file) {
  int v = 0;
  if (!parse_int(tok, v)) table_error(file);
  return v;
}

// INQUIRE(file=...,exist=...)
bool file_exists(const std::string& f) {
  std::error_code ec;
  return std::filesystem::exists(f, ec);
}

}  // namespace

// ===========================================================================
// PROGRAM JEWEL
int run_jewel(int argc, char** argv) {
  using namespace py::pyjets;
  using namespace py::pydat1;
  using namespace py::pydat3;
  using namespace py::pysubs;
  using namespace py::pypars;
  using namespace py::pydatr;
  using namespace cb::hepmcid;
  using namespace cb::np;
  using namespace cb::evrecord;
  using namespace cb::disc;
  using namespace cb::weight;
  using namespace cb::extrapolations;

  // C--local variables
  int nsimpp = 0, nsimpn = 0, nsimnp = 0, nsimnn = 0, nsimsum = 0, nsimchn = 0;
  double sumofweightstot, wdisctot, scalefac;
  std::string b1, b2;

  g_argc = argc;
  g_argv = argv;

  init();

  sumofweightstot = 0.0;
  wdisctot = 0.0;

  // C--e+ + e- event generation
  if (collider == "EEJJ") {
    b1 = "e+";
    b2 = "e-";
    io::log("");
    io::log(" ####################################################");
    io::log("");
    io::log(" generating {} events in {} + {} channel", nsim, b1, b2);
    io::log("");
    io::log(" ####################################################");
    io::log("");
    sumofweights = 0.0;
    wdisc = 0.0;
    initpythia(b1, b2);
    io::log("");
    // C--e+ + e- event loop
    for (int j = 1, j_end = nsim; j <= j_end; ++j) {
      genevent(j, b1, b2);
    }
    sumofweightstot = sumofweightstot + sumofweights;
    wdisctot = wdisctot + wdisc;
    io::log("");
    io::log(" cross section in e+ + e- channel: {} mb", PARI(1));
    io::log(" sum of event weights in e+ + e- channel: {}", sumofweights - wdisc);
    io::log("");

  } else {
    // C--hadronic event generation
    if (isochannel == "PP") {
      nsimpp = nsim;
      nsimpn = 0;
      nsimnp = 0;
      nsimnn = 0;
    } else if (isochannel == "PN") {
      nsimpp = 0;
      nsimpn = nsim;
      nsimnp = 0;
      nsimnn = 0;
    } else if (isochannel == "NP") {
      nsimpp = 0;
      nsimpn = 0;
      nsimnp = nsim;
      nsimnn = 0;
    } else if (isochannel == "NN") {
      nsimpp = 0;
      nsimpn = 0;
      nsimnp = 0;
      nsimnn = nsim;
    } else {
      nsimpp = poissonian(1.0 * nsim * ipowi<2>(nproton) / ipowi<2>(mass));
      nsimpn = poissonian(1.0 * nsim * nproton * (mass - nproton) / ipowi<2>(mass));
      nsimnp = poissonian(1.0 * nsim * nproton * (mass - nproton) / ipowi<2>(mass));
      nsimnn = poissonian(1.0 * nsim * ipowi<2>(mass - nproton) / ipowi<2>(mass));
      nsimsum = nsimpp + nsimpn + nsimnp + nsimnn;
      scalefac = nsim * 1.0 / (nsimsum * 1.0);
      nsimpp = fint(nsimpp * scalefac);
      nsimpn = fint(nsimpn * scalefac);
      nsimnp = fint(nsimnp * scalefac);
      nsimnn = fint(nsimnn * scalefac);
      nsimsum = nsimpp + nsimpn + nsimnp + nsimnn;
    }
    // C--loop over channels
    for (int kk = 1; kk <= 4; ++kk) {
      if (kk == 1) {
        b1 = "p+";
        b2 = "p+";
        nsimchn = nsimpp;
      } else if (kk == 2) {
        b1 = "p+";
        b2 = "n0";
        nsimchn = nsimpn;
      } else if (kk == 3) {
        b1 = "n0";
        b2 = "p+";
        nsimchn = nsimnp;
      } else {
        b1 = "n0";
        b2 = "n0";
        nsimchn = nsimnn;
      }
      io::log("");
      io::log(" ####################################################");
      io::log("");
      io::log(" generating {} events in {} + {} channel", nsimchn, b1, b2);
      io::log("");
      io::log(" ####################################################");
      io::log("");
      sumofweights = 0.0;
      wdisc = 0.0;
      initpythia(b1, b2);
      io::log("");
      // C--event loop
      for (int j = 1, j_end = nsimchn; j <= j_end; ++j) {
        genevent(j, b1, b2);
      }
      sumofweightstot = sumofweightstot + sumofweights;
      wdisctot = wdisctot + wdisc;
      io::log("");
      io::log(" cross section in {} + {} channel: {} mb", b1, b2, PARI(1));
      io::log(" sum of event weights in {} + {} channel: {}", b1, b2, sumofweights - wdisc);
      io::log("");
    }
  }

  // C--finish
  std::fprintf(io::hepmc(), "HepMC::IO_GenEvent-END_EVENT_LISTING\n");
  std::fprintf(io::hepmc(), "\n");  // gfortran: empty record
  std::fclose(io::hepmc());
  io::hepmc() = nullptr;

  io::log("");
  io::log(" number of extrapolations in splitting integral: {} ({}%)", noverspliti,
          (noverspliti * 1.0) / (ntotspliti * 1.0));
  io::log(" number of extrapolations in splitting partonic PDFs: {} ({}%)", noverpdf,
          (noverpdf * 1.0) / (ntotpdf * 1.0));
  io::log(" number of extrapolations in splitting cross sections: {} ({}%)", noverxsec,
          (noverxsec * 1.0) / (ntotxsec * 1.0));
  io::log(" number of extrapolations in Sudakov form factor: {} ({}%)", noversuda,
          (noversuda * 1.0) / (ntotsuda * 1.0));
  io::log("");
  io::log(" total number of discarded events: {}", ndisc);
  io::log(" number of events for which conversion to hepmc failed: {}", nstrange);
  printtime();

  io::fortran_close(logfid);

  // C--write random number generator state to file and close it
  py::pyrget(2, -1);
  io::fortran_close(2);

  return 0;
}

// ***********************************************************************
// ***	  subroutine init
// ***********************************************************************
void init_library(const std::string& paramfile) {
  static std::string prog = "jewel", par;
  static char* argv[3] = {nullptr, nullptr, nullptr};
  par = paramfile;
  argv[0] = prog.data();
  argv[1] = par.data();
  g_argc = par.empty() ? 1 : 2;
  g_argv = argv;
  cb::lib::active = true;  // no HepMC file, no random-number state files
  init();
}

void init() {
  using namespace py::pyjets;
  using namespace py::pydat1;
  using namespace py::pydat3;
  using namespace py::pysubs;
  using namespace py::pypars;
  using namespace py::pydatr;
  using namespace py::pdfas;
  using namespace py::wexpo;
  using namespace cb::pdf;
  using namespace cb::np;
  using namespace cb::param;
  using namespace cb::cohparams;
  using namespace cb::splitint;
  using namespace cb::pdfs;
  using namespace cb::xsecs;
  using namespace cb::insuda;
  using namespace cb::expint;
  using namespace cb::disc;
  using namespace cb::ftimefac;
  using namespace cb::alphasfac;
  using namespace cb::scalfac;
  using namespace cb::extrapolations;
  using namespace cb::weight;
  using namespace cb::hepmcid;
  using namespace cb::rapmax;
  using namespace cb::errline;
  using namespace cb::totnumevents;
  using namespace cb::evrecord;
  using namespace cb::storescatcen;
  using namespace cb::pythiaparams;

  // C--Variables local to this program
  int njob, intmass;
  double eovest;
  [[maybe_unused]] double r;
  std::string pdffile, xsecfile, filemed, filesplit;  // CHARACTER*80
  std::string hepmcfile, logfile;                     // CHARACTER*120
  std::string filename;                               // CHARACTER(LEN=100)
  bool pdfexist, splitiexist, xsecexist;
  bool readran;

  // (data maxnscatcen/22990/ : initialised in commons.hpp)

  hpmcfid = 4;
  logfid = 3;

  // C--default settings
  nsim = 10000;
  njob = 0;
  logfile = "out.log";
  hepmcfile = "out.hepmc";
  filesplit = "splitint.dat";
  pdffile = "pdfs.dat";
  xsecfile = "xsecs.dat";
  filemed = "medium-params.dat";
  nf = 3;
  lqcd = r4(0.4f);
  lps = r4(0.4f);
  q0 = 1.5;
  ptmin = 5.0;
  ptmax = 350.0;
  etamax = r4(3.1f);
  collider = "PPJJ";
  isochannel = "XX";
  channel = "MUON";
  sqrts = 2760;
  pdfset = 13100;
  pdfalphas = r4(0.118f);
  xsecscalefac = 1.0;
  mass = 208;
  nproton = 82;
  weighted = true;
  weightex = 5.0;
  mpi = false;
  mpifsr = true;
  angord = 1;
  aofirstsplit = false;
  allhad = false;
  hadro = true;
  hadrotype = 0;
  pi0dec = true;
  isrscat = true;
  shorthepmc = true;
  compress = true;
  writescatcen = false;
  writedummies = false;
  dosubtraction = false;
  rsub = 1.0;
  putonshell = 0;
  scatrecoil = false;
  recsoftcut = 0.0;
  rechardcut = 5.0;
  kinmode = 1;
  recmode = 0;
  cohlengthfac = 1.0;
  pcohrej = 1.0;
  cohscat = true;
  readran = false;

  if (g_argc - 1 == 0) {
    jewel::jfmt::println(" No parameter file given, will run with default settings.");
  } else {
    filename = fchar(g_argv[1], 100);
    jewel::jfmt::println(" Reading parameters from {}", filename);
    std::ifstream pfile(ftrim(filename));
    if (!pfile) goto l110;
    for (int i = 1; i <= 1000; ++i) {
      std::string line;
      if (!std::getline(pfile, line)) goto l130;
      const std::string buffer = fchar(line, 80);  // CHARACTER*80 buffer
      const char firstchar = buffer[0];
      if (firstchar == '#') continue;
      // pos=scan(buffer,' ') (1-based, 0 if no blank)
      const std::size_t bpos = buffer.find(' ');
      const std::size_t pos = bpos == std::string::npos ? 0 : bpos + 1;
      const std::string label = ftrim(buffer.substr(0, pos));
      const std::string value = buffer.substr(pos);
      if (label == "NEVENT") {
        read_list(value, nsim);
      } else if (label == "NJOB") {
        read_list(value, njob);
      } else if (label == "LOGFILE") {
        read_a(value, logfile);
      } else if (label == "HEPMCFILE") {
        read_a(value, hepmcfile);
      } else if (label == "SPLITINTFILE") {
        read_a(value, filesplit);
      } else if (label == "PDFFILE") {
        read_a(value, pdffile);
      } else if (label == "XSECFILE") {
        read_a(value, xsecfile);
      } else if (label == "MEDIUMPARAMS") {
        read_a(value, filemed);
      } else if (label == "HEAVYFLAVOUR") {
        read_list(value, cb::hf::on);
      } else if (label == "HFSTABLE") {
        read_list(value, cb::hf::stable);
      } else if (label == "HQDS") {
        read_list(value, cb::hf::hqds);
      } else if (label == "NF") {
        read_list(value, nf);
      } else if (label == "LAMBDAQCD") {
        read_list(value, lqcd);
      } else if (label == "LAMBDAQCDPS") {
        read_list(value, lps);
      } else if (label == "Q0") {
        read_list(value, q0);
      } else if (label == "PTMIN") {
        read_list(value, ptmin);
      } else if (label == "PTMAX") {
        read_list(value, ptmax);
      } else if (label == "ETAMAX") {
        read_list(value, etamax);
      } else if (label == "PROCESS") {
        read_list(value, collider, 4);
        // restricted flavour-creation processes (heavy-flavour studies): run as PPJJ with fewer subprocesses
        if ((collider == "PPQQ") || (collider == "PPCC") || (collider == "PPBB")) {
          cb::hf::hardproc = collider.substr(2, 2);
          collider = "PPJJ";
        }
      } else if (label == "ISOCHANNEL") {
        read_list(value, isochannel, 2);
      } else if (label == "CHANNEL") {
        read_list(value, channel, 4);
      } else if (label == "SQRTS") {
        read_list(value, sqrts);
      } else if (label == "PDFSET") {
        read_list(value, pdfset);
      } else if (label == "PDFALPHAS") {
        read_list(value, pdfalphas);
      } else if (label == "XSECSCALEFAC") {
        read_list(value, xsecscalefac);
      } else if (label == "MASS") {
        read_list(value, mass);
      } else if (label == "NPROTON") {
        read_list(value, nproton);
      } else if (label == "WEIGHTED") {
        read_list(value, weighted);
      } else if (label == "WEXPO") {
        read_list(value, weightex);
      } else if (label == "MPI") {
        read_list(value, mpi);
      } else if (label == "MPIFSR") {
        read_list(value, mpifsr);
      } else if (label == "ANGORD") {
        read_list(value, angord);
      } else if (label == "AOFIRSTSPLIT") {
        read_list(value, aofirstsplit);
      } else if (label == "KEEPRECOILS") {
        read_list(value, allhad);
      } else if (label == "SCATRECOIL") {
        read_list(value, scatrecoil);
      } else if (label == "HADRO") {
        read_list(value, hadro);
      } else if (label == "HADROTYPE") {
        read_list(value, hadrotype);
      } else if (label == "PI0DEC") {
        read_list(value, pi0dec);
      } else if (label == "ISRSCAT") {
        read_list(value, isrscat);
      } else if (label == "SHORTHEPMC") {
        read_list(value, shorthepmc);
      } else if (label == "COMPRESS") {
        read_list(value, compress);
      } else if (label == "WRITESCATCEN") {
        read_list(value, writescatcen);
      } else if (label == "WRITEDUMMIES") {
        read_list(value, writedummies);
      } else if (label == "DOSUBTRACTION") {
        read_list(value, dosubtraction);
      } else if (label == "RSUB") {
        read_list(value, rsub);
      } else if (label == "PUTONSHELL") {
        read_list(value, putonshell);
      } else if (label == "RECSOFTCUT") {
        read_list(value, recsoftcut);
      } else if (label == "RECHARDCUT") {
        read_list(value, rechardcut);
      } else if (label == "KINMODE") {
        read_list(value, kinmode);
      } else if (label == "RECMODE") {
        read_list(value, recmode);
      } else if (label == "COHLENGTHFAC") {
        read_list(value, cohlengthfac);
      } else if (label == "PCOHREJ") {
        read_list(value, pcohrej);
      } else if (label == "COHSCAT") {
        read_list(value, cohscat);
      } else if (label == "READRAN") {
        read_list(value, readran);
      } else {
        jewel::jfmt::println(" unknown label {}", label);
      }
    }

    // NOTE(port): as in the Fortran, a parameter file with 1000 or more lines
    // falls through to label 110 and terminates the run.
  l110:
    jewel::jfmt::println(" Unable to open parameter file, will exit the run.");
    std::exit(1);

  l130:
    pfile.close();
    jewel::jfmt::println(" ...done");
  }

  if (recmode == 2) {
    allhad = false;
    scatrecoil = false;
  }

  scalefacm = 1.0;
  ptfac = 1.0;
  ftfac = 1.0;

  if (ptmin < 3.0) ptmin = 3.0;
  if (!writescatcen) writedummies = false;
  if (dosubtraction) {
    writescatcen = false;
    writedummies = false;
  }

  if (!io::fortran_open(logfid, logfile)) {
    // NOTE(port): the Fortran OPEN has no err= branch; gfortran aborts.
    jewel::jfmt::println(stderr, "Fortran runtime error: cannot open file '{}'", logfile);
    std::exit(2);
  }
  MSTU(11) = logfid;

  printtime();
  printlogo(logfid);

  io::log("");
  io::log(" parameters of the run:");
  io::log(" NEVENT        = {}", nsim);
  io::log(" NJOB          = {}", njob);
  io::log(" LOGFILE       = {}", logfile);
  io::log(" HEPMCFILE     = {}", hepmcfile);
  io::log(" SPLITINTFILE  = {}", filesplit);
  io::log(" PDFFILE       = {}", pdffile);
  io::log(" XSECFILE      = {}", xsecfile);
  io::log(" MEDIUMPARAMS  = {}", filemed);
  io::log(" NF            = {}", nf);
  io::log(" HEAVYFLAVOUR  = {}", flog(cb::hf::on));
  if (!cb::hf::hardproc.empty())
    io::log(" hard process restricted to flavour creation: PP{} (run as PPJJ)", cb::hf::hardproc);
  io::log(" HFSTABLE      = {}", flog(cb::hf::stable));
  if (cb::hf::hqds > 0.0) io::log(" HQDS          = {}", cb::hf::hqds);
  io::log(" LAMBDAQCDPS   = {}", lps);
  io::log(" LAMBDAQCD     = {}", lqcd);
  io::log(" Q0            = {}", q0);
  io::log(" PTMIN         = {}", ptmin);
  io::log(" PTMAX         = {}", ptmax);
  io::log(" ETAMAX        = {}", etamax);
  io::log(" PROCESS       = {}", collider);
  io::log(" ISOCHANNEL    = {}", isochannel);
  io::log(" CHANNEL       = {}", channel);
  io::log(" SQRTS         = {}", sqrts);
  io::log(" PDFSET        = {}", pdfset);
  io::log(" PDFALPHAS     = {}", pdfalphas);
  io::log(" XSECSCALEFAC  = {}", xsecscalefac);
  io::log(" MASS          = {}", mass);
  io::log(" NPROTON       = {}", nproton);
  io::log(" WEIGHTED      = {}", flog(weighted));
  io::log(" WEXPO         = {}", weightex);
  io::log(" MPI           = {}", flog(mpi));
  io::log(" MPIFSR        = {}", flog(mpifsr));
  io::log(" ANGORD        = {}", angord);
  io::log(" AOFIRSTSPLIT  = {}", flog(aofirstsplit));
  io::log(" HADRO         = {}", flog(hadro));
  io::log(" HADROTYPE     = {}", hadrotype);
  io::log(" PI0DEC        = {}", flog(pi0dec));
  io::log(" ISRSCAT       = {}", flog(isrscat));
  io::log(" SHORTHEPMC    = {}", flog(shorthepmc));
  io::log(" COMPRESS      = {}", flog(compress));
  io::log(" KEEPRECOILS   = {}", flog(allhad));
  io::log(" SCATRECOIL    = {}", flog(scatrecoil));
  io::log(" RECSOFTCUT    = {}", recsoftcut);
  io::log(" RECHARDCUT    = {}", rechardcut);
  io::log(" WRITESCATCEN  = {}", flog(writescatcen));
  io::log(" WRITEDUMMIES  = {}", flog(writedummies));
  io::log(" DOSUBTRACTION = {}", flog(dosubtraction));
  io::log(" RSUB          = {}", rsub);
  io::log(" PUTONSHELL    = {}", putonshell);
  io::log(" KINMODE       = {}", kinmode);
  io::log(" RECMODE       = {}", recmode);
  io::log(" COHLENGTHFAC  = {}", cohlengthfac);
  io::log(" PCOHREJ       = {}", pcohrej);
  io::log(" COHSCAT       = {}", flog(cohscat));
  io::log("");
  io::log_flush();

  if ((collider != "PPJJ") && (collider != "EEJJ") && (collider != "PPYJ") &&
      (collider != "PPYQ") && (collider != "PPYG") && (collider != "PPZJ") &&
      (collider != "PPZQ") && (collider != "PPZG") && (collider != "PPWJ") &&
      (collider != "PPWQ") && (collider != "PPWG") && (collider != "PPDY")) {
    io::log(" Fatal error: colliding system unknown, will exit now");
    std::exit(1);
  }

  // C--initialize medium
  intmass = mass;  // int(mass)
  medium().medinit(filemed, logfid, etamax, intmass);
  medium().mednextevt();

  // heavy flavour: shower masses = PYTHIA's (also used in the hard process and hadronisation)
  if (cb::hf::on) {
    cb::hf::mc = py::pymass(4);
    cb::hf::mb = py::pymass(5);
    // testing only: scale the shower masses (e.g. 0.01 to check the m -> 0 limit of the heavy-quark code path)
    if (const char* sc = std::getenv("JEWEL_HFMASSSCALE")) {
      cb::hf::mc *= std::atof(sc);
      cb::hf::mb *= std::atof(sc);
      io::log(" JEWEL_HFMASSSCALE = {} (testing)", sc);
    }
    io::log(" heavy-flavour shower: m_c = {} GeV, m_b = {} GeV", cb::hf::mc, cb::hf::mb);
  }

  // library mode (JETSCAPE): no HepMC file, the caller collects the partons
  if (!cb::lib::active) {
    io::hepmc() = std::fopen(hepmcfile.c_str(), "w");
    if (io::hepmc() == nullptr) {
      // NOTE(port): the Fortran OPEN has no err= branch; gfortran aborts.
      jewel::jfmt::println(stderr, "Fortran runtime error: cannot open file '{}'", hepmcfile);
      std::exit(2);
    }
    std::fprintf(io::hepmc(), "\n");  // gfortran: empty record
    std::fprintf(io::hepmc(), "HepMC::Version 2.06.05\n");
    std::fprintf(io::hepmc(), "HepMC::IO_GenEvent-START_EVENT_LISTING\n");
  }

  npart = 2;

  if (ptmax > 0.0) {
    eovest = std::min(1.5 * (ptmax + 50.0) * std::cosh(etamax), sqrts / 2.0);
  } else {
    eovest = sqrts / 2.0;
  }

  eixint();
  insudaint(eovest);

  io::log("");
  splitiexist = file_exists(filesplit);
  if (splitiexist) {
    io::log(" read splitting integrals from {}", filesplit);
    std::ifstream f10(filesplit);
    {
      const auto t = table_record(f10, 3, filesplit);
      qmax = table_real(t[0], filesplit);
      zmmin = table_real(t[1], filesplit);
      npoint = table_int(t[2], filesplit);
    }
    for (int i = 1, i_end = npoint + 1; i <= i_end; ++i) {
      const auto t = table_record(f10, 2, filesplit);
      qval(i) = table_real(t[0], filesplit);
      zmval(i) = table_real(t[1], filesplit);
    }
    for (int i = 1, i_end = npoint + 1; i <= i_end; ++i) {
      for (int j = 1, j_end = npoint + 1; j <= j_end; ++j) {
        const auto t = table_record(f10, 3, filesplit);
        splitiggv(i, j) = table_real(t[0], filesplit);
        splitiqqv(i, j) = table_real(t[1], filesplit);
        splitiqgv(i, j) = table_real(t[2], filesplit);
      }
    }
    f10.close();
  } else {
    io::log(" have to integrate splitting functions, this may take some time");
    splitfncint(eovest);
    splitiexist = file_exists(filesplit);
    if (!splitiexist) {
      io::log(" write splitting integrals to {}", filesplit);
      std::ofstream f10(filesplit);
      f10 << jewel::jfmt::format(" {:.17e} {:.17e} {}\n", qmax, zmmin, npoint);
      for (int i = 1, i_end = npoint + 1; i <= i_end; ++i) {
        f10 << jewel::jfmt::format(" {:.17e} {:.17e}\n", qval(i), zmval(i));
      }
      for (int i = 1, i_end = npoint + 1; i <= i_end; ++i) {
        for (int j = 1, j_end = npoint + 1; j <= j_end; ++j) {
          f10 << jewel::jfmt::format(" {:.17e} {:.17e} {:.17e}\n", splitiggv(i, j), splitiqqv(i, j),
                             splitiqgv(i, j));
        }
      }
      f10.close();
    }
  }
  io::log("");

  pdfexist = file_exists(pdffile);
  if (pdfexist) {
    io::log(" read pdfs from {}", pdffile);
    std::ifstream f10(pdffile);
    for (int i = 1; i <= 2; ++i) {
      for (int j = 1; j <= 1000; ++j) {
        const auto t = table_record(f10, 4, pdffile);
        qinqx(i, j) = table_real(t[0], pdffile);
        ginqx(i, j) = table_real(t[1], pdffile);
        qingx(i, j) = table_real(t[2], pdffile);
        gingx(i, j) = table_real(t[3], pdffile);
      }
    }
    f10.close();
  } else {
    io::log(" have to integrate pdfs, this may take some time");
    pdfint(eovest);
    pdfexist = file_exists(pdffile);
    if (!pdfexist) {
      io::log(" write pdfs to {}", pdffile);
      std::ofstream f10(pdffile);
      for (int i = 1; i <= 2; ++i) {
        for (int j = 1; j <= 1000; ++j) {
          f10 << jewel::jfmt::format(" {:.17e} {:.17e} {:.17e} {:.17e}\n", qinqx(i, j), ginqx(i, j),
                             qingx(i, j), gingx(i, j));
        }
      }
      f10.close();
    }
  }
  io::log("");

  xsecexist = file_exists(xsecfile);
  if (xsecexist) {
    io::log(" read cross sections from {}", xsecfile);
    std::ifstream f10(xsecfile);
    for (int j = 1; j <= 1001; ++j) {
      for (int jj = 1; jj <= 101; ++jj) {
        const auto t = table_record(f10, 4, xsecfile);
        intq1(j, jj) = table_real(t[0], xsecfile);
        intq2(j, jj) = table_real(t[1], xsecfile);
        intg1(j, jj) = table_real(t[2], xsecfile);
        intg2(j, jj) = table_real(t[3], xsecfile);
      }
    }
    f10.close();
  } else {
    io::log(" have to integrate cross sections, this may take some time");
    xsecint(eovest);
    xsecexist = file_exists(xsecfile);
    if (!xsecexist) {
      io::log(" write cross sections to {}", xsecfile);
      std::ofstream f10(xsecfile);
      for (int j = 1; j <= 1001; ++j) {
        for (int jj = 1; jj <= 101; ++jj) {
          f10 << jewel::jfmt::format(" {:.17e} {:.17e} {:.17e} {:.17e}\n", intq1(j, jj), intq2(j, jj),
                             intg1(j, jj), intg2(j, jj));
        }
      }
      f10.close();
    }
  }
  io::log("");
  io::log_flush();  // CALL FLUSH(3)

  // C--initialise random number generator status
  if (njob > 0) {
    MRPY(1) = njob * 1000;
    MRPY(2) = 0;
  }

  // C--Call PYR once for initialization
  r = py::pyr(0);

  // C--read random number generator from file if desired
  // library mode (JETSCAPE): the caller seeds PYR, no random-state files
  if (!cb::lib::active) {
    if (readran) {
      if (!io::fortran_open(2, "in.ran", true, true)) {
        // NOTE(port): the Fortran OPEN (status='old') has no err= branch; gfortran aborts.
        jewel::jfmt::println(stderr, "Fortran runtime error: cannot open file 'in.ran'");
        std::exit(2);
      }
      py::pyrset(2, 0);
      io::fortran_close(2);
      io::log(" read random number generator status");
    }

    // C--write random number generator state to file
    if (!io::fortran_open(2, "out.ran", true, false)) {
      jewel::jfmt::println(stderr, "Fortran runtime error: cannot open file 'out.ran'");
      std::exit(2);
    }
    py::pyrget(2, 0);
  }

  ndisc = 0;
  ngood = 0;
  nstrange = 0;
  nevdone = 0;

  errcount = 0;
  errl = 0;

  ntotspliti = 0;
  noverspliti = 0;
  ntotpdf = 0;
  noverpdf = 0;
  ntotxsec = 0;
  noverxsec = 0;
  ntotsuda = 0;
  noversuda = 0;
}

// ***********************************************************************
// ***	  subroutine initpythia
// ***********************************************************************
void initpythia(const std::string& beam1, const std::string& beam2) {
  using namespace py::pyjets;
  using namespace py::pydat1;
  using namespace py::pydat3;
  using namespace py::pysubs;
  using namespace py::pypars;
  using namespace py::pydatr;
  using namespace py::pdfas;
  using namespace py::wexpo;
  using namespace cb::pdf;
  using namespace cb::param;
  using namespace cb::disc;
  using namespace cb::weight;
  using namespace cb::errline;
  using namespace cb::evrecord;
  using namespace cb::pythiaparams;

  // C--initialise PYTHIA
  // C--keep parton shower history in PYJETS
  if (collider != "EEJJ") MSTP(125) = 2;
  // C--multiple interactions
  if (mpi) {
    MSTP(81) = 1;
  } else {
    MSTP(81) = 0;
  }
  // C--initial state radiation
  MSTP(61) = 1;
  // C--only QCD emissions in parton shower
  MSTJ(41) = 1;
  // C--switch off final state radiation off partons emitted from space-like shower
  // !	 MSTP(63)=0
  // C--switch off final state radiation
  MSTP(71) = 0;
  // C--No hadronisation (yet)
  MSTP(111) = 0;
  // C--Lambda_QCD and Q0
  PARJ(81) = lps;
  PARJ(82) = q0;
  // C--parameter affecting treatment of string corners
  PARU(14) = 1.0;
  // C--Min shat in simulation
  CKIN(1) = 2.0;
  // C--pT-cut
  CKIN(3) = ptmin;
  CKIN(4) = ptmax;
  // !C--ISR/FSR/UE tune
  // !	 MSTP(5)=0

  // C...# Tuned flavour parameters:
  PARJ(1) = r4(0.073f);
  PARJ(2) = r4(0.2f);
  PARJ(3) = r4(0.94f);
  PARJ(4) = r4(0.032f);
  PARJ(11) = r4(0.31f);
  PARJ(12) = r4(0.4f);
  PARJ(13) = r4(0.54f);
  PARJ(25) = r4(0.63f);
  PARJ(26) = r4(0.12f);
  // C...# Fragmentation tune
  MSTJ(11) = 5;
  PARJ(21) = r4(0.325f);
  PARJ(41) = 0.5;
  PARJ(42) = r4(0.6f);
  PARJ(47) = r4(0.67f);
  // !          PARJ(81) = 0.29
  // !          PARJ(82) = 1.65

  // C--particles stable if lifetime > PARJ(71), default = 10mm
  MSTJ(22) = 2;
  // C--switch off pi0 decay
  if (pi0dec) {
    MDCY(py::pycomp(111), 1) = 1;
  } else {
    MDCY(py::pycomp(111), 1) = 0;
  }

  // C--heavy flavour: keep weakly decaying charm and bottom hadrons stable (e.g. for D0, B analyses)
  if (cb::hf::stable) {
    for (int kf : {411, 421, 431, 4122, 4132, 4232, 4332, 511, 521, 531, 541, 5122, 5132, 5232, 5332})
      MDCY(py::pycomp(kf), 1) = 0;
  }

  // C--use LHAPDF
  MSTP(52) = 2;
  // C--choose pdf: CTEQ6ll (LO fit/LO alphas) - 10042
  // C	         MSTW2008 (LO central) - 21000
  MSTP(51) = pdfset;
  if (collider == "PPYQ") {
    MSEL = 0;
    MSUB(29) = 1;
  } else if (collider == "PPYG") {
    MSEL = 0;
    MSUB(14) = 1;
    MSUB(115) = 1;
  } else if (collider == "PPYJ") {
    MSEL = 0;
    MSUB(14) = 1;
    MSUB(29) = 1;
    MSUB(115) = 1;
  } else if ((collider == "PPZJ") || (collider == "PPZQ") || (collider == "PPZG") ||
             (collider == "PPDY")) {
    MSEL = 0;
    if ((collider == "PPZJ") || (collider == "PPZQ")) MSUB(30) = 1;
    if ((collider == "PPZJ") || (collider == "PPZG")) MSUB(15) = 1;
    if (collider == "PPDY") MSUB(1) = 1;
    MDME(174, 1) = 0;  // !Z decay into d dbar',
    MDME(175, 1) = 0;  // !Z decay into u ubar',
    MDME(176, 1) = 0;  // !Z decay into s sbar',
    MDME(177, 1) = 0;  // !Z decay into c cbar',
    MDME(178, 1) = 0;  // !Z decay into b bbar',
    MDME(179, 1) = 0;  // !Z decay into t tbar',
    MDME(182, 1) = 0;  // !Z decay into e- e+',
    MDME(183, 1) = 0;  // !Z decay into nu_e nu_ebar',
    MDME(184, 1) = 0;  // !Z decay into mu- mu+',
    MDME(185, 1) = 0;  // !Z decay into nu_mu nu_mubar',
    MDME(186, 1) = 0;  // !Z decay into tau- tau+',
    MDME(187, 1) = 0;  // !Z decay into nu_tau nu_taubar',
    if (channel == "ELEC") {
      MDME(182, 1) = 1;
    } else if (channel == "MUON") {
      MDME(184, 1) = 1;
    }
  } else if ((collider == "PPWJ") || (collider == "PPWQ") || (collider == "PPWG")) {
    MSEL = 0;
    if ((collider == "PPWJ") || (collider == "PPWQ")) MSUB(31) = 1;
    if ((collider == "PPWJ") || (collider == "PPWG")) MSUB(16) = 1;
    MDME(190, 1) = 0;  // ! W+ decay into dbar u,
    MDME(191, 1) = 0;  // ! W+ decay into dbar c,
    MDME(192, 1) = 0;  // ! W+ decay into dbar t,
    MDME(194, 1) = 0;  // ! W+ decay into sbar u,
    MDME(195, 1) = 0;  // ! W+ decay into sbar c,
    MDME(196, 1) = 0;  // ! W+ decay into sbar t,
    MDME(198, 1) = 0;  // ! W+ decay into bbar u,
    MDME(199, 1) = 0;  // ! W+ decay into bbar c,
    MDME(200, 1) = 0;  // ! W+ decay into bbar t,
    MDME(202, 1) = 0;  // ! W+ decay into b'bar u,
    MDME(203, 1) = 0;  // ! W+ decay into b'bar c,
    MDME(204, 1) = 0;  // ! W+ decay into b'bar t,
    MDME(206, 1) = 0;  // ! W+ decay into e+ nu_e,
    MDME(207, 1) = 0;  // ! W+ decay into mu+ nu_mu,
    MDME(208, 1) = 0;  // ! W+ decay into tau+ nu_tau,
    MDME(209, 1) = 0;  // ! W+ decay into tau'+ nu'_tau,
    if (channel == "ELEC") {
      MDME(206, 1) = 1;
    } else if (channel == "MUON") {
      MDME(207, 1) = 1;
    }
  } else if (cb::hf::hardproc == "QQ") {
    // q qbar -> q' qbar', g g -> q qbar with light flavours: switch off the c, b, t (and 4th generation)
    // channels of the gluon, which PYTHIA uses to choose the produced flavour
    MSEL = 0;
    MSUB(12) = 1;
    MSUB(53) = 1;
    for (int i = 1; i <= MDCY(21, 3); ++i) {
      const int idc = MDCY(21, 2) + i - 1;
      if (std::abs(KFDP(idc, 1)) >= 4) MDME(idc, 1) = 0;
    }
  } else if ((cb::hf::hardproc == "CC") || (cb::hf::hardproc == "BB")) {
    MSEL = 0;
    MSUB(81) = 1;
    MSUB(82) = 1;
    MSTP(7) = (cb::hf::hardproc == "CC") ? 4 : 5;
  } else {
    // C--All QCD processes are active
    MSEL = 1;
  }
  // !	 MSEL=0
  // !	 MSUB(11)=1
  // !	 MSUB(12)=1
  // !	 MSUB(53)=1
  // !	 MSUB(13)=1
  // !	 MSUB(68)=1
  // !	 MSUB(28)=1

  // C--weighted events
  if (weighted) MSTP(142) = 1;

  // C--number of errors to be printed
  MSTU(22) = std::max(10, fint(5.0f * static_cast<float>(nsim) / 100.0f));

  // C--number of lines in event record
  MSTU(4) = 23000;
  // !	MSTU(5)=23000
  MSTU(5) = 10000;

  // C--initialisation call
  if (collider == "EEJJ") {
    offset = 9;
    py::pyinit("CMS", beam1.c_str(), beam2.c_str(), sqrts);
  } else if ((collider == "PPJJ") || (collider == "PPYJ") || (collider == "PPYG") ||
             (collider == "PPYQ")) {
    offset = 8;
    py::pyinit("CMS", beam1.c_str(), beam2.c_str(), sqrts);
  } else if ((collider == "PPWJ") || (collider == "PPZJ") || (collider == "PPWQ") ||
             (collider == "PPZQ") || (collider == "PPWG") || (collider == "PPZG")) {
    offset = 10;
    py::pyinit("CMS", beam1.c_str(), beam2.c_str(), sqrts);
  } else if (collider == "PPDY") {
    // NOTE(port): OFFSET is not set for PPDY in the Fortran (keeps its value, 0 by default).
    py::pyinit("CMS", beam1.c_str(), beam2.c_str(), sqrts);
  }
}

// ***********************************************************************
// ***	  subroutine genevent
// ***********************************************************************
void genevent(const int& j, const std::string& b1, const std::string& b2) {
  using namespace cb::hepmcid;
  using namespace py::pyjets;
  using namespace py::pydat1;
  using namespace py::pydat3;
  using namespace py::pysubs;
  using namespace py::pypars;
  using namespace py::pydatr;
  using namespace cb::param;
  using namespace cb::disc;
  using namespace cb::angor;
  using namespace cb::evoltimes;
  using namespace cb::ftimefac;
  using namespace cb::time;
  using namespace cb::colour;
  using namespace cb::weight;
  using namespace py::wexpo;
  using namespace cb::rapmax;
  using namespace cb::jetpoint;
  using namespace cb::evrecord;
  using namespace cb::totnumevents;
  using namespace cb::storescatcen;
  using namespace cb::speclines;

  // C--Variables local to this program
  // NOTE(port): PID is only assigned after the PPDY branch; there it is passed to
  // CONVERTTOHEPMC uninitialised (harmless: only used for EEJJ). Initialised to 0.
  int pid = 0;
  int tnowsec, etasec, etah, etam, etas;
  int tnow[3];  // integer*4 tnow(3), 0-based here: tnow[0]=h, [1]=m, [2]=s
  double z1, z2;
  bool printmessage;
  float time, rtime;  // real time, rtime
  char tunit, rtunit;

  N = 0;
  colmax = 600;
  discard = false;
  for (int i = 1; i <= 23000; ++i) {
    mv(i, 1) = 0.0;
    mv(i, 2) = 0.0;
    mv(i, 3) = 0.0;
    mv(i, 4) = 0.0;
    mv(i, 5) = 0.0;
    zd(i) = 0.0;
    thetad(i) = 0.0;
    qqbard(i) = false;
    cb::hf::qqbarfl(i) = 0;
    starttime(i) = -1.e6;
    medind(i) = false;
  }
  nscatcen = 0;

  medium().mednextevt();

  // C--initialisation with matrix element
  // C--production vertex
  medium().pickvtx(x0, y0);
  ltime = medium().getltimemax();

  // l99:
  py::pyevnt();
  npart = N - offset;
  evweight = PARI(10);
  sumofweights = sumofweights + evweight;
  if ((collider == "EEJJ") && (std::abs(K(8, 2)) > 6)) {
    wdisc = wdisc + evweight;
    ndisc = ndisc + 1;
    goto l102;
  } else {
    ngood = ngood + 1;
  }

  // C--DY: don't have to do anything
  if (collider == "PPDY") {
    py::pyexec();
    converttohepmc(hpmcfid, ngood, pid, b1, b2);
    goto l102;
  }

  // C--prepare event record
  interpretpythiaevent();
  hfcheck("after interpretpythiaevent");

  findspeclines();

  setcolour();

  // C--set status codes
  if (K(lme1, 1) < 11) K(lme1, 1) = 1;
  if (K(lme2, 1) < 11) K(lme2, 1) = 1;
  pid = K(lme1, 2);
  for (int ipart = offset + 1, ipart_end = offset + npart; ipart <= ipart_end; ++ipart) {
    if ((!isqorg(K(ipart, 2))) && (K(ipart, 1) < 11)) {
      K(ipart, 1) = 4;
    } else if ((ipart != lme1) && (ipart != lme2) && (K(ipart, 1) < 11)) {
      if (P(ipart, 5) > py::pymass(K(ipart, 2)) * r4(1.00001f)) {
        K(ipart, 1) = 1;
      } else {
        if ((isqorg(K(ipart, 2))) && (std::abs(py::pyp(ipart, 19)) < etamax)) {
          // heavy flavour: on-shell heavy quarks keep their mass
          if ((P(ipart, 5) > 0.0) && (!isheavy(K(ipart, 2)))) makemassless(ipart, lme1, lme2);
          K(ipart, 1) = 1;
        } else {
          K(ipart, 1) = 4;
        }
      }
    }
    // C--let (semi-)hard partons from MPI's interact in background medium
    // !	   if ((k(ipart,3).eq.0).and.(abs(k(ipart,2)).lt.22)) k(ipart,1)=1
    // C--end MPI treatment
  }

  {
    const int zero = 0;
    const bool false_ = false;
    findvirtualities(zero, zero, x0, y0, z1, z2, false_);
    hfcheck("after findvirtualities");
  }

  findmpivirtualities(x0, y0);

  // C--develop parton shower
  makecascade();
  if (discard) {
    ngood = ngood - 1;
    wdisc = wdisc + evweight;
    ndisc = ndisc + 1;
    io::log(" discard event {}", j);
    goto l102;
  }
  // heavy flavour: Langevin transport of final-state heavy quarks
  if (cb::hf::on && (cb::hf::hqds > 0.0)) hfdiffuse();

  if (!allhad) {
    for (int i = 1, i_end = N; i <= i_end; ++i) {
      if (K(i, 1) == 3) K(i, 1) = 22;
    }
  }
  if (hadro) {
    makestrings(hadrotype);
    if (discard) {
      io::log(" discard event {}", j);
      wdisc = wdisc + evweight;
      ndisc = ndisc + 1;
      ngood = ngood - 1;
      goto l102;
    }
    py::pyexec();
    if (MSTU(30) != errcount) {
      io::log(" PYTHIA discards event {}  (error number {})", j, MSTU(30));
      errcount = MSTU(30);
      wdisc = wdisc + evweight;
      ndisc = ndisc + 1;
      ngood = ngood - 1;
      goto l102;
    }
  }

  if (MSTU(30) != errcount) {
    errcount = MSTU(30);
  } else {
    if (dosubtraction) {
      // NOTE(port): the Fortran passes the CHARACTER*4 collider without a
      // terminating NUL; subtraction.C uses strcmp on it, which reads on into
      // the next COMMON member (LOGICAL hadro, little endian: 1/0 then 0,0,0).
      // The same memory image is reproduced here. (strcmp's result is also used
      // with inverted logic in subtraction.C -- kept.)
      char proc[8] = {' ', ' ', ' ', ' ', static_cast<char>(hadro ? 1 : 0), 0, 0, 0};
      for (std::size_t c = 0; c < 4 && c < collider.size(); ++c) proc[c] = collider[c];
      subtract_thmom_(reinterpret_cast<int(*)[23000]>(K.data()),
                      reinterpret_cast<double(*)[23000]>(P.data()), &N,
                      reinterpret_cast<double(*)[23000]>(scatcen.data()), &nscatcen, &rsub, proc);
      if (N == 0) {
        io::log(" event too long after subtraction ");
        io::log(" discard event {}", j);
        wdisc = wdisc + evweight;
        ndisc = ndisc + 1;
        ngood = ngood - 1;
        goto l102;
      }
      putparticlesonshell();
    }
    converttohepmc(hpmcfid, ngood, pid, b1, b2);
  }

  // C--write message to log-file
l102:
  printmessage = false;
  if (nsim <= 100) {
    printmessage = true;
  } else if (j % (nsim / 100) == 0) {
    printmessage = true;
  }

  nevdone = nevdone + 1;

  if (printmessage) {
    io::log(" done with event number {}", j);
    // !	  call printtime
    {  // call itime(tnow): local wall-clock time
      const std::time_t now = std::time(nullptr);
      std::tm lt{};
      localtime_r(&now, &lt);
      tnow[0] = lt.tm_hour;
      tnow[1] = lt.tm_min;
      tnow[2] = lt.tm_sec;
    }
    tnowsec = tnow[0] * 3600 + tnow[1] * 60 + tnow[2];
    // call etime(tarray,time): CPU time (user+system) of the process
    time = static_cast<float>(static_cast<double>(std::clock()) / CLOCKS_PER_SEC);
    rtime = static_cast<float>(nsim - nevdone) * time / static_cast<float>(nevdone);
    etasec = static_cast<int>(static_cast<float>(tnowsec) + rtime);
    etah = etasec / 3600;
    etam = (etasec - etah * 3600) / 60;
    etas = (etasec - etah * 3600 - etam * 60);
    if (etah >= 24) etah = etah - 24;
    if (time < 60.0f) {
      tunit = 's';
    } else if (time < 3600.0f) {
      tunit = 'm';
      time = time / 60.0f;
    } else {
      tunit = 'h';
      time = time / 3600.0f;
    }
    if (rtime < 60.0f) {
      rtunit = 's';
    } else if (rtime < 3600.0f) {
      rtunit = 'm';
      rtime = rtime / 60.0f;
    } else {
      rtunit = 'h';
      rtime = rtime / 3600.0f;
    }
    // 1000 format ('   ',i4,A,' ',f6.2,A3,A,f6.2,A3,A)
    io::log_line(jewel::jfmt::format("   {:4d}{} {:6.2f}{:>3}{}{:6.2f}{:>3}{}", nevdone * 100 / nsim,
                             "% of the run completed in ", time, tunit, " -> ", rtime, rtunit,
                             " remaining"));
    // 1001 format ('   ',A,' ',i2.2,':',i2.2,':',i2.2,A,' ',i2.2,':',i2.2,':',i2.2)
    io::log_line(jewel::jfmt::format("   {} {:02d}:{:02d}:{:02d}{} {:02d}:{:02d}:{:02d}", "time now is ",
                             tnow[0], tnow[1], tnow[2], "  -->  ETA: ", etah, etam, etas));
    py::pyrget(2, -1);
    io::fortran_flush(2);
  }
  io::log_flush();
}

}  // namespace jewel
