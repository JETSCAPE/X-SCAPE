// JEWEL 2.6.0 C++ port -- prototypes of all routines of the JEWEL core.
//
// Generated from jewel-2.6.0.f (tools/fsig.py) and then maintained by hand.
// One C++ function per Fortran routine, same name in lower case. Fortran passes
// every argument by reference; here arguments the routine (transitively)
// assigns are `T&`, all others `const T&`, which keeps the Fortran aliasing
// semantics. Comments give the Fortran source line of the original routine.
#pragma once

#include <string>

namespace jewel {

void init();  // jewel-2.6.0.f:281
// INIT for use as a library (JETSCAPE): reads `paramfile` like argv[1]; opens no HepMC or out.ran file
void init_library(const std::string& paramfile);
void initpythia(const std::string& beam1, const std::string& beam2);  // jewel-2.6.0.f:865
void genevent(const int& j, const std::string& b1, const std::string& b2);  // jewel-2.6.0.f:1107
void interpretpythiaevent();  // jewel-2.6.0.f:1405
void findspeclines();  // jewel-2.6.0.f:1748
void setcolour();  // jewel-2.6.0.f:1817
void findvirtualities(const int& line1, const int& line2, const double& x0, const double& y0, double& z1, double& z2, const bool& mpiscat);  // jewel-2.6.0.f:1893
void findmpivirtualities(const double& x0, const double& y0);  // jewel-2.6.0.f:2246
double getthetam(const int& line);  // jewel-2.6.0.f:2303
void makestrings(const int& which);  // jewel-2.6.0.f:2372
void makestrings_vac();  // jewel-2.6.0.f:2391
void makestrings_minl();  // jewel-2.6.0.f:2625
void cleanup(const int& nfirst);  // jewel-2.6.0.f:2877
void makecascade();  // jewel-2.6.0.f:2902
void makebranch(const int& l);  // jewel-2.6.0.f:2948
void makesplitting(const int& l);  // jewel-2.6.0.f:3288
void rejectsplitting(const int& line, const int& mother);  // jewel-2.6.0.f:3688
void makeinsplit(const int& l, double& x, const double& tsum, const double& virt, const int& typi, const double& time, const double& taurad);  // jewel-2.6.0.f:3871
void doinstatescat(int& l, double& x, int& typi, double& q, const double& tstart, double& deltat, bool& overq0, const bool& retrysplit);  // jewel-2.6.0.f:4132
void docohscattering(int& l, const double& tstart, const double& tsum, const double& splitangle, const double& splittime);  // jewel-2.6.0.f:4644
void dofistatescat(const int& l, const double& tnow, const double& dtleft, double& deltat, double& newmass, bool& overq0, double& z, bool& qqbar);  // jewel-2.6.0.f:5080
double getnewmass(const int& l, const double& q2, const double& qold2, const double& mass, const bool& in, const double& x, double& zdec, bool& qqbardec);  // jewel-2.6.0.f:5253
double getpnorad1(const int& line, const double& x, const double& y, const double& z, const double& t);  // jewel-2.6.0.f:5366
void getqvec(const int& l, const int& j, const double& dt, const double& x, const bool& recalc);  // jewel-2.6.0.f:5437
void dokinematics(const int& l, const int& lold, const int& n1, int& n2, const double& newm, bool& retrysplit, const double& time, const double& x, double& z, const bool& qqbar, const bool& forcohscat);  // jewel-2.6.0.f:5643
double getproba(const double& qi, const double& qf, const double& qaa, const double& zaa, const double& ebb, const std::string& type, const double& t1, const bool& ins2);  // jewel-2.6.0.f:6459
double getsudakov(const double& qmax1, const double& qa1, const double& qb1, const double& za1, const double& eb1, const std::string& type3, const double& t2, const bool& ins);  // jewel-2.6.0.f:6486
double getinsudakov(const double& qb, const double& qmax1, const std::string& type3);  // jewel-2.6.0.f:6545
double deriv(const double& xval, const int& w4);  // jewel-2.6.0.f:6588
double getspliti(const double& qa, const double& qb, const double& zeta, const double& eb, const std::string& type1);  // jewel-2.6.0.f:6692
double getinspliti(const double& qb, const std::string& type1);  // jewel-2.6.0.f:6837
double getpdf(const double& x, const double& q, const std::string& typ);  // jewel-2.6.0.f:6928
double getpdfxint(const double& q, const std::string& typ);  // jewel-2.6.0.f:7027
double getpdfxintexact(const double& q, const std::string& typ);  // jewel-2.6.0.f:7096
double getxsecint(const double& tm, const double& md, const std::string& typ2);  // jewel-2.6.0.f:7135
double getinsudafast(const double& q1, const double& q2, const std::string& typ);  // jewel-2.6.0.f:7247
double getinsudared(const double& q, const std::string& typ2);  // jewel-2.6.0.f:7284
double getsscat(const double& en, const double& px, const double& py, const double& pz, const double& mp, const double& lw, const char& type1, const char& type2, const double& x, const double& y, const double& z, const double& t, const int& mode);  // CHARACTER*1 args: pass first char  // jewel-2.6.0.f:7365
double getmass(const double& qbmin, const double& qbmax, const double& theta, const double& ep, const std::string& type, const double& max2, const bool& ins, double& zdec, bool& qqbardec, const double mq = 0.0, int* qqfl = nullptr);  // jewel-2.6.0.f:7465
double generatez(const double& ti, const double& ea, const double& epsi, const std::string& type);  // jewel-2.6.0.f:7632
double gett(const double& mint, const double& maxt, const double& mdeb);  // jewel-2.6.0.f:7834
double scatprimfunc(const double& t, const double& mdeb);  // jewel-2.6.0.f:7692
double intpqq(const double& z, const double& q);  // jewel-2.6.0.f:7717
double intpgglow(const double& z, const double& q);  // jewel-2.6.0.f:7740
double intpgghigh(const double& z, const double& q);  // jewel-2.6.0.f:7762
double intpqglow(const double& z, const double& q);  // jewel-2.6.0.f:7784
double intpqghigh(const double& z, const double& q);  // jewel-2.6.0.f:7809
double ei(const double& x);  // jewel-2.6.0.f:7872
double pqq(const double& z);  // jewel-2.6.0.f:7927
double pgq(const double& z);  // jewel-2.6.0.f:7938
double pgg(const double& z);  // jewel-2.6.0.f:7949
double pqg(const double& z);  // jewel-2.6.0.f:7960
double alphas(const double& t, const double& lambda);  // jewel-2.6.0.f:7971
void splitfncint(const double& emax);  // jewel-2.6.0.f:7994
void pdfint(const double& emax);  // jewel-2.6.0.f:8082
void xsecint(const double& emax);  // jewel-2.6.0.f:8142
void insudaint(const double& emax);  // jewel-2.6.0.f:8235
void eixint();  // jewel-2.6.0.f:8282
void odeint(double& ystart, const double& a, const double& b, const double& eps, const double& h1, const double& hmin, const int& w1);  // jewel-2.6.0.f:8315
void rkstepper(double& x, double& y, const double& dydx, const double& htest, double& hdid, double& hnew, const double& yscale, const double& eps, const int& w1);  // jewel-2.6.0.f:8354
double rk4step(const double& x, const double& y, const double& dydx, const double& h, const int& w1);  // jewel-2.6.0.f:8407
bool getdeltat(const int& line, const double& tstart, const double& dtmax1, double& deltat, const bool& elas);  // jewel-2.6.0.f:8423
int poissonian(const double& lambda);  // jewel-2.6.0.f:8560
void makemassless(const int& l3, const int& l1, const int& l2);  // jewel-2.6.0.f:8588
bool ishadron(const int& id);  // jewel-2.6.0.f:8638
bool isdiquark(const int& id);  // jewel-2.6.0.f:8658
bool islepton(const int& id);  // jewel-2.6.0.f:8676
bool isparton(const int& id);  // jewel-2.6.0.f:8690
bool isqorg(const int& id);  // jewel-2.6.0.f:8705
bool isprimstring(const int& l);  // jewel-2.6.0.f:8721
bool issecstring(const int& l);  // jewel-2.6.0.f:8749
bool isprimhadron(const int& l);  // jewel-2.6.0.f:8781
void getposition(const int& line, const double& time, double& posx, double& posy, double& posz);  // jewel-2.6.0.f:8807
bool compressevent(int& l1);  // jewel-2.6.0.f:8832
// library mode: origin of a line (0 shower, 1 medium recoil, 2 deleted); see utilities_hepmc.cpp
int lib_origin(int l);
void pevrec(const int& num, const bool& col);  // jewel-2.6.0.f:8919
void putparticlesonshell();  // jewel-2.6.0.f:8959
void converttohepmc(const int& j, const int& evnum, const int& pid, const std::string& beam1, const std::string& beam2);  // jewel-2.6.0.f:9112
void copyline(const int& nfr, const int& nto, const int& mode);  // jewel-2.6.0.f:9547
void printlogo(const int& fid);  // jewel-2.6.0.f:9594
void printtime();  // jewel-2.6.0.f:9683
void eix(const double& x, double& ei);  // meix.f:7 (Zhang & Jin, exponential integral)

// C-linkage entry point of subtraction.cpp (thermal momentum subtraction)
extern "C" void subtract_thmom_(int pyjetsK[5][23000], double pyjetsP[5][23000], int* pyjetsN,
                                double scatcens[5][23000], int* nscatcen, double* Rmax, char* proc);

// Program entry (PROGRAM JEWEL); argv[1] is the optional parameter file.
int run_jewel(int argc, char** argv);

}  // namespace jewel
