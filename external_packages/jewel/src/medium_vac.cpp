// JEWEL 2.6.0 C++ port -- vacuum "medium" (no medium). Translated from
// medium-vac.f:61-214.
//
// The COMMON blocks MEDPARAM and MEDIUM of medium-vac.f are private data
// members of MediumVac.
#include <memory>
#include "jewel/jfmt.hpp"
#include <string>

#include "jewel/medium.hpp"

namespace jewel {

namespace {

class MediumVac final : public Medium {
 public:
  void medinit(const std::string& file, const int& id, const double& etam, const int& mass) override;
  void mednextevt() override;
  double getcentrality() override;
  void pickvtx(double& x, double& y) override;
  void setb(const double& bval) override;
  void getscatterer(const double& x, const double& y, const double& z, const double& t, int& type,
                    double& px, double& py, double& pz, double& e, double& ms) override;
  void avscatcen(const double& x, const double& y, const double& z, const double& t, double& px,
                 double& py, double& pz, double& e, double& ms) override;
  void maxscatcen(double& px, double& py, double& pz, double& e, double& ms) override;
  double getneff(const double& x, const double& y, const double& z, const double& t) override;
  double gettemp(const double& x, const double& y, const double& z, const double& t) override;
  double gettempmax() override;
  double getmdmax() override;
  double getmdmin() override;
  double getmsmax() override;
  double getneffmax() override;
  double getnatmdmin() override;
  double getltimemax() override;
  double getmd(const double& x, const double& y, const double& z, const double& t) override;
  double getms(const double& x, const double& y, const double& z, const double& t) override;
  double medderiv(const double& xval, const int& w) override;

 private:
  // COMMON/MEDPARAM/BMIN,BMAX,CENTRMIN,CENTRMAX,BREAL,RAU
  // (only RAU is ever assigned; the other members are never referenced)
  [[maybe_unused]] double bmin = 0, bmax = 0, centrmin = 0, centrmax = 0, breal = 0;
  double rau = 0;
  // COMMON/MEDIUM/MEDIUM
  [[maybe_unused]] bool medium = false;  // DATA MEDIUM/.FALSE./
};

// NOTE(port): the Fortran MEDINIT of medium-vac.f takes no arguments although
// the core calls MEDINIT(FILEMED,logfid,etamax,intmass); the arguments are
// ignored here.
void MediumVac::medinit([[maybe_unused]] const std::string& file, [[maybe_unused]] const int& id,
                        [[maybe_unused]] const double& etam, [[maybe_unused]] const int& mass) {
  rau = 5.0;
}

void MediumVac::mednextevt() {}

double MediumVac::getcentrality() { return -1.0; }

void MediumVac::pickvtx(double& x, double& y) {
  x = 0.0;
  y = 0.0;
}

void MediumVac::setb([[maybe_unused]] const double& bval) {}

// NOTE(port): in medium-vac.f GETSCATTERER has two extra dummy arguments
// (X,Y,Z,T,TYPE,PX,PY,PZ,E,MS,MD,TEMP) that the core does not pass; none of
// the arguments is referenced. The interface signature is used and, as in the
// Fortran, no output is assigned (the caller's values stay unchanged).
void MediumVac::getscatterer([[maybe_unused]] const double& x, [[maybe_unused]] const double& y,
                             [[maybe_unused]] const double& z, [[maybe_unused]] const double& t,
                             [[maybe_unused]] int& type, [[maybe_unused]] double& px,
                             [[maybe_unused]] double& py, [[maybe_unused]] double& pz,
                             [[maybe_unused]] double& e, [[maybe_unused]] double& ms) {
  jewel::jfmt::println(" GETSCATTERER called although in vacuum");
}

// As in the Fortran, no output is assigned.
void MediumVac::avscatcen([[maybe_unused]] const double& x, [[maybe_unused]] const double& y,
                          [[maybe_unused]] const double& z, [[maybe_unused]] const double& t,
                          [[maybe_unused]] double& px, [[maybe_unused]] double& py,
                          [[maybe_unused]] double& pz, [[maybe_unused]] double& e,
                          [[maybe_unused]] double& ms) {
  jewel::jfmt::println(" AVSCATCEN called although in vacuum");
}

void MediumVac::maxscatcen(double& px, double& py, double& pz, double& e, double& ms) {
  ms = 0.5;
  e = 0.5;
  px = 0.0;
  py = 0.0;
  pz = 0.0;
}

double MediumVac::getneff([[maybe_unused]] const double& x, [[maybe_unused]] const double& y,
                          [[maybe_unused]] const double& z, [[maybe_unused]] const double& t) {
  return 0.0;
}

double MediumVac::gettemp([[maybe_unused]] const double& x, [[maybe_unused]] const double& y,
                          [[maybe_unused]] const double& z, [[maybe_unused]] const double& t) {
  return 0.0;
}

double MediumVac::gettempmax() { return 0.0; }

double MediumVac::getmdmax() { return 0.0; }

double MediumVac::getmdmin() { return 0.0; }

double MediumVac::getmsmax() { return 0.0; }

double MediumVac::getneffmax() { return 0.0; }

double MediumVac::getnatmdmin() { return 0.0; }

double MediumVac::getltimemax() { return 0.0; }

double MediumVac::getmd([[maybe_unused]] const double& x, [[maybe_unused]] const double& y,
                        [[maybe_unused]] const double& z, [[maybe_unused]] const double& t) {
  return 0.0;
}

double MediumVac::getms([[maybe_unused]] const double& x, [[maybe_unused]] const double& y,
                        [[maybe_unused]] const double& z, [[maybe_unused]] const double& t) {
  return 0.0;
}

// NOTE(port): MEDDERIV in medium-vac.f has no IMPLICIT NONE (XVAL, W are
// implicitly REAL); it only returns 0 and is never called in vacuum.
double MediumVac::medderiv([[maybe_unused]] const double& xval, [[maybe_unused]] const int& w) {
  return 0.0;
}

}  // namespace

std::unique_ptr<Medium> make_medium_vac() { return std::make_unique<MediumVac>(); }

}  // namespace jewel
