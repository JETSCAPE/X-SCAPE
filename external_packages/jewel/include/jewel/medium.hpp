// JEWEL 2.6.0 C++ port -- medium model interface.
//
// In the Fortran code the medium is a set of free routines provided by
// either medium-simple.f or medium-vac.f, selected at link time. Here it is
// an abstract class; the core calls the active model through medium().
// Argument conventions follow routines.hpp: `const T&` for inputs, `T&` for
// outputs (Fortran pass-by-reference).
#pragma once

#include <memory>
#include <string>

namespace jewel {

class Medium {
 public:
  virtual ~Medium() = default;

  // MEDINIT(FILE,id,etam,mass)
  virtual void medinit(const std::string& file, const int& id, const double& etam, const int& mass) = 0;
  virtual void mednextevt() = 0;
  virtual double getcentrality() = 0;
  virtual void pickvtx(double& x, double& y) = 0;
  virtual void setb(const double& bval) = 0;
  // GETSCATTERER(X,Y,Z,T,TYPE,PX,PY,PZ,E,MS)
  virtual void getscatterer(const double& x, const double& y, const double& z, const double& t,
                            int& type, double& px, double& py, double& pz, double& e,
                            double& ms) = 0;
  virtual void avscatcen(const double& x, const double& y, const double& z, const double& t,
                         double& px, double& py, double& pz, double& e, double& m) = 0;
  virtual void maxscatcen(double& px, double& py, double& pz, double& e, double& m) = 0;
  virtual double getmd(const double& x, const double& y, const double& z, const double& t) = 0;
  virtual double getms(const double& x, const double& y, const double& z, const double& t) = 0;
  virtual double getneff(const double& x, const double& y, const double& z, const double& t) = 0;
  // Density of scattering centres seen by a probe with lab-frame velocity (bx,by,bz), i.e. including
  // the flux factor.  The models of the JEWEL distribution put the Bjorken flux factor into GETNEFF
  // itself and do not need the velocity; media with general flow (JETSCAPE) override this.
  virtual double getneffprobe(const double& x, const double& y, const double& z, const double& t,
                              const double& bx, const double& by, const double& bz) {
    static_cast<void>(bx);
    static_cast<void>(by);
    static_cast<void>(bz);
    return getneff(x, y, z, t);
  }
  virtual double gettemp(const double& x, const double& y, const double& z, const double& t) = 0;
  virtual double gettempmax() = 0;
  virtual double getmdmax() = 0;
  virtual double getmdmin() = 0;
  virtual double getmsmax() = 0;
  virtual double getnatmdmin() = 0;
  virtual double getltimemax() = 0;
  virtual double getneffmax() = 0;
  // MEDDERIV(XVAL,W): integrand used by the core's ODEINT (via DERIV)
  virtual double medderiv(const double& xval, const int& w) = 0;
};

// The active medium model (set once in main()).
Medium& medium();
void set_medium(std::unique_ptr<Medium> m);

std::unique_ptr<Medium> make_medium_simple();  // medium-simple.f
std::unique_ptr<Medium> make_medium_vac();     // medium-vac.f

}  // namespace jewel
