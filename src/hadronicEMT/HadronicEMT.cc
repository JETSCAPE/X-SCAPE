/*******************************************************************************
 * Copyright (c) The JETSCAPE Collaboration, 2018
 *
 * Modular, task-based framework for simulating all aspects of heavy-ion collisions
 * 
 * For the list of contributors see AUTHORS.
 *
 * Report issues at https://github.com/JETSCAPE/JETSCAPE/issues
 *
 * or via email to bugs.jetscape@gmail.com
 *
 * Distributed under the GNU General Public License 3.0 (GPLv3 or later).
 * See COPYING for details.
 ******************************************************************************/

#include "HadronicEMT.h"
#include "JetScapeSignalManager.h"
#include "BulkMediaInfo.h"

namespace Jetscape {

HadronicEMT::HadronicEMT(){
  VERBOSE(8);
  dx_ = 0.3;
  dy_ = 0.3;
  dz_ = 0.3;
  InitTask(); // Get values of parameters from XML 
}

HadronicEMT::HadronicEMT(double dx, double dy, double dz){
  VERBOSE(8);
  JSINFO << "Initialize HadronicEMT ...";
  dx_ = dx;
  dy_ = dy;
  dz_ = dz;
  JSINFO << "<HadronicEMT> cell size: dx = " << dx_ << " fm, dy = " 
    << dy_ << " fm, dz = " << dz_ << " fm";
}

void HadronicEMT::InitTask(){
    JSINFO << "Initialize HadronicEMT ...";

    dx_ = JetScapeXML::Instance()->GetElementDouble({"HadronicEMT", "dx"});
    dy_ = JetScapeXML::Instance()->GetElementDouble({"HadronicEMT", "dy"});
    dz_ = JetScapeXML::Instance()->GetElementDouble({"HadronicEMT", "dz"});
}

void HadronicEMT::GetBulkInfo(Jetscape::real t, Jetscape::real x,
                              Jetscape::real y, Jetscape::real z,
                              std::unique_ptr<BulkMediaInfo> &bulk_info_ptr,
                              std::vector<Hadron> &h_list) {
  bulk_info_ptr = make_unique<BulkMediaInfo>();

  // compute the T^munu in units of GeV/fm^3 from the hadron list in the lab frame
  Jetscape::real T00, T01, T02, T03, T11, T12, T13, T22, T23, T33;
  T00 = T01 = T02 = T03 = T11 = T12 = T13 = T22 = T23 = T33 = 0.;
  const double dV = dx_*dy_*dz_;
  for (const auto ihad: h_list) {
    const FourVector r = ihad.x_in();
    // check if the hadron is in the surrounding of the requested point
    // time has to match exactly (do we want to keep this check??, if time stepped evolution there is no need for it...)
    if ((std::abs(t-r.t()) < rounding_error) && (std::abs(x-r.x()) <= dx_/2.)
        && (std::abs(y-r.y()) <= dy_/2.) && (std::abs(z-r.z()) <= dz_/2.)) {
      T00 += ihad.e();
      T01 += ihad.px();
      T02 += ihad.py();
      T03 += ihad.pz();
      T11 += ihad.px() * ihad.px() / ihad.e();
      T12 += ihad.px() * ihad.py() / ihad.e();
      T13 += ihad.px() * ihad.pz() / ihad.e();
      T22 += ihad.py() * ihad.py() / ihad.e();
      T23 += ihad.py() * ihad.pz() / ihad.e();
      T33 += ihad.pz() * ihad.pz() / ihad.e();
    }
    // fill T^{\mu\nu}
    bulk_info_ptr->tmn[0][0] = T00/dV;
    bulk_info_ptr->tmn[0][1] = T01/dV;
    bulk_info_ptr->tmn[0][2] = T02/dV;
    bulk_info_ptr->tmn[0][3] = T03/dV;
    bulk_info_ptr->tmn[1][1] = T11/dV;
    bulk_info_ptr->tmn[1][2] = T12/dV;
    bulk_info_ptr->tmn[1][3] = T13/dV;
    bulk_info_ptr->tmn[2][2] = T22/dV;
    bulk_info_ptr->tmn[2][3] = T23/dV;
    bulk_info_ptr->tmn[3][3] = T33/dV;

    // fill the rest of T^{\mu\nu} by symmetry
    bulk_info_ptr->tmn[1][0] = T01/dV;
    bulk_info_ptr->tmn[2][0] = T02/dV;
    bulk_info_ptr->tmn[2][1] = T12/dV;
    bulk_info_ptr->tmn[3][0] = T03/dV;
    bulk_info_ptr->tmn[3][1] = T13/dV;
    bulk_info_ptr->tmn[3][2] = T23/dV;

    // provide the energy density
    bulk_info_ptr->energy_density = T00/dV;
  }
}

// Get Cartesian time t from tau and eta
double HadronicEMT::get_t(double tau, double eta) const {
  return tau*cosh(eta);
}

// Get Cartesian coordinate z from tau and eta
double HadronicEMT::get_z(double tau, double eta) const {
  return tau*sinh(eta);
}

// Lorentz Transformation to get tau component of four vector
double HadronicEMT::get_ptau(double p0, double p3, double eta) const {
  return p0*cosh(eta) - p3*sinh(eta);
}

// Lorentz Transformation to get eta component of four vector
double HadronicEMT::get_peta(double p0, double p3, double eta) const {
  return p3*cosh(eta) - p0*sinh(eta);
}

};