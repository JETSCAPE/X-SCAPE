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

#ifndef HADRONICEMT_H
#define HADRONICEMT_H

#include "sigslot.h"
#include "BulkMediaInfo.h"
#include "JetScapeParticles.h"
#include <memory>

namespace Jetscape {

class HadronicEMT {
private:
  double dx_, dy_, dz_; // in [fm]

public:
  HadronicEMT();
  HadronicEMT(double dx, double dy, double dz);
  ~HadronicEMT() {};

  void InitTask();

  /// Fill in bulk media info for (t,x,y,z) from current hadron list
  void GetBulkInfo(Jetscape::real t, Jetscape::real x, Jetscape::real y, 
            Jetscape::real z, std::unique_ptr<BulkMediaInfo> &bulk_info_ptr,
            std::vector<Hadron> &h_list);

  double get_t(double tau, double eta) const;
  double get_z(double tau, double eta) const;

  double get_ptau(double px, double pz, double eta) const;
  double get_peta(double px, double pz, double eta) const;


};

}; // namespace Jetscape

#endif // HADRONICEMT_H