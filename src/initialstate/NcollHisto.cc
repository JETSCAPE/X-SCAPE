/*******************************************************************************
 * Copyright (c) The JETSCAPE Collaboration, 2018
 *
 * Modular, task-based framework for simulating all aspects of heavy-ion
 *collisions
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

#include "NcollHisto.h"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <sstream>

#include "JetScapeLogger.h"

RegisterJetScapeModule<NcollHisto> NcollHisto::reg("CustomModuleNcollHisto");

NcollHisto::NcollHisto() { SetId("NcollHisto"); }

NcollHisto::~NcollHisto() {}

void NcollHisto::InitTask() {
  auto *e = GetXMLElement({"IS", "CustomModuleNcollHisto", "file"});
  std::string file = (e && e->GetText()) ? e->GetText() : "ncoll_histo.dat";
  file.erase(0, file.find_first_not_of(" \t\n"));
  file.erase(file.find_last_not_of(" \t\n") + 1);
  std::ifstream in(file);
  std::string line;
  double sum = 0.0;
  while (std::getline(in, line)) {
    for (char &c : line)
      if (c == ',') c = ' ';
    std::istringstream is(line);
    double x, y, n;
    if (!(is >> x >> y >> n)) continue;  // header
    if (x_.size() == 1 && dx_ == 0.0 && y != y_.back()) dx_ = std::fabs(y - y_.back());
    if (n <= 0.0) {
      if (x_.empty() || dx_ == 0.0) {  // keep the first points to find the bin size
        x_.push_back(x);
        y_.push_back(y);
        cumulative_.push_back(sum);
      }
      continue;
    }
    sum += n;
    x_.push_back(x);
    y_.push_back(y);
    cumulative_.push_back(sum);
  }
  if (sum <= 0.0 || dx_ <= 0.0) {
    JSWARN << "NcollHisto: no binary collisions in " << file;
    exit(-1);
  }
  ncoll_ = sum;
  JSINFO << "NcollHisto: " << sum << " entries in " << x_.size() << " bins of "
         << dx_ << " fm from " << file;
}

void NcollHisto::SampleABinaryCollisionPoint(double &t, double &x, double &y,
                                             double &z) {
  auto &rng = *GetMt19937Generator();
  const double r = uniform_(rng) * ncoll_;
  const size_t k = std::upper_bound(cumulative_.begin(), cumulative_.end(), r) -
                   cumulative_.begin();
  const size_t i = std::min(k, x_.size() - 1);
  x = x_[i] + (uniform_(rng) - 0.5) * dx_;
  y = y_[i] + (uniform_(rng) - 0.5) * dx_;
  t = 0.0;
  z = 0.0;
}
