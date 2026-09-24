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

#ifndef COLOREDHADRONIZATIONSINGLET_H
#define COLOREDHADRONIZATIONSINGLET_H

#include "HadronizationModule.h"
#include "Pythia8/Pythia.h"
#include "JetScapeSignalManager.h"

using namespace Jetscape;

#define VORONOI_SHIFT_TYPE 1
#define MC_SHIFT_TYPE 2

class ColoredHadronizationSinglet : public HadronizationModule<ColoredHadronizationSinglet> {
public:
  //used for:
  //1) hadrons have a position and some radius
  //2) empty space is tracked by spheres marked by center and radius 
  struct Sphere {
    double x; double y; double z;
    double r;
  };

  ColoredHadronizationSinglet();
  virtual ~ColoredHadronizationSinglet();

  void InitTask();
  void DoHadronization(vector<vector<shared_ptr<Parton>>> &shower,
                       vector<shared_ptr<Hadron>> &hOut,
                       vector<shared_ptr<Parton>> &pOut);
  void WriteTask(weak_ptr<JetScapeWriter> w);

  void GetCircumSphere(std::vector<double> A, std::vector<double> B, std::vector<double> C, std::vector<double> D, double &cx, double &cy, double &cz, double &cr);
  int IsInSphere(double cx, double cy, double cz, double cr, std::vector<double> P);
  int IsInCircumsphere(std::vector<double> A, std::vector<double> B, std::vector<double> C, std::vector<double> D, std::vector<double> P);
  Sphere InitSphere(double x, double y, double z, int pid);
  double d2Sphere(struct Sphere A, struct Sphere B);
  bool isOverlap(double *x, int pid, std::vector<ColoredHadronizationSinglet::Sphere> nucleonSpheres);
  std::vector<std::vector<int>> DelaunayTriangulate3D(std::vector<std::vector<double>> coords);
  std::vector<struct Sphere> VoronoiHoles(std::vector<struct Sphere> inHadrons);
  void randUnitR3(double &x, double &y, double &z);
  
  /** @brief Pointer to the InitialState module. */
  std::shared_ptr<InitialState> ini;

  //spatial dot product of t,x,y,z vectors
  double R3dot(double* x, double *y) { return x[1]*y[1] + x[2]*y[2] + x[3]*y[3]; }

private:
  double p_fake;
  int targZ, targA;
  int shiftHadrons, shiftType;

  // Allows the registration of the module so that it is available to be used by the Jetscape framework.
  static RegisterJetScapeModule<ColoredHadronizationSinglet> reg;

protected:
  static Pythia8::Pythia pythia;
  std::uniform_real_distribution<double> ZeroOneDistribution;
};

#endif // COLOREDHADRONIZATIONSINGLET_H
