#ifndef TRANSVERSE_EPS09S_H
#define TRANSVERSE_EPS09S_H

#include "Pythia8/Pythia.h"
#include <string>
#include <fstream>
#include <iostream>
#include <algorithm>

class thickness_function {
private:
  int n_points_lins, n_points_linu;
  double ta[200];
  double tail_length;
  int a;

  static thickness_function* _instancea;
  static thickness_function* _instanceb;

protected:
  thickness_function() {}
  thickness_function(int a);
  ~thickness_function() {}

public:
  static thickness_function* getInstance(int a);
  double operator()(double s);
};

class fit_parameters {
private:
  double c_values[31][8][51][51][4];
  int order, pset;

  std::string get_filename(int eps_order, int eps_pset);
  static fit_parameters* _instance;

protected:
  fit_parameters() {}
  fit_parameters(int eps_order, int eps_pset);
  ~fit_parameters() {}

public:
  static fit_parameters* getInstance(int eps_order, int eps_pset);

  void interpolate_c_values(
      const double x,
      const double q,
      const int pset,
      double c[8][4]
  );
};

class EPS09s : public Pythia8::nPDF {

public:

  EPS09s(
      int idBeamIn = 2212,
      int iOrderIn = 1,
      int iSetIn = 1,
      Pythia8::PDFPtr protonPDFPtrIn = nullptr
  )
      : Pythia8::nPDF(idBeamIn, protonPDFPtrIn),
        ta(nullptr),
        fit_params(nullptr),
        iSet(1),
        iOrder(1),
        Tpos(0.0),
        sideLabel("")
  {
    init(iOrderIn, iSetIn);
  }

  void rUpdate(int id, double x, double Q2) override;

  void setErrorSet(int iSetIn) {
    iSet = iSetIn;
  }

  void setTpos(double TposIn) {
    Tpos = TposIn;
  }

  double getTpos() const {
    return Tpos;
  }

  void setSideLabel(const std::string& sideIn) {
    sideLabel = sideIn;
  }

private:
  thickness_function* ta;
  fit_parameters* fit_params;

  int iSet;
  int iOrder;
  double Tpos;
  std::string sideLabel;

  void init(int iOrderIn, int iSetIn);
};

#endif