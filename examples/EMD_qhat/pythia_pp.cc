// Standalone PYTHIA 8 pp reference for the JETSCAPE pp baseline.
//
//   pythia_pp out.dat.gz nEvents seed [mpi=1]
//
// pp at 5.02 TeV, HardQCD, pTHat 10-500 GeV with the same pTHat^4 bias as the
// JETSCAPE pp runs (PhaseSpace:bias2Selection, ref 10 GeV), default (Monash)
// tune, full event generation: ISR, FSR, MPI (unless mpi=0), hadronization,
// decays with tau0 < 10 mm/c as in JETSCAPE.  Final-state particles (no
// neutrinos) are written in the JETSCAPE final-state hadron format, with the
// bias weight as event weight, so that jetana_data reads them unchanged.
#include "Pythia8/Pythia.h"

#include <cstdio>
#include <cstdlib>
#include <string>

using namespace Pythia8;

int main(int argc, char **argv) {
  if (argc < 4) {
    std::fprintf(stderr, "usage: pythia_pp out.dat.gz nEvents seed [mpi=1]\n");
    return 1;
  }
  std::string out = argv[1];
  int nev = std::atoi(argv[2]);
  int seed = std::atoi(argv[3]);
  bool mpi = argc < 5 || std::atoi(argv[4]) != 0;

  Pythia pythia;
  pythia.readString("Beams:eCM = 5020.");
  pythia.readString("HardQCD:all = on");
  pythia.readString("PhaseSpace:pTHatMin = 10.");
  pythia.readString("PhaseSpace:pTHatMax = 500.");
  pythia.readString("PhaseSpace:bias2Selection = on");
  pythia.readString("PhaseSpace:bias2SelectionPow = 4.0");
  pythia.readString("PhaseSpace:bias2SelectionRef = 10.0");
  pythia.readString(std::string("PartonLevel:MPI = ") + (mpi ? "on" : "off"));
  pythia.readString("ParticleDecays:limitTau0 = on");
  pythia.readString("ParticleDecays:tau0Max = 10.0");
  pythia.readString("Random:setSeed = on");
  pythia.readString("Random:seed = " + std::to_string(seed));
  pythia.readString("Next:numberCount = 0");
  if (!pythia.init())
    return 1;

  FILE *fp = popen(("gzip -1 > '" + out + "'").c_str(), "w");
  if (!fp)
    return 1;
  std::fprintf(fp, "#\tJETSCAPE_FINAL_STATE\tv2\t|\tN\tpid\tstatus\tE\tPx\tPy\tPz\n");
  int written = 0;
  while (written < nev) {
    if (!pythia.next())
      continue;
    const Event &ev = pythia.event;
    int n = 0;
    for (int i = 0; i < ev.size(); i++) {
      int a = std::abs(ev[i].id());
      if (ev[i].isFinal() && a != 12 && a != 14 && a != 16)
        n++;
    }
    written++;
    std::fprintf(fp, "#\tEvent\t%d\tweight\t%.15g\tEPangle\t0\tN_hadrons\t%d\n", written,
                 pythia.info.weight(), n);
    int k = 0;
    for (int i = 0; i < ev.size(); i++) {
      int a = std::abs(ev[i].id());
      if (!ev[i].isFinal() || a == 12 || a == 14 || a == 16)
        continue;
      std::fprintf(fp, "%d %d 0 %.6g %.6g %.6g %.6g\n", k++, ev[i].id(), ev[i].e(), ev[i].px(),
                   ev[i].py(), ev[i].pz());
    }
  }
  pclose(fp);
  pythia.stat();
  return 0;
}
