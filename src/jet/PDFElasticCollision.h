#ifndef PDFELASTICCOLLISION_H
#define PDFELASTICCOLLISION_H

#include "PDFscat.h"
// #include "TVector3.h"
// #include <TRotation.h>
// #include <random>
// #include "FourVector.h"

class PDFElasticCollision{
    public:
        PDFElasticCollision();
		void setter(double max_energy0, double low_energy0, double energy_grid0, int index_type0, double eNucleon);
        bool elastic_kinematics(bool ProbailisticScattering, double deltaT, int &pid0, int &pid2, int &pid3, double (&pc0)[4], double (&pc2)[4], double (&pc3)[4], double &qt);
        double PerformLinearInterpolation(double EOriginal, int process_id, int type);
        double V[4];
    	//TVector3 iZ_Vector;
        PDFScat scattering_obj; 
    private:
        double obj_low_energy, obj_hig_energy, obj_grid_energy, obj_index;
        int pid_list[6]={1,2,3,-1,-2,-3};
	    std::default_random_engine generator;
	    // std::uniform_real_distribution<> uniform_rand; //between [0,1)              
};

#endif //PDFELASTICCOLLISION_H