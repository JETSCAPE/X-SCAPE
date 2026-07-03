#ifndef PDFELASTICCOLLISION_H
#define PDFELASTICCOLLISION_H

#include "PDFscat.h"
// #include "TVector3.h"
// #include <TRotation.h>
 #include <random>
// #include "FourVector.h"

class PDFElasticCollision{
    public:
        PDFElasticCollision();
		void setter(double max_energy0, double low_energy0, double energy_grid0, int index_type0);
        //bool elastic_kinematics(bool ProbailisticScattering, double deltaT, int &pid0, int &pid2, int &pid3, double (&pc0)[4], double (&pc2)[4], double (&pc3)[4], double &qt, double rho);
        double GetQhat_0(double E, int pid);
        double PerformLinearInterpolation(double EOriginal, int process_id, int type);
        void flavor(int &CT, int &KATT0, int &KATT2, int &KATT3,
                    unsigned int &max_color, unsigned int &color0,
                    unsigned int &anti_color0, unsigned int &color2,
                    unsigned int &anti_color2, unsigned int &color3,
                    unsigned int &anti_color3) ;
        void colljet22(int CT, double LambdaQCD, double p0[4], double p2[4], double p3[4], double p4[4], double &qt);
        void trans(double v[4], double p[4]);
        void transback(double v[4], double p[4]);
        void rotate(double px, double py, double pz, double pr[4], int icc);
        double V[4];
    	//TVector3 iZ_Vector;
        PDFScat scattering_obj; 
    private:
        double obj_low_energy, obj_hig_energy, obj_grid_energy, obj_index;
        //double nu_, Q2_;
        //int pid_list[6]={1,2,3,-1,-2,-3};
	    std::default_random_engine generator;
	    // std::uniform_real_distribution<> uniform_rand; //between [0,1)              
};

#endif //PDFELASTICCOLLISION_H