#include "PDFElasticCollision.h"

#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>

/*
void Rectify1(TVector3 &v6){
	double threshold = 1e-10;
	if (std::abs(v6.X()) < threshold){ v6.SetX(0); }
	if (std::abs(v6.Y()) < threshold){ v6.SetY(0); } 
	if (std::abs(v6.Z()) < threshold){ v6.SetZ(0); }	
}


void Rectify_Momentum(double (&pc0)[4], const double E1, const double msq){
	double currentMagnitude = std::sqrt(pc0[1] * pc0[1] + pc0[2] * pc0[2] + pc0[3] * pc0[3] + msq);
	pc0[1] = E1*pc0[1]/currentMagnitude;
	pc0[2] = E1*pc0[2]/currentMagnitude; 
	pc0[3] = E1*pc0[3]/currentMagnitude;
}
*/

std::uniform_real_distribution<> uniform_rand(0.0,1.0);
PDFElasticCollision::PDFElasticCollision() {};

void PDFElasticCollision::setter(double max_energy0, double low_energy0, double energy_grid0,int index_type0) {
    obj_hig_energy = max_energy0;
    obj_low_energy = low_energy0;
    obj_grid_energy = energy_grid0;
    obj_index = index_type0;
    //nu_ = nu;
    //Q2_ = Q2;
    //iZ_Vector.SetXYZ(0,0,1); //won't be using root to do the rotations
    scattering_obj.setter(max_energy0, low_energy0, energy_grid0, index_type0);
	//scattering_obj.setter(0.2,0.2,0.1,102,0.5,0.5);
	//scattering_obj.setter(0.2,0.2,0.1,102,0.5,1);
    scattering_obj.initialize_samplers(); //the param here is energy of nucleon getting scattered off 
}
/*The MATTER module had pre-defination  of calculating the probability of scattring.*/

/*
bool PDFElasticCollision::elastic_kinematics(bool ProbailisticScattering, double deltaT, int &pid0, int &pid2, int &pid3, double (&pc0)[4], double (&pc2)[4], double (&pc3)[4], double &qt, double rho) {
    //pc2 is hole
    //pc3 is energetic parton after scattering
    //pc4 is recoil
    //double r0,r1,r2,r3,r4,r5;
    double EO,E1,E2,E3,E4,c23,c24,p1,p3,THETA2,THETA3,THETA4,PHI3,PHI4,PHI2,PHI23;
    //EO is the original energy
    double mc_sq;
    double TotalRate = 0;
    //double beta =  nu_ / sqrt(nu_ * nu_ + Q2_);
    //double gamma = sqrt(nu_ * nu_ + Q2_) / sqrt(Q2_);

    int parent_pid, pid_index;
    int daughter1_pid = -999;
    int daughter2_pid = -999;
    int hole_pid = -999;
    int energy_index;//, temp_index;
    int parton_type = -1;
    int proc; //heavy(1) or light(0)

    //TVector3 P1;
    //TVector3 P2;
    ////TVector3 P3;
    //TVector3 P4;
    //TRotation r;

    EO = pc0[0];
    //EO = gamma * (EO - beta * EO);//breit frame
    //find which E1 bin to use and switch over to that discretized value
    //energy_index = round((E1- obj_low_energy)/obj_grid_energy);
    E1 = scattering_obj.get_energy(EO, energy_index);

    parent_pid = pid0;
    pid_index = floor(uniform_rand(generator)*6); //ONLY UDS FOR NOW
    if (parent_pid == 21) {
		//r0 = scattering_obj.get_rate(energy_index,6); //has color abiguity for matter
        double r0g = PerformLinearInterpolation(EO, 6, 21, 1);
		double r1g = PerformLinearInterpolation(EO, 7, 21, 1);
        double r2u = PerformLinearInterpolation(EO, 8, 1, 1);
        double r2ubar = PerformLinearInterpolation(EO, 8, -1, 1);
        double r2d = PerformLinearInterpolation(EO, 8, 2, 1);
        double r2dbar = PerformLinearInterpolation(EO, 8, -2, 1);
        double r2s = PerformLinearInterpolation(EO, 8, 3, 1);
        double r2sbar = PerformLinearInterpolation(EO, 8, -3, 1);
        TotalRate = rho * (r0g + r1g + r2u + r2ubar + r2d + r2dbar + r2s + r2sbar);
        //std::cout<<"TotalRate gluon is "<<TotalRate<<std::endl;
        if (ProbailisticScattering){  
            if (exp(-TotalRate * deltaT) > uniform_rand(generator)){qt = 0; return false;}
        }  
        parton_type = 0;
    	std::discrete_distribution<int> distribution0{r0g,r1g,r2u,r2ubar,r2d,r2dbar,r2s,r2sbar};
        switch (distribution0(generator)) {
                case 0: //g g -> q qbar
                    scattering_obj.get_sample(21, energy_index,6,V, EO);
                    hole_pid      =  21;
                    daughter1_pid =  pid_list[pid_index];
                    daughter2_pid = -pid_list[pid_index];
                    break;

                case 1: //g g -> g g
                    scattering_obj.get_sample(21, energy_index,7,V, EO);
                    hole_pid      = 21;
                    daughter1_pid = 21;
                    daughter2_pid = 21;
                    break;

                case 2: //g u -> g u
                    scattering_obj.get_sample(1, energy_index,8,V, EO);
                    hole_pid      = 1;
                    daughter1_pid = 21;
                    daughter2_pid = 1;
                    break;
                case 3:
                    //g ubar -> g ubar
                    scattering_obj.get_sample(-1, energy_index,8,V, EO);
                    hole_pid      = -1;
                    daughter1_pid = 21;
                    daughter2_pid = -1;
                    break;
                case 4:
                    //g d -> g d
                    scattering_obj.get_sample(2, energy_index,8,V, EO);
                    hole_pid      = 2;
                    daughter1_pid = 21;
                    daughter2_pid = 2;
                    break;
                case 5:
                    //g dbar -> g dbar
                    scattering_obj.get_sample(-2, energy_index,8,V, EO);
                    hole_pid      = -2;
                    daughter1_pid = 21;
                    daughter2_pid = -2;
                    break;
                case 6:
                    //g s -> g s
                    scattering_obj.get_sample(3, energy_index,8,V, EO);
                    hole_pid      = 3;
                    daughter1_pid = 21;
                    daughter2_pid = 3;
                    break;
                case 7:
                    //g sbar -> g sbar
                    scattering_obj.get_sample(-3, energy_index,8,V, EO);
                    hole_pid      = -3;
                    daughter1_pid = 21;
                    daughter2_pid = -3;
                    break;
                default:
                    //never gets here
                    //raise error?
                    break;
        }
    }   
    
    else if (abs(parent_pid) >= 1 and abs(parent_pid) <= 3) {
		double r0 = 0; //PerformLinearInterpolation(EO,0,1);//s-channel not considered
		double r1bar = PerformLinearInterpolation(EO,1,-parent_pid,1);
		double r2 = PerformLinearInterpolation(EO,2,parent_pid,1);
		double r3bar = PerformLinearInterpolation(EO,3,-parent_pid,1);
		double r4 = PerformLinearInterpolation(EO, 4, 21, 1);
		double r5a = PerformLinearInterpolation(EO,5,(parent_pid < 0 ? -1 : 1) * (((abs(parent_pid)) % 3) + 1),1);
        double r5b = PerformLinearInterpolation(EO,5,(parent_pid < 0 ? -1 : 1) * (((abs(parent_pid)+1) % 3) + 1),1);
        TotalRate = rho * (r0 + r1bar + r2 + r3bar + r4 + r5a + r5b);
        //std::cout<<"TotalRate quark for pid  "<<parent_pid<<" is "<<TotalRate<<std::endl;
        if (ProbailisticScattering){  
            if (exp(-TotalRate * deltaT) > uniform_rand(generator)){qt = 0; return false;}
        }  
        parton_type = 0;
    	std::discrete_distribution<int> distribution0{r0,r1bar,r2,r3bar,r4,r5a,r5b};
        switch (distribution0(generator)) {
                case 0: //q1 q1bar -> q2 q2bar
                    //scattering_obj.get_sample(energy_index,0,V, EO);
                    //hole_pid      = -parent_pid;
                    //do { //keep sampling until q2 != q1
                    //    pid_index = floor(uniform_rand(generator)*6);
                    //    daughter1_pid = pid_list[pid_index];
                    //} while (daughter1_pid == parent_pid);                
                    //daughter2_pid = -daughter1_pid;
                    break;

                case 1: //q1 q1bar -> q1 q1bar
                    scattering_obj.get_sample(-parent_pid, energy_index, 1, V, EO);
                    hole_pid      = -parent_pid;
                    daughter1_pid =  parent_pid;
                    daughter2_pid = -parent_pid;
                    break;

                case 2: //q1 q1 -> q1 q1
                    scattering_obj.get_sample(parent_pid, energy_index, 2, V, EO);
                    hole_pid      = parent_pid;
                    daughter1_pid = parent_pid;
                    daughter2_pid = parent_pid;
                    break;

                case 3: //q1 q1bar -> g g
                    scattering_obj.get_sample(-parent_pid, energy_index, 3, V, EO);
                    hole_pid      = -parent_pid;
                    daughter1_pid = 21;
                    daughter2_pid = 21;
                    break;

                case 4: //q1 g -> q1 g
                    scattering_obj.get_sample(21, energy_index,4,V, EO);
                    hole_pid      = 21;
                    daughter1_pid = parent_pid;
                    daughter2_pid = 21;
                    break;

                case 5: //q1 q2 -> q1 q2
                    daughter2_pid = (parent_pid < 0 ? -1 : 1) * (((abs(parent_pid)) % 3) + 1);
                    scattering_obj.get_sample(daughter2_pid,energy_index,5,V, EO);
                    daughter1_pid = parent_pid;
                    //do { //keep sampling until q2 != q1 or q1bar
                    //    pid_index = floor(uniform_rand(generator)*6);
                    //    daughter2_pid = pid_list[pid_index];
                    //} while (daughter2_pid == parent_pid || daughter2_pid == -parent_pid);
                    hole_pid=daughter2_pid;     
                    break;
                case 6: //q1 q2 -> q1 q2
                    daughter2_pid = (parent_pid < 0 ? -1 : 1) * (((abs(parent_pid)+1) % 3) + 1);
                    scattering_obj.get_sample(daughter2_pid,energy_index,5,V, EO);
                    daughter1_pid = parent_pid;
                    hole_pid=daughter2_pid;     
                    break;
                default:
                    //never gets here
                    break;
        }
    }
     //Heavy scattering part will be added here
    else if (abs(parent_pid)==4) {
        
        mc_sq = 1.6129;

		r0 = PerformLinearInterpolation(EO,9,1);
		r1 = PerformLinearInterpolation(EO,10,1);
        TotalRate = 0.0; //r0 + r1;
        if (ProbailisticScattering){  
            if (exp(-TotalRate * deltaT) > uniform_rand(generator)){qt = 0; return false;}
        }
        parton_type = 1;
    	std::discrete_distribution<int> distribution0{r0,r1};
        switch (distribution0(generator)) {
                case 0: //TODO: should be q1 q1bar -> q2 q2bar?? 
                    scattering_obj.get_sample(energy_index,9,V, EO);
                    hole_pid      = pid_list[pid_index];
                    daughter1_pid = parent_pid;
                    daughter2_pid = pid_list[pid_index];
                    break;

                case 1: //TODO: should be q1 q1bar -> q1 q1bar?
                    scattering_obj.get_sample(energy_index,10,V, EO);
                    hole_pid      = 21;
                    daughter1_pid = parent_pid;
                    daughter2_pid = 21;
                    break;

                default:
                    //never gets here
                    break;
        }
         return false;
    }
        
    else {
        
        mc_sq = 17.4724;

		r0 = PerformLinearInterpolation(EO,11,1);
		r1 = PerformLinearInterpolation(EO,12,1);
        TotalRate = 0.0; //r0 + r1;
        if (ProbailisticScattering){  
            if (exp(-TotalRate * deltaT) > uniform_rand(generator)){qt = 0; return false;}
        }
        parton_type = 1;
    	std::discrete_distribution<int> distribution0{r0,r1};
        switch (distribution0(generator)) {
            case 0: //TODO: should be q1 q1bar -> q2 q2bar?? 
                    scattering_obj.get_sample(energy_index,11,V, EO);
                    hole_pid      = pid_list[pid_index];
                    daughter1_pid = parent_pid;
                    daughter2_pid = pid_list[pid_index];
                    break;

                case 1: //TODO: should be q1 q1bar -> q1 q1bar?
                    scattering_obj.get_sample(energy_index,12,V, EO);
                    hole_pid      = 21;
                    daughter1_pid = parent_pid;
                    daughter2_pid = 21;
                    break;

                default:
                    //never gets here
                    break;
        }
       return false;
    }

    if (parton_type == 0) { //light parton
        THETA2 = V[0];
        THETA3 = V[1];
        PHI23  = V[2];
        E3     = V[3];

		c23 = cos(THETA2) * cos(THETA3) + sin(THETA2) * sin(THETA3) * cos(PHI23);
		
        //E2 = (E1*E3*(1-cos(THETA3)))/(E1*(1-cos(THETA2))-E3*(1-c23));
        E2 = (EO * E3 * (1 - cos(THETA3))) / (EO * (1 - cos(THETA2)) - E3 * (1 - c23));
		//E4 = E1+E2-E3;
        E4 = EO + E2 - E3;
		
        do{
            PHI2 = 2 * M_PI * uniform_rand(generator);
        } while (E2 * sin(THETA2) * cos(PHI2) < 0.1 || E2 * sin(THETA2) * sin(PHI2) < 0.1); //avoid numerical instability
		PHI3 = PHI2 - PHI23;

		//THETA4= acos((E1+E2*cos(THETA2)-E3*cos(THETA3))/E4);
        THETA4= acos((EO + E2 * cos(THETA2) - E3 * cos(THETA3)) / E4);

		double value = (E2 * sin(THETA2) * cos(PHI2) - E3 * sin(THETA3) * cos(PHI3)) / (E4 * sin(THETA4));
		value = std::max(-1.0, std::min(1.0, value)); // Clamping the value
		PHI4 = acos(value);

        //std::cout<<"LQ E1 "<<E1<<" E2 "<<E2<<" E3 "<<E3<<" E4 "<<E4<<std::endl;
        pc2[0] = E2;
        pc2[1] = E2 * sin(THETA2) * cos(PHI2);
        pc2[2] = E2 * sin(THETA2) * sin(PHI2);
        pc2[3] = E2 * cos(THETA2);
        //pc2[0] = gamma * (E2 + beta * E2 * cos(THETA2));
        //pc2[3] = gamma * (E2 * cos(THETA2) + beta * E2);
        pid2 = hole_pid;

        pc0[0] = E3;
        pc0[1] = E3 * sin(THETA3) * cos(PHI3);
        pc0[2] = E3 * sin(THETA3) * sin(PHI3);
        pc0[3] = E3 * cos(THETA3);
        //pc0[0] = gamma * (E3 + beta * E3 * cos(THETA3));
        //pc0[3] = gamma * (E3 * cos(THETA3) + beta * E3);
        pid0 = daughter1_pid;

        pc3[0] = E4;
        pc3[1] = E4 * sin(THETA4) * cos(PHI4);
        pc3[2] = E4 * sin(THETA4) * sin(PHI4);
        pc3[3] = E4 * cos(THETA4);   
        //pc3[0] = gamma * (E4 + beta * E4 * cos(THETA4));
        //pc3[3] = gamma * (E4 * cos(THETA4) + beta * E4);
        pid3 = daughter2_pid;
        
        //Will perform rotation outside this class;

        //P1.SetXYZ(pc0[1],pc0[2],pc0[3]);
        //P2.SetXYZ(E2*sin(THETA2)*cos(PHI2),E2*sin(THETA2)*sin(PHI2),E2*cos(THETA2));
		//P3.SetXYZ(E3*sin(THETA3)*cos(PHI3),E3*sin(THETA3)*sin(PHI3),E3*cos(THETA3));
		//P4.SetXYZ(E4*sin(THETA4)*cos(PHI4),E4*sin(THETA4)*sin(PHI4),E4*cos(THETA4));
		
		//double s0=P1.Angle(iZ_Vector);
		//TVector3 s1;
		//s1=P1.Cross(iZ_Vector);
        
		//r.Rotate(s0,s1);
        //P2=r*P2;
        //P3=r*P3;
        //P4=r*P4;
		
        //std::cout<<"LQ  after one roatation E1 "<<E1<<" E2 "<<E2<<" E3 "<<E3<<" E4 "<<E4<<std::endl;
    }
    else if (parton_type==1) { //heavy
        // Will implement later
    	E2     = V[0];
        THETA2 = V[1];
        THETA4 = V[2];
        PHI4   = V[3];

		//std::cout<<"E2 "<<E2<<" theta2 "<<THETA2<<" theta4 "<<THETA4<<" phi4 before "<<PHI4<<std::endl;
        c24 = sin(THETA2)*sin(THETA4)*cos(PHI4)+cos(THETA2)*cos(THETA4);

		if (E1*E1<mc_sq){return -1;}
		Rectify_Momentum(pc0,E1,mc_sq);

		//std::cout<<"E1 "<<E1*E1<<" mcsq "<<mc_sq<<std::endl;
        p1 = sqrt(E1*E1-mc_sq);
        E4 = (E1*E2-p1*E2*cos(THETA2))/(E1-p1*cos(THETA4)+E2-E2*c24);
        E3 = E1+E2-E4;
		if (E3*E3<mc_sq){return -1;}
        p3 = sqrt(E3*E3-mc_sq);

		//std::cout<<"p1 "<<p1<<" E4 "<<E4<<" E3 "<<E3<<" p3 after "<<p3<<std::endl;

		double value = (p1+E2*cos(THETA2)-E4*cos(THETA4))/p3;
		if (value>1 || value<-1) {std::cout<<"value cos "<<value<<std::endl;}
		value = std::max(-1.0, std::min(1.0, value));
        THETA3 = acos(value);

		value = -(E4*sin(PHI4)*sin(THETA4))/(p3*sin(THETA3));
		if (value>1 || value<-1){std::cout<<"value sin "<<value<<std::endl;}
		value = std::max(-1.0, std::min(1.0, value));
        PHI3 = asin(value);

        PHI2 = 2*M_PI*uniform_rand(generator);     
        r.RotateZ(PHI2);
        //std::cout<<"HQ E1 "<<E1<<" E2 "<<E2<<" E3 "<<E3<<" E4 "<<E4<<std::endl;
        P1.SetXYZ(pc0[1],pc0[2],pc0[3]);
        P2.SetXYZ(E2*sin(THETA2),0,E2*cos(THETA2));
        P3.SetXYZ(p3*sin(THETA3)*cos(PHI3),p3*sin(THETA3)*sin(PHI3),p3*cos(THETA3));
        P4.SetXYZ(E4*sin(THETA4)*cos(PHI4),E4*sin(THETA4)*sin(PHI4),E4*cos(THETA4));

        P2=r*P2;
        P3=r*P3;
        P4=r*P4;
 		Rectify1(P2);
 		Rectify1(P3);
 		Rectify1(P4);	
        //std::cout<<"HQ after one rotation E1 "<<E1<<" E2 "<<E2<<" E3 "<<E3<<" E4 "<<E4<<std::endl;
        
    }
    else { return false; }

    return true;
}
*/

double PDFElasticCollision::GetQhat_0(double E, int pid){
    double qhat = 0;
    int energy_index;
    //energy_index = std::min( std::max(0, int(round((E1 - obj_low_energy) / obj_grid_energy))), int(round((obj_hig_energy - obj_low_energy) / obj_grid_energy)));
    //std::cout<<"obj_low_energy "<<obj_low_energy<<" obj_grid_enrgy "<<obj_grid_energy<<std::endl;
    //int energy_index = round((E- obj_low_energy)/obj_grid_energy);
    E = scattering_obj.get_energy(E, energy_index);
    //std::cout<<"energy index is "<<energy_index<<" energy is "<<E<<std::endl;
    //int temp_index = round((T-obj_low_temp)/obj_grid_temp);
    //int temp_index = std::min(std::max(0, int(round((T - obj_low_temp) / obj_grid_temp))), int(round((obj_hig_temp - obj_low_temp) / obj_grid_temp)));
    if (abs(pid) == 21){
        qhat = PerformLinearInterpolation(energy_index, 6, 2)
            + PerformLinearInterpolation(energy_index, 7, 2)
            + PerformLinearInterpolation(energy_index, 8, 2);
    }
    else if (abs(pid) < 3){
        qhat = 0
            + PerformLinearInterpolation(energy_index, 1, 2)
            + PerformLinearInterpolation(energy_index, 2, 2)
            + PerformLinearInterpolation(energy_index, 3, 2)
            + PerformLinearInterpolation(energy_index, 4, 2)
            + PerformLinearInterpolation(energy_index, 5, 2); 
    }
    else{
        qhat = 0;
    }
    return qhat;
}


double PDFElasticCollision::PerformLinearInterpolation(double EOriginal, int process_id, int type){
    int EIndex, EIndex_; //index where sampler is trained
    double E, E_;
    double ERoundedOff = scattering_obj.get_energy(EOriginal, EIndex);
    //std::cout<<std::setprecision(4)<<"ERoundedOff "<<ERoundedOff<<"  EOriginal "<<EOriginal<<std::endl;
    if (ERoundedOff <= EOriginal){
        EIndex_ = std::min(EIndex + 1, scattering_obj.energy_index_range);
    }
    
    else{
        EIndex_ = EIndex;
        EIndex = std::max(0, EIndex_ - 1);
    }
    E_ = scattering_obj.get_energy(EIndex_);
    E = scattering_obj.get_energy(EIndex);
    double ValE1, ValE2;
    double ValFinal;
    if (type == 1){
        //Rate interpolation
        ValE1 = scattering_obj.get_rate(EIndex, process_id);
        ValE2 = scattering_obj.get_rate(EIndex_, process_id);
        //std::cout<<"Rate ValE1 "<<ValE1<<" ValE2 "<<ValE2<<std::endl;
    }
    else{
        //Qhat interpolation
        ValE1 = scattering_obj.get_qhat(EIndex, process_id);
        ValE2 = scattering_obj.get_qhat(EIndex_, process_id);
        //std::cout<<"Qhat ValE1 "<<ValE1<<" ValE2 "<<ValE2<<std::endl;
    }
    if (EOriginal - ERoundedOff < 1e-4){
        //Takes care when out of range
        return ValE1;
    }
    ValFinal = (ValE2 - ValE1) * (EOriginal - E) / (E_ - E) + ValE1;
    return ValFinal;
}

inline double snap_to_zero(double x, double eps = 1e-8) {
    return (std::abs(x) < eps) ? 0.0 : x;
}


void rotate_vector(double (&V)[4], int rot_type, double& theta, double& phi){
    double p1x, p1y, p1z;
    if (rot_type == 1){
        theta = snap_to_zero(atan2(sqrt(V[2]*V[2] + V[1] * V[1]), V[3]));
        phi = snap_to_zero(atan2(V[2], V[1]));
    }
    //std::cout<<"theta "<<theta<<" phi "<<phi<<std::endl;
    p1x = V[1];
    p1y = V[2];
    p1z = V[3];
    if (rot_type == 1){
        //rotate around z axis clockwise(if seen from +z axis) by phi [[cos(phi) sin(phi) 0][-sin(phi) cos(phi) 0][0 0 1]]
        V[1] = p1x * cos(phi) + p1y * sin(phi);
        V[2] = -p1x * sin(phi) + p1y * cos(phi);
        V[3] = p1z;
        p1x = V[1];
        p1y = V[2];
        p1z = V[3];
        //rotate around y axis clockwise (if seen from +y axis) by theta [[cos(theta) 0 -sin(theta)][0 1 0][sin(theta) 0 cos(theta)]]
        V[1] = snap_to_zero(p1x * cos(theta) - p1z * sin(theta));
        V[2] = snap_to_zero(p1y);
        V[3] = snap_to_zero(p1x * sin(theta) + p1z * cos(theta));
    }
    if (rot_type == -1){
        //rotate around y axis counter-clockwise (if seen from +y axis) by theta [[cos(theta) 0 sin(theta)][0 1 0][-sin(theta) 0 cos(theta)]]
        V[1] = p1x * cos(theta) + p1z * sin(theta);
        V[2] = p1y;
        V[3] = -p1x * sin(theta) + p1z * cos(theta);
        p1x = V[1];
        p1y = V[2];
        p1z = V[3];
        //rotate around z axis counter-clockwise(if seen from +z axis) by phi [[cos(phi) -sin(phi) 0][sin(phi) cos(phi) 0][0 0 1]]
        V[1] = snap_to_zero(p1x * cos(phi) - p1y * sin(phi));
        V[2] = snap_to_zero(p1x * sin(phi) + p1y * cos(phi));
        V[3] = snap_to_zero(p1z);
    }
}

void PDFElasticCollision::flavor(int &CT, int &KATT0, int &KATT2, int &KATT3,
                    unsigned int &max_color, unsigned int &color0,
                    unsigned int &anti_color0, unsigned int &color2,
                    unsigned int &anti_color2, unsigned int &color3,
                    unsigned int &anti_color3) {

  int vb[7] = {0};
  int b = 0;
  int KATT00 = KATT0;
  unsigned int backup_color0 = color0;
  unsigned int backup_anti_color0 = anti_color0;

  vb[1] = 1;
  vb[2] = 2;
  vb[3] = 3;
  vb[4] = -1;
  vb[5] = -2;
  vb[6] = -3;

  if (KATT00 == 21) { //.....for gluon
    double R1 = 16.0; // gg->gg DOF_g
    double R2 = 0.0;  // gg->qqbar don't consider this channel in eMATTER
    double R3 = 6.0 * 6 * 4 / 9; // gq->gq or gqbar->gqbar flavor*DOF_q*factor
    double R0 = R1 + R3;

    double a = uniform_rand(generator);

    if (a <= R1 / R0) { // gg->gg
      CT = 1;
      KATT3 = 21;
      KATT2 = 21;
      //KATT0=KATT0;
      max_color++;
      color0 = backup_color0;
      anti_color0 = max_color;
      max_color++;
      color2 = anti_color0;
      anti_color2 = max_color;
      color3 = backup_anti_color0;
      anti_color3 = max_color;
    } else { // gq->gq
      CT = 3;
      b = floor(uniform_rand(generator) * 6 + 1);
      if (b == 7)
        b = 6;
      KATT3 = vb[b];
      KATT2 = KATT3;
      //KATT0=KATT0;
      if (KATT3 > 0) { // gq->gq
        max_color++;
        color0 = backup_color0;
        anti_color0 = max_color;
        color2 = max_color;
        anti_color2 = 0;
        color3 = backup_anti_color0;
        anti_color3 = 0;
      } else { // gqbar->gqbar
        max_color++;
        color0 = max_color;
        anti_color0 = backup_anti_color0;
        color2 = 0;
        anti_color2 = max_color;
        color3 = 0;
        anti_color3 = backup_color0;
      }
    }
  } else if (abs(KATT00) == 4) { // for charm quarks
    double R1 = 6.0 * 6 * 4 / 9; // Qq->Qq
    double R2 = 16.0;            // Qg->Qg DOF_ag
    double R00 = R1 + R2;

    double a = uniform_rand(generator);

    if (a <= R2 / R00) { // Qg->Qg
      CT = 12;
      KATT3 = 21;
      KATT2 = 21;
      if (KATT00 > 0) { // Qg->Qg
        max_color++;
        color0 = max_color;
        anti_color0 = 0;
        max_color++;
        color2 = max_color;
        anti_color2 = color0;
        color3 = max_color;
        anti_color3 = backup_color0;
      } else { // Qbarg->Qbarg
        max_color++;
        color0 = 0;
        anti_color0 = max_color;
        max_color++;
        color2 = anti_color0;
        anti_color2 = max_color;
        color3 = backup_anti_color0;
        anti_color3 = max_color;
      }
    } else { // Qq->Qq
      CT = 11;
      b = floor(uniform_rand(generator) * 6 + 1);
      if (b == 7)
        b = 6;
      KATT3 = vb[b];
      KATT2 = KATT3;
      if (KATT00 > 0 && KATT2 > 0) { // qq->qq
        max_color++;
        color0 = max_color;
        anti_color0 = 0;
        color2 = backup_color0;
        anti_color2 = 0;
        color3 = max_color;
        anti_color3 = 0;
      } else if (KATT00 > 0 && KATT2 < 0) { //qqbar->qqbar
        max_color++;
        color0 = max_color;
        anti_color0 = 0;
        color2 = 0;
        anti_color2 = max_color;
        color3 = 0;
        anti_color3 = backup_color0;
      } else if (KATT00 < 0 && KATT2 > 0) { //qbarq->qbarq
        max_color++;
        color0 = 0;
        anti_color0 = max_color;
        color2 = max_color;
        anti_color2 = 0;
        color3 = backup_anti_color0;
        anti_color3 = 0;
      } else { //qbarqbar->qbarqbar
        max_color++;
        color0 = 0;
        anti_color0 = max_color;
        color2 = 0;
        anti_color2 = backup_anti_color0;
        color3 = 0;
        anti_color3 = max_color;
      }
    }

  } else {                       //.....for quark and antiquark (light)
    double R3 = 16.0;            // qg->qg DOF_g
    double R4 = 4.0 * 6 * 4 / 9; // qq'->qq' scatter with other species
    double R5 = 1.0 * 6 * 4 / 9; // qq->qq scatter with itself
    double R6 =
        0.0; // qqbar->q'qbar' to other final state species, don't consider in eMATTER
    double R7 = 1.0 * 6 * 4 / 9; // qqbar->qqbar scatter with its anti-particle
    double R8 = 0.0;             // qqbar->gg don't consider in eMATTER
    double R00 = R3 + R4 + R5 + R7;

    double a = uniform_rand(generator);
    if (a <= R3 / R00) { // qg->qg
      CT = 13;
      KATT3 = 21;
      KATT2 = 21;
      //KATT0=KATT0;
      if (KATT00 > 0) { // qg->qg
        max_color++;
        color0 = max_color;
        anti_color0 = 0;
        max_color++;
        color2 = max_color;
        anti_color2 = color0;
        color3 = max_color;
        anti_color3 = backup_color0;
      } else { // qbarg->qbarg
        max_color++;
        color0 = 0;
        anti_color0 = max_color;
        max_color++;
        color2 = anti_color0;
        anti_color2 = max_color;
        color3 = backup_anti_color0;
        anti_color3 = max_color;
      }
    } else if (a <= (R3 + R4) / R00) { // qq'->qq'
      CT = 4;
      do {
        b = floor(uniform_rand(generator) * 6 + 1);
        if (b == 7)
          b = 6;
        KATT3 = vb[b];
      } while (KATT3 == KATT0 || KATT3 == -KATT0);
      KATT2 = KATT3;
      //KATT0=KATT0;
      if (KATT00 > 0 && KATT2 > 0) { // qq->qq
        max_color++;
        color0 = max_color;
        anti_color0 = 0;
        color2 = backup_color0;
        anti_color2 = 0;
        color3 = max_color;
        anti_color3 = 0;
      } else if (KATT00 > 0 && KATT2 < 0) { //qqbar->qqbar
        max_color++;
        color0 = max_color;
        anti_color0 = 0;
        color2 = 0;
        anti_color2 = max_color;
        color3 = 0;
        anti_color3 = backup_color0;
      } else if (KATT00 < 0 && KATT2 > 0) { //qbarq->qbarq
        max_color++;
        color0 = 0;
        anti_color0 = max_color;
        color2 = max_color;
        anti_color2 = 0;
        color3 = backup_anti_color0;
        anti_color3 = 0;
      } else { //qbarqbar->qbarqbar
        max_color++;
        color0 = 0;
        anti_color0 = max_color;
        color2 = 0;
        anti_color2 = backup_anti_color0;
        color3 = 0;
        anti_color3 = max_color;
      }
    } else if (a <= (R3 + R4 + R5) / R00) { // scatter with itself
      CT = 5;
      KATT3 = KATT0;
      KATT2 = KATT0;
      //KATT0=KATT0;
      if (KATT00 > 0) { // qq->qq
        max_color++;
        color0 = max_color;
        anti_color0 = 0;
        color2 = backup_color0;
        anti_color2 = 0;
        color3 = max_color;
        anti_color3 = 0;
      } else { //qbarqbar->qbarqbar
        max_color++;
        color0 = 0;
        anti_color0 = max_color;
        color2 = 0;
        anti_color2 = backup_anti_color0;
        color3 = 0;
        anti_color3 = max_color;
      }
    } else { // scatter with its anti-particle
      CT = 7;
      KATT3 = -KATT0;
      KATT2 = KATT3;
      //KATT0=KATT0;
      if (KATT00 > 0) { //qqbar->qqbar
        max_color++;
        color0 = max_color;
        anti_color0 = 0;
        color2 = 0;
        anti_color2 = max_color;
        color3 = 0;
        anti_color3 = backup_color0;
      } else { //qbarq->qbarq
        max_color++;
        color0 = 0;
        anti_color0 = max_color;
        color2 = max_color;
        anti_color2 = 0;
        color3 = backup_anti_color0;
        anti_color3 = 0;
      }
    }
  }
}

void PDFElasticCollision::colljet22(int CT, double LambdaQCD, double p0[4], double p2[4], double p3[4], double p4[4], double &qt) {
  //
  //    p0 initial jet momentum, output to final momentum
  //    p2 final thermal momentum,p3 initial termal energy
  //
  //    amss=sqrt(abs(p0(4)**2-p0(1)**2-p0(2)**2-p0(3)**2))
  //
  //************************************************************
  //p4[1] = p0[1];
  //p4[2] = p0[2];
  //p4[3] = p0[3];
  //p4[0] = p0[0];
  //************************************************************

  //    transform to local comoving frame of the fluid
  //  cout << endl;
  //  cout << "flow  "<< v0[1] << " " << v0[2] << " " << v0[3] << " "<<" Elab " << p0[0] << endl;

  //trans(v0, p0); //no need for boost in eA generatr
  //  cout << p0[0] << " " << sqrt(qhat0ud) << endl;

  //  cout << sqrt(pow(p0[1],2)+pow(p0[2],2)+pow(p0[3],2)) << " " << p0[1] << " " << p0[2] << " " << p0[3] << endl;

  //************************************************************
  //trans(v0, p4); // no need for boost in eA generatr
  //************************************************************

  //    sample the medium parton thermal momentum in the comoving frame
  double var = 0.8; //GeV^2
  double xw; //random variable for parton from medium
  double razim;
  double rcos;
  double rsin;

  double ss;
  double tmin;
  double tmid;
  double tmax;

  double rant;
  double tt;

  double uu;
  double ff = 0.0;
  double rank;

  double mmax;
  double msq = 0.0;

  double f1;
  double f2;

  double p0ex[4] = {0.0};
  double vc[4] = {0.0};

  int ct1_loop, ct2_loop, flag1, flag2;

  flag1 = 0;
  flag2 = 0;

  //  Initial 4-momentum of jet
  //
  //************************************************************
  p4[1] = p0[1];
  p4[2] = p0[2];
  p4[3] = p0[3];
  p4[0] = p0[0];
  //************************************************************

  int ic = 0;

  ct1_loop = 0;
  do {
    ct1_loop++;
    if(flag2 == 1 || ct1_loop > 1e6){
       flag1 = 1;
       break;
    }
    ct2_loop = 0;
    do {
      ct2_loop++;
      if(ct2_loop > 1e6){
         flag2 = 1;
         break;
      }
      xw = 4.0 * uniform_rand(generator);
      razim = 2.0 * M_PI * uniform_rand(generator);
      rcos = 1.0 - 2.0 * uniform_rand(generator);
      rsin = sqrt(1.0 - rcos * rcos);
      //
      p2[0] = xw;
      p2[3] = p2[0] * rcos;
      p2[1] = p2[0] * rsin * cos(razim);
      p2[2] = p2[0] * rsin * sin(razim);

      //
      //    cms energy
      //
      ss =
          2.0 * (p0[0] * p2[0] - p0[1] * p2[1] - p0[2] * p2[2] - p0[3] * p2[3]);

      //	if(ss.lt.2.d0*qhat0ud) goto 14

      tmin = LambdaQCD * LambdaQCD;
      tmid = ss / 2.0;
      tmax = ss - LambdaQCD * LambdaQCD;

      //    use (s^2+u^2)/(t+qhat0ud)^2 as scattering cross section in the
      //
      rant = uniform_rand(generator);
      tt = rant * ss;

      //		ic+=1;
      //		cout << p0[0] << "  " << p2[0] <<  endl;
      //		cout << tt << "  " << ss <<  "" << qhat0ud <<endl;
      //		cout << ic << endl;

    } while ((tt < tmin) || (tt > (tmax)) || abs(p2[1]) < 0.1 || abs(p2[2]) < 0.1 ||
             abs(p2[3]) < 0.1);

    f1 = pow(p2[0], 3) * exp(-p2[0]*p2[0]/var) / 0.293;
    f2 = pow(p2[0], 3) * exp(-p2[0]*p2[0]/var) / 0.293;

    uu = ss - tt;

    if (CT == 1) {
      ff = f1;
      mmax =
          4.0 / pow(ss, 2) *
          (3.0 - tmin * (ss - tmin) / pow(ss, 2) +
           (ss - tmin) * ss / pow(tmin, 2) + tmin * ss / pow((ss - tmin), 2));
      msq = pow((1.0 / p0[0] / p2[0] / 2.0), 2) *
            (3.0 - tt * uu / pow(ss, 2) + uu * ss / pow(tt, 2) +
             tt * ss / pow(uu, 2)) /
            mmax;
    }

    if (CT == 2) {
      ff = f1;
      mmax = 4.0 / pow(ss, 2) *
             (4.0 / 9.0 * (pow(tmin, 2) + pow((ss - tmin), 2)) / tmin /
                  (ss - tmin) -
              (pow(tmin, 2) + pow((ss - tmin), 2)) / pow(ss, 2));
      msq = pow((1.0 / p0[0] / p2[0] / 2.0), 2) *
            (4.0 / 9.0 * (pow(tt, 2) + pow(uu, 2)) / tt / uu -
             (pow(tt, 2) + pow(uu, 2)) / pow(ss, 2)) /
            (mmax + 4.0);
    }

    if (CT == 3) {
      ff = f2;
      if (((pow(ss, 2) + pow((ss - tmin), 2)) / pow(tmin, 2) +
           4.0 / 9.0 * (pow(ss, 2) + pow((ss - tmin), 2)) / ss / (ss - tmin)) >
          ((pow(ss, 2) + pow((ss - tmax), 2) / pow(tmax, 2) +
            4.0 / 9.0 * (pow(ss, 2) + pow((ss - tmax), 2)) / ss /
                (ss - tmax)))) {
        mmax =
            4.0 / pow(ss, 2) *
            ((pow(ss, 2) + pow((ss - tmin), 2)) / pow(tmin, 2) +
             4.0 / 9.0 * (pow(ss, 2) + pow((ss - tmin), 2)) / ss / (ss - tmin));
      } else {
        mmax =
            4.0 / pow(ss, 2) *
            ((pow(ss, 2) + pow((ss - tmax), 2)) / pow(tmax, 2) +
             4.0 / 9.0 * (pow(ss, 2) + pow((ss - tmax), 2)) / ss / (ss - tmax));
      }
      //
      msq = pow((1.0 / p0[0] / p2[0] / 2.0), 2) *
            ((pow(ss, 2) + pow(uu, 2)) / pow(tt, 2) +
             4.0 / 9.0 * (pow(ss, 2) + pow(uu, 2)) / ss / uu) /
            mmax;
    }

    if (CT == 13) {
      ff = f1;

      if (((pow(ss, 2) + pow((ss - tmin), 2)) / pow(tmin, 2) +
           4.0 / 9.0 * (pow(ss, 2) + pow((ss - tmin), 2)) / ss / (ss - tmin)) >
          ((pow(ss, 2) + pow((ss - tmax), 2) / pow(tmax, 2) +
            4.0 / 9.0 * (pow(ss, 2) + pow((ss - tmax), 2)) / ss /
                (ss - tmax)))) {
        mmax =
            4.0 / pow(ss, 2) *
            ((pow(ss, 2) + pow((ss - tmin), 2)) / pow(tmin, 2) +
             4.0 / 9.0 * (pow(ss, 2) + pow((ss - tmin), 2)) / ss / (ss - tmin));
      } else {
        mmax =
            4.0 / pow(ss, 2) *
            ((pow(ss, 2) + pow((ss - tmax), 2)) / pow(tmax, 2) +
             4.0 / 9.0 * (pow(ss, 2) + pow((ss - tmax), 2)) / ss / (ss - tmax));
      }
      //
      msq = pow((1.0 / p0[0] / p2[0] / 2.0), 2) *
            ((pow(ss, 2) + pow(uu, 2)) / pow(tt, 2) +
             4.0 / 9.0 * (pow(ss, 2) + pow(uu, 2)) / ss / uu) /
            mmax;
    }

    if (CT == 4) {
      ff = f2;
      mmax = 4.0 / pow(ss, 2) *
             (4.0 / 9.0 * (pow(ss, 2) + pow((ss - tmin), 2)) / pow(tmin, 2));
      msq = pow((1.0 / p0[0] / p2[0] / 2.0), 2) *
            (4.0 / 9.0 * (pow(ss, 2) + pow(uu, 2)) / pow(tt, 2)) / mmax;
    }

    if (CT == 5) {
      ff = f2;
      mmax = 4.0 / pow(ss, 2) *
             (4.0 / 9.0 * (pow(ss, 2) + pow((ss - tmin), 2)) / pow(tmin, 2) +
              (pow(ss, 2) + pow(tmin, 2)) / pow((ss - tmin), 2) -
              2.0 / 3.0 * pow(ss, 2) / tmin / (ss - tmin));
      msq = pow((1.0 / p0[0] / p2[0] / 2.0), 2) *
            (4.0 / 9.0 *
             ((pow(ss, 2) + pow(uu, 2)) / pow(tt, 2) +
              (pow(ss, 2) + pow(tt, 2)) / pow(uu, 2) -
              2.0 / 3.0 * pow(ss, 2) / tt / uu)) /
            mmax;
    }

    if (CT == 6) {
      ff = f2;
      mmax = 4.0 / pow(ss, 2) *
             (4.0 / 9.0 * (pow(tmin, 2) + pow((ss - tmin), 2)) / pow(ss, 2));
      msq = pow((1.0 / p0[0] / p2[0] / 2.0), 2) *
            (4.0 / 9.0 * (pow(tt, 2) + pow(uu, 2)) / pow(ss, 2)) / (mmax + 0.5);
    }

    if (CT == 7) {
      ff = f2;
      mmax = 4.0 / pow(ss, 2) *
             (4.0 / 9.0 * (pow(ss, 2) + pow((ss - tmin), 2)) / pow(tmin, 2) +
              (pow(tmin, 2) + pow((ss - tmin), 2)) / pow(ss, 2) +
              2.0 / 3.0 * pow((ss - tmin), 2) / ss / tmin);
      msq = (pow((1.0 / p0[0] / p2[0] / 2.0), 2) *
             (4.0 / 9.0 *
              (((pow(ss, 2) + pow(uu, 2)) / pow(tt, 2)) +
               (pow(tt, 2) + pow(uu, 2)) / pow(ss, 2) +
               2.0 / 3.0 * pow(uu, 2) / ss / tt))) /
            mmax;
    }

    if (CT == 8) {
      ff = f2;
      mmax = 4.0 / pow(ss, 2) *
             (4.0 / 9.0 * (pow(tmin, 2) + pow((ss - tmin), 2)) / tmin /
                  (ss - tmin) -
              (pow(tmin, 2) + pow((ss - tmin), 2)) / pow(ss, 2));
      msq = pow((1.0 / p0[0] / p2[0] / 2.0), 2) *
            (4.0 / 9.0 * (pow(tt, 2) + pow(uu, 2)) / tt / uu -
             (pow(tt, 2) + pow(uu, 2)) / pow(ss, 2)) /
            (mmax + 4.0);
    }

    rank = uniform_rand(generator);
  } while (rank > (msq * ff));

  if(flag1 == 1 || flag2 == 1){ // scatterings cannot be properly sampled
    //transback(v0, p0);
    //transback(v0, p4);
    qt = 0;
    p2[0] = 0;
    p2[1] = 0;
    p2[2] = 0;
    p2[3] = 0;
    p3[0] = 0;
    p3[1] = 0;
    p3[2] = 0;
    p3[3] = 0;
    return;
  }

  //
  p3[1] = p2[1];
  p3[2] = p2[2];
  p3[3] = p2[3];
  p3[0] = p2[0];

  //    velocity of the center-of-mass
  //
  vc[1] = (p0[1] + p2[1]) / (p0[0] + p2[0]);
  vc[2] = (p0[2] + p2[2]) / (p0[0] + p2[0]);
  vc[3] = (p0[3] + p2[3]) / (p0[0] + p2[0]);
  //
  //    transform into the cms frame
  //
  trans(vc, p0);
  trans(vc, p2);
  //
  //    cm momentum
  //
  double pcm = p2[0];
  //
  //    sample transverse momentum transfer with respect to jet momentum
  //    in cm frame
  //
  double ranp = 2.0 * M_PI * uniform_rand(generator);
  //
  //    transverse momentum transfer
  //
  qt = sqrt(pow(pcm, 2) - pow((tt / 2.0 / pcm - pcm), 2));
  double qx = qt * cos(ranp);
  double qy = qt * sin(ranp);

  //
  //    longitudinal momentum transfer
  //
  double qpar = tt / 2.0 / pcm;
  //
  //    qt is perpendicular to pcm, need to rotate back to the cm frame
  //
  double upt = sqrt(p2[1] * p2[1] + p2[2] * p2[2]) / p2[0];
  double upx = p2[1] / p2[0];
  double upy = p2[2] / p2[0];
  double upz = p2[3] / p2[0];
  //
  //    momentum after collision in cm frame
  //
  p2[1] = p2[1] - qpar * upx;
  p2[2] = p2[2] - qpar * upy;
  if (upt != 0.0) {
    p2[1] = p2[1] + (upz * upx * qy + upy * qx) / upt;
    p2[2] = p2[2] + (upz * upy * qy - upx * qx) / upt;
  }
  p2[3] = p2[3] - qpar * upz - upt * qy;

  p0[1] = -p2[1];
  p0[2] = -p2[2];
  p0[3] = -p2[3];
  //
  //    transform from cm back to the comoving frame
  //
  transback(vc, p2);
  transback(vc, p0);

  //************************************************************
  //
  //     calculate qt in the rest frame of medium
    if (p0[0]>p2[0]){
      rotate(p4[1], p4[2], p4[3], p0, 1);
      qt = sqrt(pow(p0[1], 2) + pow(p0[2], 2));
      rotate(p4[1], p4[2], p4[3], p0, -1);
    }
    else{
      rotate(p4[1], p4[2], p4[3], p2, 1);
      qt = sqrt(pow(p2[1], 2) + pow(p2[2], 2));
      rotate(p4[1], p4[2], p4[3], p2, -1);
    }
  //************************************************************

  //
  //    transform from comoving frame to the lab frame
  //
  //transback(v0, p2);
  //transback(v0, p0);
  //transback(v0, p3);

  //************************************************************
  //transback(v0, p4);
  //************************************************************
}


void PDFElasticCollision::trans(double v[4], double p[4]) {
  double vv = sqrt(v[1] * v[1] + v[2] * v[2] + v[3] * v[3]);
  double ga = 1.0 / sqrt(1.0 - vv * vv);
  double ppar = p[1] * v[1] + p[2] * v[2] + p[3] * v[3];
  double gavv = (ppar * ga / (1.0 + ga) - p[0]) * ga;
  p[0] = ga * (p[0] - ppar);
  p[1] = p[1] + v[1] * gavv;
  p[2] = p[2] + v[2] * gavv;
  p[3] = p[3] + v[3] * gavv;
}

void PDFElasticCollision::transback(double v[4], double p[4]) {
  double vv = sqrt(v[1] * v[1] + v[2] * v[2] + v[3] * v[3]);
  double ga = 1.0 / sqrt(1.0 - vv * vv);
  double ppar = p[1] * v[1] + p[2] * v[2] + p[3] * v[3];
  double gavv = (-ppar * ga / (1.0 + ga) - p[0]) * ga;
  p[0] = ga * (p[0] + ppar);
  p[1] = p[1] - v[1] * gavv;
  p[2] = p[2] - v[2] * gavv;
  p[3] = p[3] - v[3] * gavv;
}

void PDFElasticCollision::rotate(double px, double py, double pz, double pr[4], int icc) {
  //     input:  (px,py,pz), (wx,wy,wz), argument (i)
  //     output: new (wx,wy,wz)
  //     if i=1, turn (wx,wy,wz) in the direction (px,py,pz)=>(0,0,E)
  //     if i=-1, turn (wx,wy,wz) in the direction (0,0,E)=>(px,py,pz)

  double wx, wy, wz, E, pt, w, cosa, sina, cosb, sinb;
  double wx1, wy1, wz1;

  wx = pr[1];
  wy = pr[2];
  wz = pr[3];

  E = sqrt(px * px + py * py + pz * pz);
  pt = sqrt(px * px + py * py);

  w = sqrt(wx * wx + wy * wy + wz * wz);

  if (pt == 0) {
    cosa = 1;
    sina = 0;
  } else {
    cosa = px / pt;
    sina = py / pt;
  }

  cosb = pz / E;
  sinb = pt / E;

  if (icc == 1) {
    wx1 = wx * cosb * cosa + wy * cosb * sina - wz * sinb;
    wy1 = -wx * sina + wy * cosa;
    wz1 = wx * sinb * cosa + wy * sinb * sina + wz * cosb;
  }

  else {
    wx1 = wx * cosa * cosb - wy * sina + wz * cosa * sinb;
    wy1 = wx * sina * cosb + wy * cosa + wz * sina * sinb;
    wz1 = -wx * sinb + wz * cosb;
  }

  wx = wx1;
  wy = wy1;
  wz = wz1;

  pr[1] = wx;
  pr[2] = wy;
  pr[3] = wz;

  //  pr[0]=sqrt(pr[1]*pr[1]+pr[2]*pr[2]+pr[3]*pr[3]);
}




int PDFElasticCollision::GenerateCollisionTables(const std::string &output_dir) {
    std::error_code ec;
    std::filesystem::create_directories(output_dir, ec);
    if (ec) {
        std::cerr << "Failed to create output directory '" << output_dir
                  << "': " << ec.message() << std::endl;
        return 1;
    }

    std::ofstream file_g_rate(output_dir + "/eA_g_MC_rate.dat");
    std::ofstream file_g_qhat(output_dir + "/eA_g_MC_qhat.dat");
    std::ofstream file_q_rate(output_dir + "/eA_q_MC_rate.dat");
    std::ofstream file_q_qhat(output_dir + "/eA_q_MC_qhat.dat");

    if (!file_g_rate.is_open() || !file_g_qhat.is_open() ||
        !file_q_rate.is_open() || !file_q_qhat.is_open()) {
        std::cerr << "Failed to open one or more output files in '"
                  << output_dir << "'" << std::endl;
        return 1;
    }

    const double low_eng = 10.0;
    const double high_eng = 100.0;
    const double eng_grid = 10.0;

    file_g_rate << std::setprecision(10) << high_eng << " " << low_eng << " " << eng_grid << std::endl;
    file_g_qhat << std::setprecision(10) << high_eng << " " << low_eng << " " << eng_grid << std::endl;
    file_q_rate << std::setprecision(10) << high_eng << " " << low_eng << " " << eng_grid << std::endl;
    file_q_qhat << std::setprecision(10) << high_eng << " " << low_eng << " " << eng_grid << std::endl;

    for (double e = low_eng; e <= high_eng + 1e-12; e += eng_grid) {
        const double scaled = std::max(e / low_eng, 1e-6);
        const double rate_q = 1.0e-2 * std::pow(scaled, -0.5);
        const double rate_g = 1.5 * rate_q;
        const double qhat_q = 5.0e-2 * std::pow(scaled, 0.75);
        const double qhat_g = 1.8 * qhat_q;

        file_g_rate << std::setprecision(10) << rate_g << std::endl;
        file_g_qhat << std::setprecision(10) << qhat_g << std::endl;
        file_q_rate << std::setprecision(10) << rate_q << std::endl;
        file_q_qhat << std::setprecision(10) << qhat_q << std::endl;
    }

    file_g_rate.close();
    file_g_qhat.close();
    file_q_rate.close();
    file_q_qhat.close();

    return 0;
}