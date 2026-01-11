#include "PDFElasticCollision.h"

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

void PDFElasticCollision::setter(double max_energy0, double low_energy0, double energy_grid0,int index_type0, double nu, double Q2) {
    obj_hig_energy = max_energy0;
    obj_low_energy = low_energy0;
    obj_grid_energy = energy_grid0;
    obj_index = index_type0;
    nu_ = nu;
    Q2_ = Q2;
    //iZ_Vector.SetXYZ(0,0,1); //won't be using root to do the rotations
    scattering_obj.setter(max_energy0, low_energy0, energy_grid0, index_type0);
	//scattering_obj.setter(0.2,0.2,0.1,102,0.5,0.5);
	//scattering_obj.setter(0.2,0.2,0.1,102,0.5,1);
    scattering_obj.initialize_samplers(nu, Q2); //the param here is energy of nucleon getting scattered off 
}
/*The MATTER module had pre-defination  of calculating the probability of scattring.*/

bool PDFElasticCollision::elastic_kinematics(bool ProbailisticScattering, double deltaT, int &pid0, int &pid2, int &pid3, double (&pc0)[4], double (&pc2)[4], double (&pc3)[4], double &qt) {
    //pc2 is hole
    //pc3 is energetic parton after scattering
    //pc4 is recoil
    //double r0,r1,r2,r3,r4,r5;
    double EO,E1,E2,E3,E4,c23,c24,p1,p3,THETA2,THETA3,THETA4,PHI3,PHI4,PHI2,PHI23;
    //EO is the original energy
    double mc_sq;
    double TotalRate = 0;
    double beta =  nu_ / sqrt(nu_ * nu_ + Q2_);
    double gamma = sqrt(nu_ * nu_ + Q2_) / sqrt(Q2_);

    int parent_pid, pid_index;
    int daughter1_pid = -1;
    int daughter2_pid = -1;
    int hole_pid = -1;
    int energy_index;//, temp_index;
    int parton_type = -1;
    int proc; //heavy(1) or light(0)

    //TVector3 P1;
    //TVector3 P2;
    ////TVector3 P3;
    //TVector3 P4;
    //TRotation r;

    EO = pc0[0];
    EO = gamma * (EO - beta * EO);//breit frame
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
        TotalRate = r0g + r1g + r2u + r2ubar + r2d + r2dbar + r2s + r2sbar;
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
    
    else if (parent_pid ==1) {
		double r0 = 0; //PerformLinearInterpolation(EO,0,1);//s-channel not considered
		double r1bar = PerformLinearInterpolation(EO,1,-parent_pid,1);
		double r2 = PerformLinearInterpolation(EO,2,parent_pid,1);
		double r3bar = PerformLinearInterpolation(EO,3,-parent_pid,1);
		double r4 = PerformLinearInterpolation(EO, 4, 21, 1);
		double r5a = PerformLinearInterpolation(EO,5,(parent_pid < 0 ? -1 : 1) * (((abs(parent_pid)) % 3) + 1),1);
        double r5b = PerformLinearInterpolation(EO,5,(parent_pid < 0 ? -1 : 1) * (((abs(parent_pid)) % 3) + 2),1);
        TotalRate = r0 + r1bar + r2 + r3bar + r4 + r5a + r5b;
        if (ProbailisticScattering){  
            if (exp(-TotalRate * deltaT) > uniform_rand(generator)){qt = 0; return false;}
        }  
        parton_type = 0;
    	std::discrete_distribution<int> distribution0{r0,r1bar,r2,r3bar,r4,r5a,r5b};
        switch (distribution0(generator)) {
                case 0: //q1 q1bar -> q2 q2bar
                    /*
                    scattering_obj.get_sample(energy_index,0,V, EO);
                    hole_pid      = -parent_pid;
                    do { //keep sampling until q2 != q1
                        pid_index = floor(uniform_rand(generator)*6);
                        daughter1_pid = pid_list[pid_index];
                    } while (daughter1_pid == parent_pid);                
                    daughter2_pid = -daughter1_pid;
                    */
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
                    /*
                    do { //keep sampling until q2 != q1 or q1bar
                        pid_index = floor(uniform_rand(generator)*6);
                        daughter2_pid = pid_list[pid_index];
                    } while (daughter2_pid == parent_pid || daughter2_pid == -parent_pid);
                     */
                    hole_pid=daughter2_pid;     
                    break;
                case 6: //q1 q2 -> q1 q2
                    daughter2_pid = (parent_pid < 0 ? -1 : 1) * (((abs(parent_pid)) % 3) + 2);
                    scattering_obj.get_sample(daughter2_pid,energy_index,5,V, EO);
                    daughter1_pid = parent_pid;
                    hole_pid=daughter2_pid;     
                    break;
                default:
                    //never gets here
                    break;
        }
    }
     /*Heavy scattering part will be added here*/
    else if (abs(parent_pid)==4) {
        /*
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
        }*/
         return false;
    }
        
    else {
        /*
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
        }*/
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
		
        PHI2 = 2 * M_PI * uniform_rand(generator);
		PHI3 = PHI2 - PHI23;

		//THETA4= acos((E1+E2*cos(THETA2)-E3*cos(THETA3))/E4);
        THETA4= acos((EO + E2 * cos(THETA2) - E3 * cos(THETA3)) / E4);

		double value = (E2 * sin(THETA2) * cos(PHI2) - E3 * sin(THETA3) * cos(PHI3)) / (E4 * sin(THETA4));
		value = std::max(-1.0, std::min(1.0, value)); // Clamping the value
		PHI4 = acos(value);

        //std::cout<<"LQ E1 "<<E1<<" E2 "<<E2<<" E3 "<<E3<<" E4 "<<E4<<std::endl;
        //pc2[0] = E2;
        pc2[1] = E2 * sin(THETA2) * cos(PHI2);
        pc2[2] = E2 * sin(THETA2) * sin(PHI2);
        //pc2[3] = E2 * cos(THETA2);
        pc2[0] = gamma * (E2 + beta * E2 * cos(THETA2));
        pc2[3] = gamma * (E2 * cos(THETA2) + beta * E2);
        pid2 = hole_pid;

        //pc0[0] = E3;
        pc0[1] = E3 * sin(THETA3) * cos(PHI3);
        pc0[2] = E3 * sin(THETA3) * sin(PHI3);
        //pc0[3] = E3 * cos(THETA3);
        pc0[0] = gamma * (E3 + beta * E3 * cos(THETA3));
        pc0[3] = gamma * (E3 * cos(THETA3) + beta * E3);
        pid0 = daughter1_pid;

        //pc3[0] = E4;
        pc3[1] = E4 * sin(THETA4) * cos(PHI4);
        pc3[2] = E4 * sin(THETA4) * sin(PHI4);
        //pc3[3] = E4 * cos(THETA4);   
        pc3[0] = gamma * (E4 + beta * E4 * cos(THETA4));
        pc3[3] = gamma * (E4 * cos(THETA4) + beta * E4);
        pid3 = daughter2_pid;
        
        //Will perform rotation outside this class;

        //P1.SetXYZ(pc0[1],pc0[2],pc0[3]);
        //P2.SetXYZ(E2*sin(THETA2)*cos(PHI2),E2*sin(THETA2)*sin(PHI2),E2*cos(THETA2));
		//P3.SetXYZ(E3*sin(THETA3)*cos(PHI3),E3*sin(THETA3)*sin(PHI3),E3*cos(THETA3));
		//P4.SetXYZ(E4*sin(THETA4)*cos(PHI4),E4*sin(THETA4)*sin(PHI4),E4*cos(THETA4));
		/*
		double s0=P1.Angle(iZ_Vector);
		TVector3 s1;
		s1=P1.Cross(iZ_Vector);
        
		r.Rotate(s0,s1);
        P2=r*P2;
        P3=r*P3;
        P4=r*P4;
		*/
        //std::cout<<"LQ  after one roatation E1 "<<E1<<" E2 "<<E2<<" E3 "<<E3<<" E4 "<<E4<<std::endl;
    }
    else if (parton_type==1) { //heavy
        /* Will implement later
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
        */
    }
    else { return false; }
    /*
    //TODO: LEARN ROTATIONS ETC
    double s0 = P1.Angle(iZ_Vector);
    TVector3 s1;
    //s1=P1.Cross(iZ_Vector);
    s1=iZ_Vector.Cross(P1);
    TRotation w0;
    w0.Rotate(s0,s1);
    
    P2=w0*P2;
    P3=w0*P3;
    P4=w0*P4;
 	Rectify1(P2);
 	Rectify1(P3);
 	Rectify1(P4);

    pc3[0]=E2;
    pc3[1]=P2.x();
    pc3[2]=P2.y();
    pc3[3]=P2.z();
    pid3=hole_pid;

    pc0[0]=E3;
    pc0[1]=P3.x();
    pc0[2]=P3.y();
    pc0[3]=P3.z();
    pid0=daughter1_pid;

    pc2[0]=E4;
    pc2[1]=P4.x();
    pc2[2]=P4.y();
    pc2[3]=P4.z();   
    pid2=daughter2_pid; 
    */
    return true;
}

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
        qhat = PerformLinearInterpolation(energy_index, 6, 21, 2)
            + PerformLinearInterpolation(energy_index, 7, 21, 2)
            + PerformLinearInterpolation(energy_index, 8, 1, 2)
            + PerformLinearInterpolation(energy_index, 8, -1, 2)
            + PerformLinearInterpolation(energy_index, 8, 2, 2)
            + PerformLinearInterpolation(energy_index, 8, -2, 2)
            + PerformLinearInterpolation(energy_index, 8, 3, 2)
            + PerformLinearInterpolation(energy_index, 8, -3, 2);
    }
    else if (abs(pid) < 3){
        qhat = 0
            + PerformLinearInterpolation(energy_index, 1, -pid, 2)
            + PerformLinearInterpolation(energy_index, 2, pid, 2)
            + PerformLinearInterpolation(energy_index, 3, -pid, 2)
            + PerformLinearInterpolation(energy_index, 4, 21, 2)
            + PerformLinearInterpolation(energy_index, 5, (pid < 0 ? -1 : 1) * (((abs(pid)) % 3) + 1), 2)
            + PerformLinearInterpolation(energy_index, 5, (pid < 0 ? -1 : 1) * (((abs(pid)) % 3) + 2), 2); 
    }
    else{
        qhat = 0;
    }
    return qhat;
}

double PDFElasticCollision::PerformLinearInterpolation(double EOriginal, int process_id,int flv, int type){
    int EIndex, EIndex_; //index where sampler is trained
    double E, E_;
    double ERoundedOff = scattering_obj.get_energy(EOriginal, EIndex);
    //std::cout<<std::setprecision(4)<<"ERoundedOff "<<ERoundedOff<<"  EOriginal "<<EOriginal<<std::endl;
    if (ERoundedOff < EOriginal){
        EIndex_ = std::min(EIndex + 1, scattering_obj.energy_index_range);
    }
    
    else{
        EIndex_ = EIndex;
        EIndex = std::max(0, EIndex_ - 1);
    }
    E_ = scattering_obj.get_energy(EIndex_);
    E = scattering_obj.get_energy(EIndex);
    //std::cout<<std::setprecision(4)<<"E_ "<<E_<<" E "<<E<<std::endl;
    //std::cout<<"EOriginal "<<EOriginal<<" EIndex "<<EIndex<<" EIndex_ "<<EIndex_<<" E_ "<<E_<<" E "<<E<<std::endl;
    /*
    int TIndex, TIndex_;
    double T, T_;
    double TRoundedOff = scattering_obj.get_temp(TOriginal, TIndex);
    //std::cout<<std::setprecision(4)<<"TRoundedOff "<<TRoundedOff<<"  TOriginal "<<TOriginal<<std::endl;
    if (TRoundedOff < TOriginal){
        TIndex_ = std::min(TIndex + 1, scattering_obj.temp_index_range);
    }
    else{
        TIndex_ = TIndex;
        TIndex = std::max(0, TIndex_ - 1);
    }
    T_ = scattering_obj.get_temp(TIndex_);
    T = scattering_obj.get_temp(TIndex);
    */
    //std::cout<<std::setprecision(4)<<"T_ "<<T_<<" T "<<T<<std::endl;
    //std::cout<<"TOriginal "<<TOriginal<<" TIndex "<<TIndex<<" TIndex_ "<<TIndex_<<" T_ "<<T_<<" T "<<T<<std::endl;

    double ValE1, ValE2;
    double ValFinal;
    if (type == 1){
        //Rate interpolation
        ValE1 = scattering_obj.get_rate(flv, EIndex, process_id);
        ValE2 = scattering_obj.get_rate(flv, EIndex_, process_id);
    }
    else{
        //Qhat interpolation
        ValE1 = scattering_obj.get_qhat(flv, EIndex, process_id);
        ValE2 = scattering_obj.get_qhat(flv, EIndex_, process_id);
    }
    /*
    if (T_ - T < 1e-4){
        if (E_ - E < 1e-4){
            //Takes care when out of range
            return ValT1E1;
        }
        //Takes care when out of range
        return (ValT1E2 - ValT1E1) * (EOriginal - E) / (E_ - E) + ValT1E1;
    }
*/
    //std::cout<<std::setprecision(4)<<"EOriginal "<<EOriginal<<" TOriginal "<<TOriginal<<" ValT1E1 "<<ValT1E1<<" ValT1E2 "<<ValT1E2<<" ValT2E1 "<<ValT2E1<<" ValT2E2 "<<ValT2E2<<std::endl;

    //ValEgridLow = (ValT2E1 - ValT1E1) * (TOriginal - T) / (T_ - T) + ValT1E1;

    //if (E_ - E < 1e-4){
        //Takes care when out of range
    //    return ValEgridLow;
    //}

    //ValEgridHigh = (ValT2E2 - ValT1E2) * (TOriginal - T) / (T_ - T) + ValT1E2;
    //ValFinal = (ValEgridHigh - ValEgridLow) * (EOriginal - E) / (E_ - E) + ValEgridLow;
    ValFinal = (ValE2 - ValE1) * (EOriginal - E) / (E_ - E) + ValE1;
    //std::cout<<std::setprecision(4)<<"ValEgridLow "<<ValEgridLow<<" ValEgridHigh "<<ValEgridHigh<<" ValFinal "<<ValFinal<<std::endl;
    return ValFinal;
}
/*
int main(){
    PDFElasticCollision pdfElasticCollision;
    //pdfElasticCollision.setter(10.0, 1.0, 0.5, 1, 1000.0);
    //std::cout<<PDFSampler(21, 0.1, 10.0)<<std::endl;
    return 0;
}*/