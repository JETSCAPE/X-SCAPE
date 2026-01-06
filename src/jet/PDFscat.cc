#include "PDFscat.h"

//======================================================================== constructor/destructor/setters
PDFScat::PDFScat() {
    Pythia8::Logger logger;
    pythiaPDF = Pythia8::make_shared<Pythia8::LHAGrid1>(2212, "20", "/home/eric/Documents/pythiainst/pythia8315/share/Pythia8/pdfdata", &logger);
}

PDFScat::~PDFScat() {}

void PDFScat::setter(double max_energy0, double low_energy0, double energy_grid0, int index_type0){
    obj_hig_energy = max_energy0;
    obj_low_energy = low_energy0;
    obj_grid_energy = energy_grid0;
    obj_index = index_type0;
}

//======================================================================== getters
double PDFScat::get_energy(double E1, int& energy_index){
        if (obj_index == 1){
                //energy_index = findClosestIndexSorted(energy_marker, E1);
                energy_index = std::round((E1 - obj_low_energy)/obj_grid_energy);
        }
        else if (obj_index == 2){
                energy_index = std::round(std::sqrt(E1));
        }
        else if (obj_index == 3){
                energy_index = std::round(std::pow(E1, 2.0/3.0));
        }
        else{
                energy_index = std::round(10.0 * std::log(E1)) + 10;
        }
        energy_index = std::min(std::max(energy_index, 0), energy_index_range);
        return energy_marker[energy_index];
}

double PDFScat::get_energy(int energy_index){
        return energy_marker[energy_index];
}


double PDFScat::get_rate(int energy_index, int process_index){
    if (energy_index > energy_index_range) { energy_index=energy_index_range; }
    if (energy_index < 0) { energy_index = 0; }
    return rates[process_index][energy_index];
}

void PDFScat::get_sample(int energy_index, int process_index, double (&V)[4]){
    samplers[process_index][energy_index]->Sample(V);
}

long double PDFSampler(int i, double x, double Q2) {
    // Pythia8::PDFPtr pythiaPDF = Pythia8::make_shared<Pythia8::LHAGrid1>(2212, "20", "/home/eric/Documents/pythiainst/pythia8315/share/Pythia8/pdfdata", &logger);
    return pythiaPDF->xf(i, x, Q2);
}

void PDFScat::GetEIndexRange(){
        if (obj_index == 1)
        {
                energy_index_range = std::ceil((obj_hig_energy - obj_low_energy)/ obj_grid_energy);
                for(int i=0; i <= energy_index_range; i++){
                        std::cout<<i<<" pushed energy is "<<i * obj_grid_energy + obj_low_energy<<std::endl;
                        energy_marker.push_back(i * obj_grid_energy + obj_low_energy);
                }
        }
        else if (obj_index == 2)
        {
                energy_index_range = std::ceil(std::sqrt(obj_hig_energy));
                for(int i=0; i <= energy_index_range; i++){
                        energy_marker.push_back(std::pow(i, 2.0));
                }
        }
        else if (obj_index == 3)
        {
                energy_index_range = std::ceil(std::pow(obj_hig_energy, 2.0/3.0));
                for(int i=0; i <= energy_index_range; i++){
                        energy_marker.push_back(std::pow(i, 3.0/2.0));
                }
        }
        else{
	        energy_index_range = 10.0 * std::ceil(std::log(obj_hig_energy)) + 10.0;
                //std::cout<<"energy_index_range is "<<energy_index_range<<std::endl;
                for(int i=0; i <= energy_index_range; i++){
                        double e_temp = std::exp((i - 10.0) / 10.0);
                        std::cout<<i<<" pushed energy is "<<e_temp<<std::endl;
                        energy_marker.push_back(e_temp);
                }
        }
}

//======================================================================== sampler initialization

void PDFScat::initialize_samplers(double eNucleon){
    const double chm_massSq = 1.6129;
    const double btm_massSq = 17.4724;
    // const double proc_[13] = {1,2,3,4,5,6,7,8,9,1,2,1,2};
    GetEIndexRange();

    //energy_index_range = int((obj_hig_energy - obj_low_energy) / obj_grid_energy);


    double local_energy, temp_rate;
    TF1 *f1;
    ROOT::Math::DistSampler *sampler;
    double xmin[4], xmax[4], par0[];
    bool ret;

    for (int pi=0; pi<13; pi++){
        std::vector<double> rates_row;
        std::vector<ROOT::Math::DistSampler*> samplers_row;                     
        for (int ei=0; ei<=energy_index_range; ei++) {    
            local_energy = energy_marker[ei];

            if (pi<9) { //massless partons
                temp_rate = Integrator_LQ(local_energy, pi, eNucleon);
                f1 = new TF1("myfunc", functionToIntegrate_LQ, 0, 1, 3); //TODO: WHAT ARE THESE ARGS
                sampler = ROOT::Math::Factory::CreateDistSampler("Foam");

                //xmin = {0.0, 0.0, 0.0, 0.0}; //TODO: WHAT'S THE ORDER FOR THESE
                xmin = {M_PI - 1e-4, 0.0, 0.0, 0.0};
                //xmax = {M_PI, M_PI, 2 * M_PI, local_energy + 1}; //TODO: UPPER LIMIT FOR ENERGY?
                xmax = {M_PI, M_PI, 2.0 * M_PI, local_energy * 1.2};
                par0 = {local_energy, (double)pi, eNucleon}; //TODO: CHANGE THE PARAMS OBV
            }
            else if (pi>=9 && pi<11) { //charm quark
                temp_rate = Integrator_HQ(local_energy, pi, eNucleon, chm_massSq);
                f1 = new TF1("myfunc", functionToIntegrate_HQ, 0, 1, 4); //TODO: WHAT ARE THESE ARGS
                sampler = ROOT::Math::Factory::CreateDistSampler("Foam");

                //xmin = {0.0, 0.0, 0.0, 0.0}; //TODO: WHAT'S THE ORDER FOR THESE
                xmin = {0.0, M_PI - 1e-4, 0.0, 0.0};
                //xmax = {local_temp*15, M_PI, M_PI, 2.0*M_PI}; //TODO: UPPER LIMIT FOR ENERGY?
                xmax = {local_energy, M_PI, M_PI, 2.0*M_PI};
                par0 = {local_energy, (double)pi, eNucleon, chm_massSq}; //TODO: CHANGE THE PARAMS OBV                                    
            }
            else { //bottom quark
                temp_rate = Integrator_HQ(local_energy, pi, eNucleon, btm_massSq);
                f1 = new TF1("myfunc", functionToIntegrate_HQ, 0, 1, 4); //TODO: WHAT ARE THESE ARGS
                sampler = ROOT::Math::Factory::CreateDistSampler("Foam");
                
                xmin = {0.0, M_PI - 1e-4, 0.0, 0.0}; //TODO: WHAT'S THE ORDER FOR THESE
                xmax = {local_energy, M_PI, M_PI, 2.0*M_PI}; //TODO: UPPER LIMIT FOR ENERGY?
                par0 = {local_energy, (double)pi, eNucleon, btm_massSq}; //TODO: CHANGE THE PARAMS OBV                               
            }

            f1->SetParameters(par0);
            sampler->SetFunction(*f1, 4); //4 is num arguments
            sampler->SetRange(xmin, xmax);
            
            ret = sampler->Init();

            samplers_row.push_back(sampler);
            rates_row.push_back(temp_rate);
        }
        rates.push_back(rates_row);
        samplers.push_back(samplers_row);
    }

    delete f1; //?
}


//======================================================================== vegas integrators

double PDFScat::Integrator_LQ(double E, int pi, double eNucleon){
    TF1 *f1 = new TF1("myfunc", functionToIntegrate_LQ, 0, 1, 3); //args are min, max, #args -- min and max are overwritten anyway

    double xmin[] = {M_PI - 1e-4, 0.0, 0.0, 0.0}; //TODO: SAME AS
    //double xmax[] = {M_PI, M_PI, 2*M_PI, E*1.2};
    double xmax[] = {M_PI, M_PI, 2.0 * M_PI, E*1.2}; //p2 is close to colliniear in our refence frame
    //double par0[] = {E, (double)pi, eNucleon, pzNucleon}; 
    double par0[] = {E, (double)pi, eNucleon}; 

    f1->SetParameters(par0);
    ROOT::Math::IntegratorMultiDim ig(ROOT::Math::IntegrationMultiDim::kVEGAS, 1.E-10, 1.E-10, 50000); //abs error, rel error, #runs
    ig.SetFunction(*f1, 4); //arg is num dims

    double int_result = ig.Integral(xmin,xmax);
    delete f1;

    return int_result;
}

double PDFScat::Integrator_HQ(double E, int pi, double eNucleon, double msq){
    TF1 *f1 = new TF1("myfunc", functionToIntegrate_HQ, 0, 1, 4); //extra arg because of msq

    double xmin[] = {0.0, M_PI - 1e-4, 0.0, 0.0}; //TODO: SAME AS
    //double xmax[] = {15*t, M_PI, M_PI, 2.0*M_PI};
    double xmax[] = {E, M_PI, M_PI, 2.0*M_PI}; //p2 is close to colliniear in our refence frame
    double par0[] = {E, (double)pi, eNucleon, msq}; 

    f1->SetParameters(par0);
    ROOT::Math::IntegratorMultiDim ig(ROOT::Math::IntegrationMultiDim::kVEGAS, 1.E-10, 1.E-10, 50000);
    ig.SetFunction(*f1,4);
  
    double int_result = ig.Integral(xmin,xmax);
    delete f1;

    return int_result;
}


//======================================================================== cross sections

//TODO: MAKE SURE THESE ARE ALL CORRECT

inline double q1q1b_to_q2q2b(double s, double t, double u) {
    return (4.0/9)*(pow(t,2)+pow(u,2))/pow(s,2);
}

inline double q1bq1_to_q2bq2(double s, double t, double u) {
    return (4.0/9)*(pow(t,2)+pow(u,2))/pow(s,2);
}

inline double q1q2_to_q1q2(double s, double t, double u) {
    return (4.0/9)*(pow(s,2)+pow(u,2))/pow(t,2);
}

inline double q1bq2b_to_q1bq2b(double s, double t, double u) {
    return (4.0/9)*(pow(s,2)+pow(u,2))/pow(t,2);
}

inline double q1q1b_to_q1q1b(double s, double t, double u) {
    return (4.0/9)*((pow(s,2)+pow(u,2))/pow(t,2)+(pow(t,2)+pow(u,2))/pow(s,2)-2*pow(u,2)/(3*s*t));
}

inline double q1bq1_to_q1bq1(double s, double t, double u) {
    return (4.0/9)*((pow(s,2)+pow(u,2))/pow(t,2)+(pow(t,2)+pow(u,2))/pow(s,2)-2*pow(u,2)/(3*s*t));
}

inline double q1q1_to_q1q1(double s, double t, double u) {
    return (4.0/9)*((pow(u,2)+pow(s,2))/pow(t,2)+(pow(t,2)+pow(s,2))/pow(u,2)-2*pow(s,2)/(3*u*t));
}

inline double q1bq1b_to_q1bq1b(double s, double t, double u) {
    return (4.0/9)*((pow(u,2)+pow(s,2))/pow(t,2)+(pow(t,2)+pow(s,2))/pow(u,2)-2*pow(s,2)/(3*u*t));
}

inline double q1q1b_to_gg(double s, double t, double u) {
    return (32.0/27)*(u/t + t/u -(9.0/4)*(pow(t,2)+pow(u,2))/pow(s,2));
}

inline double q1bq1_to_gg(double s, double t, double u) {
    return (32.0/27)*(u/t + t/u -(9.0/4)*(pow(t,2)+pow(u,2))/pow(s,2));
}

inline double q1g_to_q1g(double s, double t, double u) {
    return (4.0/9)*(((-1*u)/s)+((-1*s)/u)+(9.0/4)*(pow(s,2)+pow(u,2))/pow(t,2));
}

inline double q1bg_to_q1bg(double s, double t, double u) {
    return (4.0/9)*(((-1*u)/s)+((-1*s)/u)+(9.0/4)*(pow(s,2)+pow(u,2))/pow(t,2));
}

inline double gq1_to_gq1(double s, double t, double u) {
    return (4.0/9)*(((-1*u)/s)+((-1*s)/u)+(9.0/4)*(pow(s,2)+pow(u,2))/pow(t,2));
}

inline double gg_to_q1q1b(double s, double t, double u) {
    return (1/6.0)*(u/t+t/u-(9.0/4)*(pow(t,2)+pow(u,2))/pow(s,2));
}

inline double gg_to_gg(double s, double t, double u) {
    return (9.0/2)*(3-t*u/pow(s,2)-s*u/pow(t,2)-s*t/pow(u,2));
}

inline double cq_to_cq(double s,double t,double u,double mc_sq){
        return ((4.0/9.0)*(pow((mc_sq-u),2)+pow(s-mc_sq,2)+2*mc_sq*t)/pow(t,2));
}

inline double cg_to_cg(double s,double t,double u,double mc_sq){
        return ((2.0*(s-mc_sq)*(mc_sq-u))/pow(t,2)+
                (4.0/9.0)*(((s-mc_sq)*(mc_sq-u)+2.0*mc_sq*(s+mc_sq))/(pow(s-mc_sq,2))+((s-mc_sq)*(mc_sq-u)+2.0*mc_sq*(u+mc_sq))/(pow(mc_sq-u,2))+
        (mc_sq*(4.0*mc_sq-t))/(4.0*(s-mc_sq)*(mc_sq-u)))
                +(1.0)*(((s-mc_sq)*(mc_sq-u)+mc_sq*(s-u))/(t*(s-mc_sq))-
        ((s-mc_sq)*(mc_sq-u)-mc_sq*(s-u))/(t*(mc_sq-u))));
}


//======================================================================== integrands


double functionToIntegrate_LQ(double *x, double *params) {
    double E1     = params[0];
    int proc = (int)params[1];
    double eNucleon = params[2]; //isn't our coordinate system such that nucleon is moving along z axis?
    //double pzNucleon = params[3];

    //return 0.0;
    //theta2, theta3 , phi23, E3
    double c23 = cos(x[0])*cos(x[1])+sin(x[0])*sin(x[1])*cos(x[2]); //cos(theta_23)
    double E2 = E1*x[3]*(1-cos(x[1]))/(E1*(1-cos(x[0]))-x[3]*(1-c23)); //nuclear parton energy

    double s = 2*E1*E2*(1-cos(x[0]));
    double t = -2*E1*x[3]*(1-cos(x[1]));
    double u = -s-t;

    //find alphas
    // double EScale2 = std::max(std::pow(std::abs(t),2.0),std::pow(M_PI*Temp,2.0));
    double EScale2 = std::pow(std::abs(t),2.0);
    double Nf = 3.0;
    double QCDScale2 = std::pow(0.2,2.0); //MS bas scheme 0.33; 0.2 LBT value
    double alphaS = 0.3;//4.0*M_PI/(11.0-2.0*Nf/3.0)*1/std::log(EScale2/QCDScale2);
    double lambdaQCD = 0.2; //GeV
    // double debye = alphaS*pow(Temp,2.0)*(4.0*M_PI)*1.5;
     if(s <= 2*lambdaQCD*lambdaQCD || t <= -s+lambdaQCD*lambdaQCD || t >= -lambdaQCD*lambdaQCD){
             return 0;
     }

    double matrix_element, stat;
    int g_b; //generacy factor
    switch (proc) {
        case 0:
            matrix_element = 0.5*q1q1b_to_q2q2b(s,t,u); // 0.5*q1q1b_to_q2q2b(s,t,u);     //problem
            stat = 0.;
            g_b = 6;
            break;
        
        case 1:
            matrix_element = 0.5*q1q1b_to_q1q1b(s,t,u); //0.5*q1q1b_to_q1q1b(s,t,u);     //problem
            stat = 0.;
            g_b = 6;
            break;

        case 2:
            matrix_element = 0.5*q1q1_to_q1q1(s,t,u);       //correct
            stat = 0.;
            g_b = 6;
            break;
        
        case 3:
            matrix_element = 0.5*q1q1b_to_gg(s,t,u); //correct
            stat = 0.;
            g_b = 6;
            break;
            
        case 4:
            matrix_element = q1g_to_q1g(s,t,u);     //correct
            stat = 1.;
            g_b = 16;
            break;
        
        case 5:
            matrix_element = 4*q1q2_to_q1q2(s,t,u); //4 q and qb contribution of both considered    //correct
            stat = 0.;
            g_b = 6;
            break;
            
        case 6:
            matrix_element = gg_to_q1q1b(s,t,u);            //correct
            stat = 1.;
            g_b = 16;
            break;
        
        case 7:
            matrix_element = 0.5*gg_to_gg(s,t,u);           //correct
            stat = 1.;
            g_b = 16;
            break;
        
        case 8:
            matrix_element = 6*gq1_to_gq1(s,t,u);        //both q and qb contribution //correct
            stat = 0.;
            g_b = 6;
            break;
       }

    // double thermal_dist = 0;
    // thermal_dist = ThermalDistribution(E2,Temp,Stat,b_val);

    //x = (p3-p1)^2 / 2 Pnuc dot (p3-p1)
    //p1 and p3 onshell so (p3-p1)^2 = 2p1 dot p3
    //p1 and Pnuc are approximately along the z axis
    double p3z = -x[3]*c23;
    double x = t / (eNucleon*(x[3] - p3z))/2.0);
    //double x = 0.;
    double pdf = PDFSampler(1,x,t);

    double ans = std::pow(alphaS,2.0)*pow(4*M_PI,2)*g_b/(pow(2*M_PI,4)*16.0*E1)*pdf*matrix_element*pow(E2,2)*2*x[3]*sin(x[0])*sin(x[1])/abs(t);
    return ans;
}

double functionToIntegrate_HQ(double *x, double *params) {
    double E1     = params[0];
    int proc = (int)params[1];
    double eNucleon = params[2];
    //double pzNucleon = params[3];
    double msq    = params[3];

    //return 0.0;
    //x[0]=E2; x[1]=theta2; x[2]=theta4; x[3]=phi4
    //on-shell
    double c24 = sin(x[1])*sin(x[2])*cos(x[3])+cos(x[1])*cos(x[2]);
    double p1  = sqrt(E1*E1-msq);
    double E4  = (E1*x[0]- p1*x[0]*cos(x[1]))/(E1-p1*cos(x[2])+x[0]-x[0]*c24);
        
    double s = msq+ 2*(E1*x[0]-p1*x[0]*cos(x[1]));
    double u = msq-2*(E1*E4-p1*E4*cos(x[2]));
    double t = 2*msq-s-u;

    //find alphas
    double EScale2 = std::pow(std::abs(t),2.0);
    double Nf = 3.0;
    double QCDScale2 = std::pow(0.2,2.0); //MS bas scheme 0.33; 0.2 LBT value
    double alphaS = 0.3;//4.0*M_PI/(11.0-2.0*Nf/3.0)*1/std::log(EScale2/QCDScale2);
`   double lambdaQCD = 0.2; //GeV
    // double debye = 0.3*pow(Temp,2)*(4*M_PI)*(4.5/3.0);
    // //std::cout<<"Values fom scattering.cpp s "<<s<<" t "<<t<<" u "<<u<<std::endl;
     if(s <= 2*lambdaQCD*lambdaQCD || t <= -s+lambdaQCD*lambdaQCD || t >= -lambdaQCD*lambdaQCD){
    //         //std::cout<<"RETURN 0!";
             return 0;
    }

    double matrix_element, stat;
    int g_b;
    switch (proc) {
        case 9: case 11: //first process for charm, bottom
            matrix_element = 6*cq_to_cq(s,t,u,msq);
            stat = 0.;
            g_b = 6;
            break;

        case 10: case 12: //second process for charm, bottom
            matrix_element = cg_to_cg(s,t,u,msq);
            stat = 1.;
            g_b = 16;
            break;  
    }

    //TODO: what is extra_stat etc
    long double thermal_dist = ThermalDistribution(x[0],Temp,Stat,b_val);
    long double extra_stat = ThermalDistribution(E4,Temp,Stat,b_val);
    if (stat==0){ extra_stat=extra_stat*(-1.0); }

    double p3z = -x[3]*c23;
    double x = t / (eNucleon*(x[3] - p3z))/2.0);
    //double x = 0.;
    double pdf = PDFSampler(1,x,t);

    double ans = pow(alphaS,2)*pow(4*M_PI,2)*g_b/(pow(2*M_PI,4)*16.0*E1)*pdf*matrix_element*(1+extra_stat)*(x[0]*E4*sin(x[1])*sin(x[2]))/(E1-p1*cos(x[2])+x[0]-x[0]*c24);
    return ans;
}