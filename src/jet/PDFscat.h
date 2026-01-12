#ifndef PDFSCAT_H
#define PDFSCAT_H

#include <cmath>
 
#include "Math/DistSampler.h"
#include "Math/Factory.h"
#include "Math/IntegratorMultiDim.h"
#include "TF1.h"
#include "Pythia8/Pythia.h"
//#include "Pythia8/Logger.h" //for now

class PDFScat {
	public:
    	PDFScat();
    	~PDFScat();
		void setter(double max_energy0, double low_energy0, double energy_grid0, int index_type0);

		double get_energy(int energy_index);
		double get_energy(double E1, int& energy_index);
		double get_rate(int flv, int energy_index, int process_index);
		void get_sample(int flv, int energy_index, int process_index, double (&V)[4], double EO);
		double get_qhat(int flv, int energy_index, int process_index);
		void GetEIndexRange();
		void initialize_samplers();
		double Integrator_LQ(double E, int pi, int flv);
		double Integrator_LQ_QHat(double E, int pi,int flv);
		double Integrator_HQ(double E, int proc, double eNucleon, double msq);
		double Integrator_HQ_QHat(double E, int pi,double eNucleon, double msq);
		// double Integrator_(double E, double proc);
		// void Sampler_();
		std::vector<double> energy_marker;
		int energy_index_range;

    private:
		double obj_low_energy, obj_hig_energy, obj_grid_energy;
		int obj_index;
		std::vector<std::vector<std::vector<double>>> rates;
		std::vector<std::vector<std::vector<double>>> qhat;
		std::vector<std::vector<std::vector<ROOT::Math::DistSampler*>>> samplers;
		// std::vector<double>energy_marker;
};

double functionToIntegrate_HQ(double *x, double *params);
double functionToIntegrate_LQ(double *x, double *params);
double functionToIntegrate_LQ_QHat(double *x, double *params);
double functionToIntegrate_HQ_QHat(double *x, double *params);

long double PDFSampler(int i, double x, double Q2);

double q1q1b_to_q2q2b(double s, double t, double u);
double q1bq1_to_q2bq2(double s, double t, double u);
double q1q2_to_q1q2(double s, double t, double u);
double q1bq2b_to_q1bq2b(double s, double t, double u);
double q1q1b_to_q1q1b(double s, double t, double u);
double q1bq1_to_q1bq1(double s, double t, double u);
double q1q1_to_q1q1(double s, double t, double u);
double q1bq1b_to_q1bq1b(double s, double t, double u);
double q1q1b_to_gg(double s, double t, double u);
double q1bq1_to_gg(double s, double t, double u);
double q1g_to_q1g(double s, double t, double u);
double q1bg_to_q1bg(double s, double t, double u);
double gq1_to_gq1(double s, double t, double u);
double gg_to_q1q1b(double s, double t, double u);
double gg_to_gg(double s, double t, double u);
double cq_to_cq(double s,double t,double u,double mc_sq);
double cg_to_cg(double s,double t,double u,double mc_sq);

#endif // PDFSCAT_H
