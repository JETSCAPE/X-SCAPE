#include <iostream>    //....C++ headers
#include <fstream>
#include <sstream>
#include <iomanip>
#include <vector>
#include <string>
#include <cstdio>    //....C headers
#include <cstdlib>
#include <cmath>
#include <ctime>
//#include "Pythia8/Pythia.h"    //....PYTHIA8 headers
#include "fastjet/ClusterSequence.hh"    //....FASTJET headers

/*
    The script to calculate the jet pt-spectra.
    Input file is the X-SCAPE output of hard hadrons.
*/

using namespace std;
using namespace fastjet;

int main(int argc, char* argv[]) {

    const int Nevent = 400;
    const int length_pt = 25;
    // selection for final particles which are used to reconstruct jet
    double absetamax = 0.9;
    double particle_ptmin = 0.15; // ALICE cut, arXiv: 1909.09718
    // parameter setting
    const double R_plus_eta = 0.7; //  ALICE cut, arXiv: 1909.09718
    const double jet_ptmin = 5.0;
    double R, jet_absetamax;
    vector<fastjet::PseudoJet> input_particles;
    double dndptdyjet[length_pt][6] = {0.0};
    double pT_array_temp[] = {5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 12.0, 14., 16., 18, 20., 22., 24, 26, 28, 30., 34., 36, 40., 50., 60.,
                         80., 100., 120, 160., 200.};
    vector<double> pT_array;
    double pT_diff[] = {0.0};
    for (int ipt = 0; ipt < length_pt; ipt ++) {
        pT_diff[ipt] = pT_array_temp[ipt+1] - pT_array_temp[ipt];
        pT_array.push_back(pT_array_temp[ipt]);
    }
    pT_array.push_back(pT_array_temp[length_pt]);
    char inputfile[128];
    sprintf(inputfile, "test_out_final_state_hadrons.dat");
    FILE* infile;
    infile = fopen(inputfile,"r");
    char* stemp1;
    char** stemp2;
    int total_number_of_particles, pid, event_id, int_temp, status;
    double px, py, pz, energy, mass, dummpx, dummpy, dummpz, dummpt, weight;
    double sigmaGen, sigmaErr, pThat, sigmaGen_ev;
    fscanf(infile,"%s %s %s %s %s %s %s %s %s %s %s\n",&stemp1, &stemp2,
           &stemp2, &stemp2, &stemp2, &stemp2, &stemp2, &stemp2, &stemp2, &stemp2,
           &stemp2);
    double weight_sum = 0.0;
    int even_loop_flag = 1;
    int count_event_number = 0;
    for (int iev = 0; iev < Nevent; iev ++) {
        if(feof(infile)) {
            even_loop_flag = 0;
            cout << " End the event loop ~~~ " << endl;
            break;
        }
        fscanf(infile,"%s %s %d %s %lf %s %d %s %d %lf\n",&stemp1, &stemp2,
               &event_id, &stemp2, &weight, &stemp2, &int_temp, &stemp2, 
               &total_number_of_particles, &sigmaGen_ev);
        //weight_array.push_back(weight);
        input_particles.clear();
        for (auto i=0; i<total_number_of_particles; i++) {
            if(feof(infile)) {
                even_loop_flag = 0;
                cout << " End the event loop, and drop last event ~~~ " << endl;
                break;
            }
            fscanf(infile,"%d %d %d %lf %lf %lf %lf\n",&int_temp, &pid,
                   &status, &energy, &px, &py, &pz);
            if (status == 11) continue; // don't count hydro hadrons.
            fastjet::PseudoJet particle = PseudoJet(px, py, pz, energy);
            particle.set_user_index(pid);
            input_particles.push_back(particle);
        }
        if(even_loop_flag == 0) {
            cout << " End the event loop and drop last event ~~~ " << endl;
            break;
        }
        count_event_number++;
        weight_sum = weight_sum + weight;

        fastjet::Selector particle_selector = fastjet::SelectorAbsEtaMax(absetamax) && fastjet::SelectorPtMin( particle_ptmin );
        // a jet algorithm with a given radius parameter
        for (int ir = 1; ir < 7; ir ++) {
            R = ir * 0.1;
            jet_absetamax = R_plus_eta - R;
            fastjet::JetDefinition jet_def(fastjet::antikt_algorithm, R);
            // select jet
            fastjet::Selector jet_selector = fastjet::SelectorAbsEtaMax( jet_absetamax ) && fastjet::SelectorPtMin( jet_ptmin );
            input_particles = particle_selector(input_particles);
            fastjet::ClusterSequence clust_seq(input_particles, jet_def);
            // get the resulting jets ordered in pt
            vector<fastjet::PseudoJet> inclusive_jets = sorted_by_pt( jet_selector(clust_seq.inclusive_jets()) );

            // calculate the jet pT-spectra
            for (unsigned int i = 0; i < inclusive_jets.size(); i++) {
                if (abs(inclusive_jets[i].pseudorapidity()) < jet_absetamax) {
                    for (int kk = 0; kk < length_pt; kk++) {
                        if ((inclusive_jets[i].pt() > pT_array[kk]) && (inclusive_jets[i].pt() <= pT_array[kk+1])) {
                            dndptdyjet[kk][ir-1] = dndptdyjet[kk][ir-1] + weight / jet_absetamax / 2.0 / pT_diff[kk];
                        }
                    }
                }
            }
        }
    }

    if (even_loop_flag == 1) {
        fscanf(infile,"%s %s %lf %s %lf %s %lf %s %lf\n",&stemp1, &stemp2,
                       &sigmaGen, &stemp1, &sigmaErr, &stemp1, &weight,
                       &stemp1, &pThat);
    } else {
        sigmaGen = sigmaGen_ev;
    }
    fclose(infile);

    //output the results to file
    char output_filename[128];
    sprintf(output_filename,"Hard_jet_yield");
    ofstream output(output_filename);
    if( ! output.is_open() ) {
        cout << "cannot open output file:"<< endl
             << output_filename << endl;
        return -1;
    }
    output << " # " << "pT lower, pT upper, center pT, dN/dpTdeta, R = 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7 " << " sigmaGen " << " Wight_sum " << " Nevent " << endl; 
    for (int i = 0; i < length_pt; ++i) {
        output << pT_array[i] << "  " << pT_array[i+1] << "  " << (pT_array[i]+pT_array[i+1])/2 << "  " 
               << dndptdyjet[i][0] << "  " << dndptdyjet[i][1] << "  " << dndptdyjet[i][2] << "  " 
               << dndptdyjet[i][3] << "  " << dndptdyjet[i][4] << "  " << dndptdyjet[i][5] << "  " 
               << sigmaGen << "  " << weight_sum << "  " << count_event_number
               << endl;
    }
    output.close();
    return 0;
}
