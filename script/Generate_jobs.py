#!/usr/bin/env python3

"""This script generate xml files."""

from subprocess import call
import sys
from os import path, mkdir
import shutil
import subprocess
import argparse
from math import ceil
from glob import glob
import random

# X-SCAPE
xscape_dict = {
    'eCM': 5020.,
    'initial_mode': 13,
    'Q0': 2.0,
    'pTHatMin': 4.0,
    'pTmin': 3.0,
    'alphaSvalue': 0.16,
    'vir_factor': 0.25,
    'ylossParam4At2': 1.5,
    'ylossParam4At4': 1.5,
    'ylossParam4At6': 2.0,
    'ylossParam4At10': 6.0,
    'ylossParam4var': 0.5,
    'string_source_sigma_x': 0.5,
    'string_source_sigma_eta': 0.6,
    'stringPreEqFlowFactor': 0.18,
    'shear_viscosity_3_T_kink': 0.15,
    'bulk_viscosity_3_max': 0.1,
    'esw': 0.3
}

def generate_xml(random_seed):
    xml_files = '''<?xml version="1.0"?>

<jetscape>

  <enableAutomaticTaskListDetermination> false </enableAutomaticTaskListDetermination>

  <setReuseHydro> false </setReuseHydro>
  <!-- <nReuseHydro> 10 </nReuseHydro> -->

  <vlevel> 4 </vlevel>

  <nEvents> 10 </nEvents>

  <Random>
    <seed> {seed} </seed>
  </Random>

  <JetScapeWriterAscii> on </JetScapeWriterAscii>

  <!-- Hard Process -->
  <Hard>
    <PythiaGun>
      <pTHatMin>{pTHatMin}</pTHatMin>
      <pTHatMax>-1</pTHatMax>
      <eCM>{eCM}</eCM>
      <LinesToRead>
        PhaseSpace:bias2Selection = on
        PhaseSpace:bias2SelectionPow = 4
        PhaseSpace:bias2SelectionRef = 10
        SigmaProcess:alphaSvalue = {alphaSvalue}
        MultipartonInteractions:alphaSvalue = {alphaSvalue}
        MultipartonInteractions:pTmin = {pTmin}
      </LinesToRead>
      <useHybridHad>1</useHybridHad>
    </PythiaGun>
  </Hard>

  <!--Preequilibrium Dynamics Module -->
  <Preequilibrium>
    <NullPreDynamics> </NullPreDynamics>
  </Preequilibrium>

  <!-- Initial condition  Module  -->
  <IS>
    <MCGlauber> 
        <ylossParam4At2> {yloss2} </ylossParam4At2>
        <ylossParam4At4> {yloss4} </ylossParam4At4>
        <ylossParam4At6> {yloss6} </ylossParam4At6>
        <ylossParam4At10> {yloss10} </ylossParam4At10>
        <ylossParam4var> {ylossvar} </ylossParam4var>
    </MCGlauber>
  </IS>

  <!-- Hydro  Module  -->
  <Hydro>
        <MUSIC>
            <InitialProfile> {initial_mode} </InitialProfile>
            <string_source_sigma_x> {sigma_x} </string_source_sigma_x>
            <string_source_sigma_eta> {sigma_eta} </string_source_sigma_eta>
            <stringPreEqFlowFactor> {preflow} </stringPreEqFlowFactor>
            <output_evolution_to_file>0</output_evolution_to_file>
            <store_hydro_info_in_memory>0</store_hydro_info_in_memory>
            <T_dependent_Shear_to_S_ratio> 3 </T_dependent_Shear_to_S_ratio>
            <temperature_dependent_bulk_viscosity>3</temperature_dependent_bulk_viscosity>
            <bulk_viscosity_3_max> {zeta_s} </bulk_viscosity_3_max>
            <shear_viscosity_3_at_kink> {eta_s} </shear_viscosity_3_at_kink>
            <freezeout_temperature> {esw} </freezeout_temperature>
        </MUSIC>
  </Hydro>

  <!--Eloss Modules -->
  <Eloss>
    <maxT>50</maxT>
    <Matter>
      <Q0> {Q0} </Q0>
      <in_vac> 1 </in_vac>
      <useHybridHad>1</useHybridHad>
      <vir_factor> {vir_factor} </vir_factor>
      <recoil_on> 0 </recoil_on>
      <broadening_on> 0 </broadening_on>
      <brick_med> 0 </brick_med>
    </Matter>
  </Eloss>

  <!-- Jet Hadronization Module -->
  <JetHadronization>
    <name>colorless</name>
  </JetHadronization>

  <!-- Particlization Module  -->
  <SoftParticlization>
    <iSS>
        <number_of_repeated_sampling> 1000 </number_of_repeated_sampling>
        <Perform_resonance_decays> 1 </Perform_resonance_decays>
    </iSS>
  </SoftParticlization>

</jetscape>
'''.format(seed = random_seed, eCM = xscape_dict['eCM'], initial_mode = xscape_dict['initial_mode'], Q0 = xscape_dict['Q0'],
           pTHatMin = xscape_dict['pTHatMin'], pTmin = xscape_dict['pTmin'], alphaSvalue = xscape_dict['alphaSvalue'], 
           vir_factor = xscape_dict['vir_factor'], yloss2 = xscape_dict['ylossParam4At2'], 
           yloss4 = xscape_dict['ylossParam4At4'], yloss6 = xscape_dict['ylossParam4At6'], yloss10 = xscape_dict['ylossParam4At10'],
           ylossvar = xscape_dict['ylossParam4var'], sigma_x = xscape_dict['string_source_sigma_x'], 
           sigma_eta = xscape_dict['string_source_sigma_eta'], preflow = xscape_dict['stringPreEqFlowFactor'],
           eta_s = xscape_dict['shear_viscosity_3_T_kink'], zeta_s = xscape_dict['bulk_viscosity_3_max'], esw = xscape_dict['esw'])
    xml_name = "jetscape_user_iMATTERMCGlauberMUSIC_test.xml"
    with open(xml_name, 'w') as fout:
        fout.write(xml_files)
    #call(['mv', xml_name, './'])
def update_parameters_bayesian(bayes_file):
    parfile = open(bayes_file, "r")
    for line in parfile:
        key, val = line.split()
        if key in xscape_dict.keys():
            xscape_dict[key] = float(val)

if __name__=='__main__':
    parser = argparse.ArgumentParser(
        description='\U0000269B Welcome to the X-SCAPE framework',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-seed',
                        '--random_seed',
                        metavar='',
                        type=int,
                        default='-1',
                        help='Random Seed (-1: according to system time)')
    parser.add_argument('-b',
                        '--bayes_file',
                        metavar='',
                        type=str,
                        default='',
                        help='parameters from bayesian analysis')
    args = parser.parse_args()
    if args.bayes_file != "":
        args.bayes_file = path.join(path.abspath("."), args.bayes_file)
        update_parameters_bayesian(args.bayes_file)
        #shutil.copy(args.bayes_file, working_folder_name)
    generate_xml(args.random_seed)
