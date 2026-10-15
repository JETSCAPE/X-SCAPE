#!/usr/bin/env python3
"""JETSCAPE user-XML templates for the EMD qhat study (bricks and hydro).

Used by make_configs.py, which names every configuration of the study.
"""
import sys

# JETSCAPE MATTER+LBT tune: virtuality-dependent HTL qhat (Type 5, A = 10, B = 100,
# alpha_s = 0.3, Q0 = 2 GeV), from the JETSCAPE Bayesian analysis, PRC 104, 024905 (2021)
MATTER_TUNE = """    <Matter>
      <in_vac> {in_vac} </in_vac> <brick_med> {brick} </brick_med> <brick_length> {L} </brick_length>
      <recoil_on> 0 </recoil_on> <Q0> {Q0} </Q0> <QS> {QS} </QS> <T0> 0.16 </T0> <hydro_Tc> 0.16 </hydro_Tc>
      <QhatParametrizationType> 5 </QhatParametrizationType> <qhat0> -2.0 </qhat0>
      <alphas> 0.3 </alphas> <qhatA> 10.0 </qhatA> <qhatB> 100.0 </qhatB>
    </Matter>
    <Lbt> <name> Lbt </name> <Q0> {Q0} </Q0> <in_vac> 0 </in_vac> <only_leading> 0 </only_leading>
      <hydro_Tc> 0.16 </hydro_Tc> <alphas> 0.3 </alphas> <run_alphas> 1 </run_alphas> </Lbt>
"""
MATTER_EMD = """    <Matter>
      <in_vac> {in_vac} </in_vac> <brick_med> {brick} </brick_med> <brick_length> {L} </brick_length>
      <recoil_on> 0 </recoil_on> <Q0> {Q0} </Q0> <QS> {QS} </QS> <T0> 0.16 </T0> <hydro_Tc> 0.16 </hydro_Tc>
      <QhatParametrizationType> 8 </QhatParametrizationType> <qhat0> 2.0 </qhat0>
      <alphas> 0.3 </alphas> <lambda_t> {lam} </lambda_t> <qhat_K> 1.0 </qhat_K> <conformal_Tref> {ctref} </conformal_Tref>
    </Matter>
    <MyEloss> <name> MyEloss </name> <Q0> {Q0} </Q0> <lambda_t> {lam} </lambda_t> <conformal_Tref> {ctref} </conformal_Tref>
      <radiation_on> {rad} </radiation_on> <broadening_on> {broad} </broadening_on> </MyEloss>
"""
MATTER_EMDV = MATTER_EMD.replace("<QhatParametrizationType> 8 </QhatParametrizationType>",
    "<QhatParametrizationType> 9 </QhatParametrizationType> <qhatA> 10.0 </qhatA> <qhatB> 100.0 </qhatB>")
# Option 1 ("basics"): stock MATTER+LBT with only the radiative qhat(T) swapped
# for EMD: MATTER Type 9 (EMD x Type-5 virtuality factor) + LBT qhat_model=1
# (radiative qhat from EMD; LBT elastic scattering unchanged)
LBT_EMD_BRICK = MATTER_EMDV.split("    <MyEloss>")[0] + """    <Lbt> <name> Lbt </name> <Q0> {Q0} </Q0> <in_vac> 0 </in_vac> <only_leading> 0 </only_leading>
      <hydro_Tc> 0.16 </hydro_Tc> <alphas> 0.3 </alphas> <run_alphas> 1 </run_alphas>
      <qhat_model> 1 </qhat_model> <lambda_t> {lam} </lambda_t> <qhat_K> 1.0 </qhat_K> </Lbt>
"""
# Route A with JETSCAPE's HTL qhat in both stages: MATTER Type 5 and MyEloss
# qhat_model=1, both times qhat_K = K (the HTL counterpart of MATTER_EMDV)
MATTER_HTLV = MATTER_EMDV.replace(
    "<QhatParametrizationType> 9 </QhatParametrizationType> <qhatA> 10.0 </qhatA> <qhatB> 100.0 </qhatB> <qhat0> 2.0 </qhat0>",
    "<QhatParametrizationType> 5 </QhatParametrizationType> <qhatA> 10.0 </qhatA> <qhatB> 100.0 </qhatB> <qhat0> -2.0 </qhat0>").replace(
    "<lambda_t> {lam} </lambda_t> <qhat_K> 1.0 </qhat_K> <conformal_Tref> {ctref} </conformal_Tref>",
    "<lambda_t> 1.0 </lambda_t> <qhat_K> {K} </qhat_K> <conformal_Tref> 0.0 </conformal_Tref>").replace(
    "<MyEloss> <name> MyEloss </name> <Q0> {Q0} </Q0> <lambda_t> {lam} </lambda_t> <conformal_Tref> {ctref} </conformal_Tref>",
    "<MyEloss> <name> MyEloss </name> <Q0> {Q0} </Q0> <lambda_t> 1.0 </lambda_t> <qhat_model> 1 </qhat_model> <qhat_K> {K} </qhat_K> <htl_alphas> 0.3 </htl_alphas>")
assert MATTER_HTLV.count("{K}") == 2
VAC = """    <Matter> <in_vac> 1 </in_vac> <Q0> 1.0 </Q0> <QS> 1.0 </QS> <recoil_on> 0 </recoil_on> </Matter>
"""

def eloss(model, **kw):
    kw.setdefault("in_vac", 0); kw.setdefault("QS", 1.0); kw.setdefault("ctref", 0.0); kw.setdefault("Q0", 2.0); kw.setdefault("lam", 1.0); kw.setdefault("K", 1.0)
    kw.setdefault("rad", 1); kw.setdefault("broad", 1); kw.setdefault("brick", 0); kw.setdefault("L", 5.0)
    return {"vac": VAC, "tune": MATTER_TUNE, "emd": MATTER_EMD, "emdv": MATTER_EMDV, "htlv": MATTER_HTLV, "lbtemd": LBT_EMD_BRICK}[model].format(**kw)

def brick(model, out, T=0.3, L=4.0, pid=1, pt=100, nev=200, seed=1, bjorken=False, t0=0.0, **kw):
    return f"""<?xml version="1.0"?>
<jetscape>
  <nEvents> {nev} </nEvents> <vlevel> 0 </vlevel> <nEvents_printout> 50 </nEvents_printout>
  <outputFilename>{out}</outputFilename>
  <JetScapeWriterAscii> off </JetScapeWriterAscii>
  <JetScapeWriterFinalStatePartonsAscii> on </JetScapeWriterFinalStatePartonsAscii>
  <Random> <seed>{seed}</seed> </Random>
  <Hard> <PGun> <name>PGun</name> <pT>{pt}</pT> <parID>{pid}</parID> </PGun> </Hard>
  <Hydro> <Brick bjorken_expansion_on="{"true" if bjorken else "false"}" start_time="{t0}"> <name>Brick</name>
    <T>{T}</T> <v_x>0.0</v_x> <v_y>0.0</v_y> <v_z>0.0</v_z> </Brick> </Hydro>
  <Eloss> <maxT>30</maxT> <tStart> {t0} </tStart>
{eloss(model, brick=1, L=L, **kw)}  </Eloss>
</jetscape>
"""

if __name__ == "__main__":
    print(brick(sys.argv[1], sys.argv[2]))

# pT-hat biased sampling: events carry weights ~ (pT-hat/10)^-4 compensation
BIAS = """
    <LinesToRead>
      PhaseSpace:bias2Selection = on
      PhaseSpace:bias2SelectionPow = 4.0
      PhaseSpace:bias2SelectionRef = 10.0
    </LinesToRead>"""

AA_MEDIUM = """  <IS>
    <grid_max_x> 15 </grid_max_x> <grid_max_y> 15 </grid_max_y> <grid_max_z> 0 </grid_max_z>
    <grid_step_x> {dx} </grid_step_x> <grid_step_y> {dx} </grid_step_y>
    <Trento>
      <PhysicsInputs projectile='Pb' target='Pb' sqrts='5020' cross-section='6.7' normalization='18.38'> </PhysicsInputs>
      <CutInputs centrality-low='{clow}' centrality-high='{chigh}'> </CutInputs>
      <TransInputs reduced-thickness='0.007' fluctuation="0.918" nucleon-width="0.956" nucleon-min-dist="1.27"> </TransInputs>
      <LongiInputs mean-coeff="1.0" std-coeff="3.0" skew-coeff="0.0" skew-type="1" jacobian="0.8"> </LongiInputs>
    </Trento>
  </IS>
  <Preequilibrium>
    <tau0>0.0</tau0> <tauj>0.1</tauj> <taus>1.16</taus> <evolutionInMemory>1</evolutionInMemory>
    <FreestreamMilne> <name>FreestreamMilne</name> <freestream_input_file>freestream_input</freestream_input_file> <E_DEP_FS>0</E_DEP_FS> </FreestreamMilne>
  </Preequilibrium>
  <Hydro>
    <MUSIC>
      <beastMode>0</beastMode> <output_evolution_to_file>0</output_evolution_to_file>
      <output_evolution_to_memory>1</output_evolution_to_memory> <surface_in_memory>1</surface_in_memory>
      <Initial_time_tau_0>1.16</Initial_time_tau_0>
      <T_dependent_Shear_to_S_ratio>2</T_dependent_Shear_to_S_ratio>
      <shear_viscosity_2_curv>-0.48</shear_viscosity_2_curv> <shear_viscosity_2_slope>1.11</shear_viscosity_2_slope>
      <shear_viscosity_2_min>0.081</shear_viscosity_2_min>
      <temperature_dependent_bulk_viscosity>2</temperature_dependent_bulk_viscosity>
      <bulk_viscosity_2_normalisation>0.052</bulk_viscosity_2_normalisation>
      <bulk_viscosity_2_width_in_GeV>0.022</bulk_viscosity_2_width_in_GeV>
      <bulk_viscosity_2_peak_in_GeV>0.183</bulk_viscosity_2_peak_in_GeV>
      <freezeout_temperature>0.151</freezeout_temperature>
    </MUSIC>
  </Hydro>
"""

def collision(model, out, pthat, nev, seed, nreuse=None, clow=0, chigh=10, dx=0.3, bias=True, **kw):
    """pp (model='vac') or Pb-Pb (model in tune/emd/emdv) jets at 5.02 TeV."""
    aa = model != "vac"
    reuse = f"<setReuseHydro> true </setReuseHydro> <nReuseHydro> {nreuse} </nReuseHydro>" if aa else ""
    return f"""<?xml version="1.0"?>
<jetscape>
  <nEvents> {nev} </nEvents> <vlevel> 0 </vlevel> <nEvents_printout> 10 </nEvents_printout>
  {reuse}
  <outputFilename>{out}</outputFilename>
  <JetScapeWriterAscii> off </JetScapeWriterAscii>
  <JetScapeWriterFinalStateHadronsAscii> on </JetScapeWriterFinalStateHadronsAscii>
  <Random> <seed>{seed}</seed> </Random>
{AA_MEDIUM.format(clow=clow, chigh=chigh, dx=dx) if aa else ""}
  <Hard> <PythiaGun> <name>PythiaGun</name> <FSR_on>0</FSR_on> <pTHatMin>{pthat[0]}</pTHatMin> <pTHatMax>{pthat[1]}</pTHatMax> <eCM>5020</eCM>{BIAS if bias else ""} </PythiaGun> </Hard>
  <Eloss> <tStart> 0.6 </tStart>
{eloss(model, **kw)}  </Eloss>
  <JetHadronization> <name>colorless</name> </JetHadronization>
</jetscape>
"""


# --------------------------------------------------------------------------
# Published-JETSCAPE-style runs (arXiv:2204.01163 options) for LHC and RHIC.
# All AA models share: maxT 150, colorless hadronization with take_recoil and
# eCMforHadronization = sqrt(s)/2, pi0 not decayed, EPS09 nPDF at the LHC
# (Pythia ships EPS09 for Pb only, so none for Au).

SYSTEMS = {
    "LHC": dict(ecm=5020, pthat=(10, 500), bias_ref=10.0, npdf_beam=100822080, medium=AA_MEDIUM),
    "RHIC": dict(ecm=200, pthat=(7, 100), bias_ref=7.0, npdf_beam=None, medium="""  <IS>
    <grid_max_x> 15 </grid_max_x> <grid_max_y> 15 </grid_max_y> <grid_max_z> 0 </grid_max_z>
    <grid_step_x> {dx} </grid_step_x> <grid_step_y> {dx} </grid_step_y>
    <Trento>
      <PhysicsInputs projectile='Au' target='Au' sqrts='200' cross-section='4.2' normalization='5.73'> </PhysicsInputs>
      <CutInputs centrality-low='{clow}' centrality-high='{chigh}'> </CutInputs>
      <TransInputs reduced-thickness='0.063' fluctuation="1.05" nucleon-width="1.12" nucleon-min-dist="1.4"> </TransInputs>
      <LongiInputs mean-coeff="1.0" std-coeff="3.0" skew-coeff="0.0" skew-type="1" jacobian="0.8"> </LongiInputs>
    </Trento>
  </IS>
  <Preequilibrium>
    <tau0>0.0</tau0> <tauj>0.1</tauj> <taus>1.46</taus> <evolutionInMemory>1</evolutionInMemory>
    <FreestreamMilne> <name>FreestreamMilne</name> <freestream_input_file>freestream_input</freestream_input_file> </FreestreamMilne>
  </Preequilibrium>
  <Hydro>
    <MUSIC>
      <beastMode>0</beastMode> <output_evolution_to_file>0</output_evolution_to_file>
      <output_evolution_to_memory>1</output_evolution_to_memory> <surface_in_memory>1</surface_in_memory>
      <Initial_time_tau_0>1.46</Initial_time_tau_0>
      <T_dependent_Shear_to_S_ratio>3</T_dependent_Shear_to_S_ratio>
      <!-- JETSCAPE Au-Au 200 GeV MAP (arXiv:2011.01430), v2.1 parameter names -->
      <shear_viscosity_3_T_kink_in_GeV>0.223</shear_viscosity_3_T_kink_in_GeV> <shear_viscosity_3_low_T_slope_in_GeV>-0.776</shear_viscosity_3_low_T_slope_in_GeV>
      <shear_viscosity_3_high_T_slope_in_GeV>0.37</shear_viscosity_3_high_T_slope_in_GeV> <shear_viscosity_3_at_kink>0.096</shear_viscosity_3_at_kink>
      <temperature_dependent_bulk_viscosity>3</temperature_dependent_bulk_viscosity>
      <bulk_viscosity_3_max>0.133</bulk_viscosity_3_max> <bulk_viscosity_3_T_peak_in_GeV>0.12</bulk_viscosity_3_T_peak_in_GeV>
      <bulk_viscosity_3_width_in_GeV>0.072</bulk_viscosity_3_width_in_GeV> <bulk_viscosity_3_lambda_asymm>-0.122</bulk_viscosity_3_lambda_asymm>
      <freezeout_temperature>0.136</freezeout_temperature>
    </MUSIC>
  </Hydro>
"""),
}

MATTER5 = """    <Matter>
      <in_vac> 0 </in_vac> <brick_med> 0 </brick_med> <recoil_on> {mrecoil} </recoil_on>
      <Q0> 2.0 </Q0> <QS> 1.0 </QS> <T0> 0.16 </T0> <hydro_Tc> 0.16 </hydro_Tc> <vir_factor> 0.25 </vir_factor>
      <QhatParametrizationType> {mtype} </QhatParametrizationType> <qhat0> {qhat0} </qhat0>
      <alphas> 0.3 </alphas> <qhatA> 10.0 </qhatA> <qhatB> 100.0 </qhatB>
      <lambda_t> {lam} </lambda_t> <qhat_K> {K} </qhat_K>
    </Matter>
"""
LBT = """    <Lbt> <name> Lbt </name> <Q0> 2.0 </Q0> <in_vac> 0 </in_vac> <only_leading> 0 </only_leading>
      <hydro_Tc> 0.16 </hydro_Tc> <alphas> 0.3 </alphas> <run_alphas> 1 </run_alphas> </Lbt>
"""
LBT_EMD = """    <Lbt> <name> Lbt </name> <Q0> 2.0 </Q0> <in_vac> 0 </in_vac> <only_leading> 0 </only_leading>
      <hydro_Tc> 0.16 </hydro_Tc> <alphas> 0.3 </alphas> <run_alphas> 1 </run_alphas>
      <qhat_model> 1 </qhat_model> <lambda_t> {lam} </lambda_t> <qhat_K> {K} </qhat_K> </Lbt>
"""
MYELOSS = """    <MyEloss> <name> MyEloss </name> <Q0> 2.0 </Q0> <qhat_model> {qmodel} </qhat_model>
      <lambda_t> {lam} </lambda_t> <qhat_K> {K} </qhat_K> <htl_alphas> 0.3 </htl_alphas> </MyEloss>
"""


HAD_PUB = """<JetHadronization> <name>colorless</name> <take_recoil>1</take_recoil> <eCMforHadronization>{ecm:g}</eCMforHadronization>
    <LinesToRead>
      111:mayDecay = off
    </LinesToRead> </JetHadronization>"""
HAD_DEFAULT = """<JetHadronization> <name>colorless</name> <eCMforHadronization>{ecm:g}</eCMforHadronization> </JetHadronization>"""


def eloss2(model, lam=9.0, K=1.0):
    """tune: MATTER+LBT as in the first runs (MATTER recoil off);
    tunepub: published MATTER+LBT (MATTER recoil on); emdv: MATTER(Type 9)+
    MyEloss(EMD); lbtemd: MATTER(Type 9)+LBT with EMD radiative qhat
    (option 1, everything else as in tune); htlv: MATTER(Type 5)+MyEloss(HTL). lam/K scale qhat in both
    stages."""
    if model == "vac":
        return VAC
    if model == "tune":
        return MATTER5.format(mrecoil=0, mtype=5, qhat0=-2.0, lam=1.0, K=1.0) + LBT
    if model == "tunepub":
        return MATTER5.format(mrecoil=1, mtype=5, qhat0=-2.0, lam=1.0, K=1.0) + LBT
    if model == "emdv":
        return MATTER5.format(mrecoil=0, mtype=9, qhat0=2.0, lam=lam, K=K) + MYELOSS.format(qmodel=0, lam=lam, K=K)
    if model == "lbtemd":
        return MATTER5.format(mrecoil=0, mtype=9, qhat0=2.0, lam=lam, K=K) + LBT_EMD.format(lam=lam, K=K)
    if model == "htlv":
        return MATTER5.format(mrecoil=0, mtype=5, qhat0=-2.0, lam=1.0, K=K) + MYELOSS.format(qmodel=1, lam=1.0, K=K)
    raise ValueError(model)


def collision2(system, model, out, nev, seed, nreuse=None, clow=0, chigh=10, dx=0.3,
               maxT=20, npdf=False, published_hadronization=False, **kw):
    """Defaults reproduce the settings of the first production runs (maxT 20,
    no nPDF, default hadronization); the published arXiv:2204.01163 options
    are maxT=150, npdf=True, published_hadronization=True."""
    S = SYSTEMS[system]
    aa = model != "vac"
    lines = [f"PhaseSpace:bias2Selection = on", "PhaseSpace:bias2SelectionPow = 4.0",
             f"PhaseSpace:bias2SelectionRef = {S['bias_ref']}"]
    if aa and npdf and S["npdf_beam"]:
        b = S["npdf_beam"]
        lines += ["PDF:useHardNPDFA = on", "PDF:useHardNPDFB = on", "PDF:nPDFSetA = 1", "PDF:nPDFSetB = 1",
                  f"PDF:nPDFBeamA = {b}", f"PDF:nPDFBeamB = {b}"]
    lt = "\n".join("      " + l for l in lines)
    ecm_had = S["ecm"] / 2 if aa else S["ecm"]
    reuse = f"<setReuseHydro> true </setReuseHydro> <nReuseHydro> {nreuse} </nReuseHydro>" if aa else ""
    medium = S["medium"].format(clow=clow, chigh=chigh, dx=dx) if aa else ""
    return f"""<?xml version="1.0"?>
<jetscape>
  <nEvents> {nev} </nEvents> <vlevel> 0 </vlevel> <nEvents_printout> 10 </nEvents_printout>
  {reuse}
  <outputFilename>{out}</outputFilename>
  <JetScapeWriterAscii> off </JetScapeWriterAscii>
  <JetScapeWriterFinalStateHadronsAscii> on </JetScapeWriterFinalStateHadronsAscii>
  <Random> <seed>{seed}</seed> </Random>
{medium}
  <Hard> <PythiaGun> <name>PythiaGun</name> <FSR_on>0</FSR_on> <pTHatMin>{S['pthat'][0]}</pTHatMin> <pTHatMax>{S['pthat'][1]}</pTHatMax> <eCM>{S['ecm']}</eCM>
    <LinesToRead>
{lt}
    </LinesToRead> </PythiaGun> </Hard>
  <Eloss> <tStart> 0.6 </tStart> <maxT> {maxT} </maxT>
{eloss2(model, **kw)}  </Eloss>
  {HAD_PUB.format(ecm=ecm_had) if published_hadronization else HAD_DEFAULT.format(ecm=S["ecm"])}
</jetscape>
"""
