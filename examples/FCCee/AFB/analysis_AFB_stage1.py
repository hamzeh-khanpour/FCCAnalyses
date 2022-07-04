# Forward-Backward-Asymmetry June 2022

#Mandatory: List of processes
processList = {
    'p8_ee_Zbb_ecm91':{'fraction':0.001},#Run the full statistics in one output file named <outputDir>/p8_ee_Zbb_ecm91.root
#    'p8_ee_WW_ecm240':{'fraction':0.5, 'chunks':2}, #Run 50% of the statistics in two files named <outputDir>/p8_ee_WW_ecm240/chunk<N>.root
#    'p8_ee_ZH_ecm240':{'fraction':0.2, 'output':'p8_ee_ZH_ecm240_out'} #Run 20% of the statistics in one file named <outputDir>/p8_ee_ZH_ecm240_out.root (example on how to change the output name)
}

#Mandatory: Production tag when running over EDM4Hep centrally produced events, this points to the yaml files for getting sample statistics
prodTag     = "FCCee/spring2021/IDEA/" 

#Optional: output directory, default is local running directory 
outputDir   = "outputs/FCCee/AFB/stage1" 

#Optional: analysisName, default is ""
#analysisName = "My Analysis"

#Optional: ncpus, default is 4
nCPUS       = 8

#Optional running on HTCondor, default is False
#runBatch    = False

#Optional batch queue name when running on HTCondor, default is workday
#batchQueue = "longlunch"

#Optional computing account when running on HTCondor, default is group_u_FCC.local_gen
#compGroup = "group_u_FCC.local_gen"

#Optional test file
#testFile ="root://eospublic.cern.ch//eos/experiment/fcc/ee/generation/DelphesEvents/spring2021/IDEA/p8_ee_ZH_ecm240/events_101027117.root"

#Mandatory: RDFanalysis class where the use defines the operations on the TTree
class RDFanalysis():

    #__________________________________________________________
    #Mandatory: analysers funtion to define the analysers to process, please make sure you return the last dataframe, in this example it is df2
    def analysers(df):
        df2 = (
            df
               .Alias("Electron0", "Electron#0.index")
               .Define("electrons", "ReconstructedParticle::get(Electron0, ReconstructedParticles)")               
               
               .Alias("Muon0", "Muon#0.index")
               .Define("muons", "ReconstructedParticle::get(Muon0, ReconstructedParticles)")



               .Define("n_electrons", "ReconstructedParticle::get_n(electrons)")
               .Define("electron_theta", "ReconstructedParticle::get_theta(electrons)")           
               .Define("electron_e", "ReconstructedParticle::get_e(electrons)")
               .Define("electron_px", "ReconstructedParticle::get_px(electrons)")
               .Define("electron_py", "ReconstructedParticle::get_py(electrons)")
               .Define("electron_pz", "ReconstructedParticle::get_pz(electrons)")
               .Define("electron_pt", "ReconstructedParticle::get_pt(electrons)")
               .Define("electron_eta", "ReconstructedParticle::get_eta(electrons)")
               .Define("electron_phi", "ReconstructedParticle::get_phi(electrons)")
               .Define("electron_mass", "ReconstructedParticle::get_mass(electrons)")
               .Define("electron_charge", "ReconstructedParticle::get_charge(electrons)")               
               

               .Define("n_muons", "ReconstructedParticle::get_n(muons)")
               .Define("muon_theta", "ReconstructedParticle::get_theta(muons)")
               .Define("muon_e", "ReconstructedParticle::get_e(muons)")
               .Define("muon_px", "ReconstructedParticle::get_px(muons)")
               .Define("muon_py", "ReconstructedParticle::get_py(muons)")
               .Define("muon_pz", "ReconstructedParticle::get_pz(muons)")
               .Define("muon_pt", "ReconstructedParticle::get_pt(muons)")
               .Define("muon_eta", "ReconstructedParticle::get_eta(muons)")
               .Define("muon_phi", "ReconstructedParticle::get_phi(muons)")
               .Define("muon_mass", "ReconstructedParticle::get_mass(muons)")
               .Define("muon_charge", "ReconstructedParticle::get_charge(muons)")               



               .Define("leptons", "ReconstructedParticle::merge(electrons, muons)")
               .Define("n_leptons", "ReconstructedParticle::get_n(leptons)")
               .Define("lepton_theta", "ReconstructedParticle::get_theta(leptons)")
               .Define("lepton_e", "ReconstructedParticle::get_e(leptons)")
               .Define("lepton_px", "ReconstructedParticle::get_px(leptons)")
               .Define("lepton_py", "ReconstructedParticle::get_py(leptons)")
               .Define("lepton_pz", "ReconstructedParticle::get_pz(leptons)")
               .Define("lepton_pt", "ReconstructedParticle::get_pt(leptons)")
               .Define("lepton_eta", "ReconstructedParticle::get_eta(leptons)")
               .Define("lepton_phi", "ReconstructedParticle::get_phi(leptons)")
               .Define("lepton_mass", "ReconstructedParticle::get_mass(leptons)")
               .Define("lepton_charge", "ReconstructedParticle::get_charge(leptons)") 



#################--------------------------------------------------------------------------


               .Alias("MCRecoAssociations0", "MCRecoAssociations#0.index")
               .Alias("MCRecoAssociations1", "MCRecoAssociations#1.index")

               .Alias("Particle0", "Particle#0.index")
               .Alias("Particle1", "Particle#1.index")
               
               .Define("MC_px",         "MCParticle::get_px(Particle)")
               .Define("MC_py",         "MCParticle::get_py(Particle)")
               .Define("MC_pz",         "MCParticle::get_pz(Particle)")
               .Define("MC_p",          "MCParticle::get_p(Particle)")
               .Define("MC_e",          "MCParticle::get_e(Particle)")
               .Define("MC_pdg",        "MCParticle::get_pdg(Particle)")
               .Define("MC_theta",      "MCParticle::get_theta(Particle)")
               .Define("MC_charge",     "MCParticle::get_charge(Particle)")               
               .Define("MC_mass",       "MCParticle::get_mass(Particle)")
               .Define("MC_status",     "MCParticle::get_genStatus(Particle)")
               .Define("MC_vertex_x",   "MCParticle::get_vertex_x(Particle)")
               .Define("MC_vertex_y",   "MCParticle::get_vertex_y(Particle)")
               .Define("MC_vertex_z",   "MCParticle::get_vertex_z(Particle)")
#               .Define("MC_vertex_z",   "FCCAnalyses::MCParticle::get_vertex_z(Particle)")
               
               

               .Define("MC_PDG", "FCCAnalyses::MCParticle::get_pdg(Particle)")
               .Define("MC_n",   "int(MC_PDG.size())")
               .Define("MC_M1",  "myUtils::get_MCMother1(Particle,Particle0)")
               .Define("MC_M2",  "myUtils::get_MCMother2(Particle,Particle0)")
               .Define("MC_D1",  "myUtils::get_MCDaughter1(Particle,Particle1)")
               .Define("MC_D2",  "myUtils::get_MCDaughter2(Particle,Particle1)")



#################--------------------------------------------------------------------------



                .Define("selected_jets", "ReconstructedParticle::sel_p(0.0)(Jet)")          # select only jets with a pT > 0 GeV  
                .Define("n_jets",        "ReconstructedParticle::get_n(selected_jets)")     # count how many jets are in the event in total
                .Define("seljet_p",      "ReconstructedParticle::get_p(selected_jets)")     # momentum p
		        .Define("seljet_theta",  "ReconstructedParticle::get_theta(selected_jets)") # theta
                .Define("seljet_e",   "ReconstructedParticle::get_e(selected_jets)")        # Energy of jet
                .Define("seljet_px",  "ReconstructedParticle::get_px(selected_jets)")       # create branch with jet px
                .Define("seljet_py",  "ReconstructedParticle::get_py(selected_jets)")       # create branch with jet py
                .Define("seljet_pz",  "ReconstructedParticle::get_pz(selected_jets)")       # create branch with jet pz 
                .Define("seljet_eta", "ReconstructedParticle::get_eta(selected_jets)")      # create branch with jet eta
                
 

               .Alias("Jet3","Jet#3.index") 
               .Define("JET_btag", "ReconstructedParticle::getJet_btag(Jet3, ParticleIDs, ParticleIDs_0)")
               .Define("EVT_nbtag", "ReconstructedParticle::getJet_ntags(JET_btag)")



               .Define("RecoParticles_wo_electrons",  "ReconstructedParticle::remove(ReconstructedParticles,electrons)")
               .Define("RecoParticles_wo_leptons",  "ReconstructedParticle::remove(RecoParticles_wo_electrons,muons)")


               .Define("RP_px",         "ReconstructedParticle::get_px(RecoParticles_wo_leptons)")
               .Define("RP_py",         "ReconstructedParticle::get_py(RecoParticles_wo_leptons)")
               .Define("RP_pz",         "ReconstructedParticle::get_pz(RecoParticles_wo_leptons)")
               .Define("RP_p",          "ReconstructedParticle::get_p(RecoParticles_wo_leptons)")
               .Define("RP_e",          "ReconstructedParticle::get_e(RecoParticles_wo_leptons)")
               .Define("RP_charge",     "ReconstructedParticle::get_charge(RecoParticles_wo_leptons)")
               .Define("RP_mass",       "ReconstructedParticle::get_mass(RecoParticles_wo_leptons)")
               
               .Define("pseudo_jets",   "JetClusteringUtils::set_pseudoJets_xyzm(RP_px, RP_py, RP_pz, RP_mass)")



#################--------------------------------------------------------------------------




             #=================================================================================
                    # re-clustering the jet, Generalised-kt for e+e- cambridge (clustering_ee_genkt_E Scheme)-Exclusive jet selection exactly 2 jets
             #===================================================================================
              
               .Define("tw_FCCAnalysesJets_ee_genkt_ES","JetClustering::clustering_ee_genkt(0.5, 2, 2, 1, 0, 0)(pseudo_jets)")  #get the jets out of the struct
               .Define("tw_jets_ee_genkt_ES",           "JetClusteringUtils::get_pseudoJets(tw_FCCAnalysesJets_ee_genkt_ES)")   #get the jets constituents out of the struct
               .Define("tw_jetconstituents_ee_genkt_ES","JetClusteringUtils::get_constituents(tw_FCCAnalysesJets_ee_genkt_ES)")
               .Define("tw_jets_ee_genkt_ES_e",         "JetClusteringUtils::get_e(tw_jets_ee_genkt_ES)")
               .Define("tw_jets_ee_genkt_ES_px",        "JetClusteringUtils::get_px(tw_jets_ee_genkt_ES)")
               .Define("tw_jets_ee_genkt_ES_py",        "JetClusteringUtils::get_py(tw_jets_ee_genkt_ES)")
               .Define("tw_jets_ee_genkt_ES_pz",        "JetClusteringUtils::get_pz(tw_jets_ee_genkt_ES)")
               .Define("tw_jets_ee_genkt_ES_theta",     "JetClusteringUtils::get_theta(tw_jets_ee_genkt_ES)")
               .Define("tw_jets_ee_genkt_ES_charge",    "JetClusteringUtils::get_charge(tw_jets_ee_genkt_ES)")               
               .Define("tw_jets_ee_genkt_ES_flavour",   "JetTaggingUtils::get_flavour(tw_jets_ee_genkt_ES, Particle)")
               .Define("tw_jets_ee_genkt_ES_btag_true", "JetTaggingUtils::get_btag(tw_jets_ee_genkt_ES_flavour, 1.0)")
               .Define("tw_jets_ee_genkt_ES_btag",      "JetTaggingUtils::get_btag(tw_jets_ee_genkt_ES_flavour, 0.80)")
               .Define("tw_jets_ee_genkt_ES_ctag",      "JetTaggingUtils::get_ctag(tw_jets_ee_genkt_ES_flavour, 0.10)") 
               
               

               
               #====================================================================================================================
                    # re-clustering the jet, Generalised-kt for e+e- cambridge (clustering_ee_genkt_E Scheme)-Exclusive jet selection exactly 3 jets
               #======================================================================================================================              
              
               .Define("th_FCCAnalysesJets_ee_genkt_ES","JetClustering::clustering_ee_genkt(0.5, 3, 3, 1, 0, 0)(pseudo_jets)")  #get the jets out of the struct
               .Define("th_jets_ee_genkt_ES",           "JetClusteringUtils::get_pseudoJets(th_FCCAnalysesJets_ee_genkt_ES)")   #get the jets constituents out of the struct
               .Define("th_jetconstituents_ee_genkt_ES","JetClusteringUtils::get_constituents(th_FCCAnalysesJets_ee_genkt_ES)")
               .Define("th_jets_ee_genkt_ES_e",         "JetClusteringUtils::get_e(th_jets_ee_genkt_ES)")
               .Define("th_jets_ee_genkt_ES_px",        "JetClusteringUtils::get_px(th_jets_ee_genkt_ES)")
               .Define("th_jets_ee_genkt_ES_py",        "JetClusteringUtils::get_py(th_jets_ee_genkt_ES)")
               .Define("th_jets_ee_genkt_ES_pz",        "JetClusteringUtils::get_pz(th_jets_ee_genkt_ES)")
               .Define("th_jets_ee_genkt_ES_theta",     "JetClusteringUtils::get_theta(th_jets_ee_genkt_ES)")
               .Define("th_jets_ee_genkt_ES_charge",    "JetClusteringUtils::get_charge(th_jets_ee_genkt_ES)")
               .Define("th_jets_ee_genkt_ES_flavour",   "JetTaggingUtils::get_flavour(th_jets_ee_genkt_ES, Particle)")
               .Define("th_jets_ee_genkt_ES_btag_true", "JetTaggingUtils::get_btag(th_jets_ee_genkt_ES_flavour, 1.0)")
               .Define("th_jets_ee_genkt_ES_btag",      "JetTaggingUtils::get_btag(th_jets_ee_genkt_ES_flavour, 0.80)")
               .Define("th_jets_ee_genkt_ES_ctag",      "JetTaggingUtils::get_ctag(th_jets_ee_genkt_ES_flavour, 0.10)")
              



#################--------------------------------------------------------------------------



                 #MET
                 #=====================================================
                .Define("MET",   "ReconstructedParticle::get_p(MissingET)")  #absolute value of MET
		        .Define("MET_x", "ReconstructedParticle::get_px(MissingET)") #x-component of MET
		        .Define("MET_y", "ReconstructedParticle::get_py(MissingET)") #y-component of MET
                .Define("MET_z", "ReconstructedParticle::get_pz(MissingET)") #z-component of MET



#################--------------------------------------------------------------------------



               .Define("RP_px_all_objects",         "ReconstructedParticle::get_px(ReconstructedParticles)")
               .Define("RP_py_all_objects",         "ReconstructedParticle::get_py(ReconstructedParticles)")
               .Define("RP_pz_all_objects",         "ReconstructedParticle::get_pz(ReconstructedParticles)")
               .Define("RP_e_all_objects",          "ReconstructedParticle::get_e(ReconstructedParticles)")


               .Define('EVT_thrust',      'Algorithms::minimize_thrust("Minuit2","Migrad")(RP_px_all_objects, RP_py_all_objects, RP_pz_all_objects)')
               .Define('EVT_thrust_val',  'EVT_thrust.at(0)')
               .Define('EVT_thrust_x',    'EVT_thrust.at(1)')
               .Define('EVT_thrust_x_err','EVT_thrust.at(2)')
               .Define('EVT_thrust_y',    'EVT_thrust.at(3)')
               .Define('EVT_thrust_y_err','EVT_thrust.at(4)')
               .Define('EVT_thrust_z',    'EVT_thrust.at(5)')
               .Define('EVT_thrust_z_err','EVT_thrust.at(6)')              
               
#               .Define('EVT_thrust_costheta',  'EVT_thrust.at(7)')



               #.Define("EVT_thrustNP",         'Algorithms::minimize_thrust("Minuit2","Migrad")(RP_px_all_objects, RP_py_all_objects, RP_pz_all_objects)')
               #.Define("RP_thrustangleNP",     'Algorithms::getAxisCosTheta(EVT_thrustNP, RP_px_all_objects, RP_py_all_objects, RP_pz_all_objects)')
               #.Define("EVT_thrust",           'Algorithms::getThrustPointing(1)(RP_thrustangleNP, RP_e_all_objects, EVT_thrustNP)')

               #.Define('EVT_thrust_costheta',  'EVT_thrust.at(1)')
               

               .Define('EVT_sphericity',      'Algorithms::minimize_sphericity("Minuit2","Migrad")(RP_px_all_objects, RP_py_all_objects, RP_pz_all_objects)')
               .Define('EVT_sphericity_val',  'EVT_sphericity.at(0)')
               .Define('EVT_sphericity_x',    'EVT_sphericity.at(1)')
               .Define('EVT_sphericity_x_err','EVT_sphericity.at(2)')
               .Define('EVT_sphericity_y',    'EVT_sphericity.at(3)')
               .Define('EVT_sphericity_y_err','EVT_sphericity.at(4)')
               .Define('EVT_sphericity_z',    'EVT_sphericity.at(5)')
               .Define('EVT_sphericity_z_err','EVT_sphericity.at(6)')
               

               .Define('RP_thrustangle', 'Algorithms::getAxisCosTheta(EVT_thrust, RP_px_all_objects, RP_py_all_objects, RP_pz_all_objects)')
               .Define('RP_sphericityangle', 'Algorithms::getAxisCosTheta(EVT_sphericity, RP_px_all_objects, RP_py_all_objects, RP_pz_all_objects)') 
               


#################--------------------------------------------------------------------------

        )
        return df2

    #__________________________________________________________
    #Mandatory: output function, please make sure you return the branchlist as a python list
    def output():
        branchList = [
                "leptons", "n_leptons", "lepton_theta", "lepton_e", "lepton_px", 
                "lepton_py", "lepton_pz", "lepton_pt", "lepton_eta", "lepton_phi",
                "lepton_mass", "lepton_charge",
                "tw_jetconstituents_ee_genkt_ES", 
                "tw_jets_ee_genkt_ES_e",
                "tw_jets_ee_genkt_ES_px",
                "tw_jets_ee_genkt_ES_py",
                "tw_jets_ee_genkt_ES_pz",
                "tw_jets_ee_genkt_ES_theta",
                "tw_jets_ee_genkt_ES_charge",
                "tw_jets_ee_genkt_ES_flavour",
                "tw_jets_ee_genkt_ES_btag",
                "tw_jets_ee_genkt_ES_btag_true",
                "tw_jets_ee_genkt_ES_ctag",
                "RP_thrustangle", 
                "RP_sphericityangle",
                "EVT_thrust_x", "EVT_thrust_y", "EVT_thrust_z", "EVT_thrust_val",
                "EVT_sphericity_x", "EVT_sphericity_y", "EVT_sphericity_z", "EVT_sphericity_val"
#               "EVT_thrust_costheta"
        ]
        return branchList
