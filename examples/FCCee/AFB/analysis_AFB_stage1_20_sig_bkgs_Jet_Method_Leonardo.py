# Forward-Backward-Asymmetry June 2022 - sig + bkgs

#Mandatory: List of processes
processList = {
'p8_ee_Zbb_ecm91':{'fraction':0.1,  'chunks':8}, #Run the full statistics in one output file named <outputDir>/p8_ee_Zbb_ecm91.root
'p8_ee_Zcc_ecm91':{'fraction':0.1,  'chunks':8}, #Run 50% of the statistics in two files named <outputDir>/p8_ee_WW_ecm240/chunk<N>.root
'p8_ee_Zuds_ecm91':{'fraction':0.1, 'chunks':8}, #Run 20% of the statistics in one file named <outputDir>/p8_ee_ZH_ecm240_out.root (example on how to change the output name)
'p8_ee_Zmumu_ecm91':{'fraction':1,  'chunks':1}, #Run 20% of the statistics in one file named <outputDir>/p8_ee_ZH_ecm240_out.root (example on how to change the output name)
}

#Mandatory: Production tag when running over EDM4Hep centrally produced events, this points to the yaml files for getting sample statistics
prodTag     = "FCCee/spring2021/IDEA/" 

#Optional: output directory, default is local running directory 
outputDir   = "/eos/experiment/fcc/ee/tmp/Hamzeh-Khanpour-FCNC-FCCee/AFB_FCCee_Udine_ICTP_20_sig_bkgs_eeKT/" 

#Optional: analysisName, default is ""
#analysisName = "My Analysis"

#Optional: ncpus, default is 4
nCPUS       = 8

#Optional running on HTCondor, default is False
runBatch    = False
#runBatch    = True

#Optional batch queue name when running on HTCondor, default is workday
#batchQueue = "longlunch"

#Optional computing account when running on HTCondor, default is group_u_FCC.local_gen
#compGroup = "group_u_FCC.local_gen"
#compGroup = "group_u_THEORY.u_t3" 

#Optional test file
#testFile ="root://eospublic.cern.ch//eos/experiment/fcc/ee/generation/DelphesEvents/spring2021/IDEA/p8_ee_ZH_ecm240/events_101027117.root"

#Mandatory: RDFanalysis class where the use defines the operations on the TTree
class RDFanalysis():

    #__________________________________________________________
    #Mandatory: analysers funtion to define the analysers to process, please make sure you return the last dataframe, in this example it is df2
    def analysers(df):
        df2 = (
            df

#################--------------------------------------------------------------------------


               .Alias("MCRecoAssociations0", "MCRecoAssociations#0.index")
               .Alias("MCRecoAssociations1", "MCRecoAssociations#1.index")

               .Alias("Particle0", "Particle#0.index")
               .Alias("Particle1", "Particle#1.index")
               
               .Define("MC_pdg",        "MCParticle::get_pdg(Particle)")
               .Define("MC_size",       "int(MC_pdg.size())")
               .Define("MC_px",         "MCParticle::get_px(Particle)")
               .Define("MC_py",         "MCParticle::get_py(Particle)")
               .Define("MC_pz",         "MCParticle::get_pz(Particle)")
               .Define("MC_p",          "MCParticle::get_p(Particle)")
               .Define("MC_e",          "MCParticle::get_e(Particle)")
               .Define("MC_theta",      "MCParticle::get_theta(Particle)")
               .Define("MC_charge",     "MCParticle::get_charge(Particle)")               
               .Define("MC_mass",       "MCParticle::get_mass(Particle)")
               .Define("MC_status",     "MCParticle::get_genStatus(Particle)")
               .Define("MC_vertex_x",   "MCParticle::get_vertex_x(Particle)")
               .Define("MC_vertex_y",   "MCParticle::get_vertex_y(Particle)")
               .Define("MC_vertex_z",   "MCParticle::get_vertex_z(Particle)")
#               .Define("MC_vertex_z",   "FCCAnalyses::MCParticle::get_vertex_z(Particle)")
               
  

               .Define("MC_PDG_ID", "FCCAnalyses::MCParticle::get_pdg(Particle)")
               .Define("MC_n",   "int(MC_PDG_ID.size())")
               .Define("MC_M1",  "myUtils::get_MCMother1(Particle,Particle0)")
               .Define("MC_M2",  "myUtils::get_MCMother2(Particle,Particle0)")
               .Define("MC_D1",  "myUtils::get_MCDaughter1(Particle,Particle1)")
               .Define("MC_D2",  "myUtils::get_MCDaughter2(Particle,Particle1)")



#################--------------------------------------------------------------------------



                .Define("selected_jets", "ReconstructedParticle::sel_p(0.0)(Jet)")           # select only jets with a pT > 0 GeV   
                .Define("n_jets",        "ReconstructedParticle::get_n(selected_jets)")      # count how many jets are in the event in total
                .Define("seljet_p",      "ReconstructedParticle::get_p(selected_jets)")      # momentum p
                .Define("seljet_eta",    "ReconstructedParticle::get_eta(selected_jets)")    # create branch with jet eta
		        .Define("seljet_theta",  "ReconstructedParticle::get_theta(selected_jets)")  # theta
		        .Define("seljet_charge", "ReconstructedParticle::get_charge(selected_jets)") # charge
                .Define("seljet_e",      "ReconstructedParticle::get_e(selected_jets)")      # Energy of jet
                .Define("seljet_px",     "ReconstructedParticle::get_px(selected_jets)")     # create branch with jet px
                .Define("seljet_py",     "ReconstructedParticle::get_py(selected_jets)")     # create branch with jet py
                .Define("seljet_pz",     "ReconstructedParticle::get_pz(selected_jets)")     # create branch with jet pz 


               .Alias("Jet3","Jet#3.index") 
               .Define("JET_btag",  "ReconstructedParticle::getJet_btag(Jet3, ParticleIDs, ParticleIDs_0)")
               .Define("EVT_nbtag", "ReconstructedParticle::getJet_ntags(JET_btag)")

               .Define("RP_px",          "ReconstructedParticle::get_px(ReconstructedParticles)")
               .Define("RP_py",          "ReconstructedParticle::get_py(ReconstructedParticles)")
               .Define("RP_pz",          "ReconstructedParticle::get_pz(ReconstructedParticles)")               
               .Define("RP_mass",        "ReconstructedParticle::get_mass(ReconstructedParticles)")
               .Define("RP_e",           "ReconstructedParticle::get_e(ReconstructedParticles)")               
               .Define("RP_charge",      "ReconstructedParticle::get_charge(ReconstructedParticles)")

             
               .Define("pseudo_jets",   "JetClusteringUtils::set_pseudoJets_xyzm(RP_px, RP_py, RP_pz, RP_mass)")

#################--------------------------------------------------------------------------


              
               #====================================================================================================================
                    # re-clustering the jet, eeKT Durham (clustering_eeKT_E Scheme)-Exclusive jet selection exactly 2 jets
               #======================================================================================================================              
              
               .Define("eeKT_FCCAnalysesJets_ee_ES","JetClustering::clustering_ee_kt(2, 2, 1, 0)(pseudo_jets)")     # eeKT
                 
               .Define("eeKT_jets_ee_ES",           "JetClusteringUtils::get_pseudoJets(eeKT_FCCAnalysesJets_ee_ES)")   #get the jets constituents out of the struct
               .Define("eeKT_jetconstituents_ee_ES","JetClusteringUtils::get_constituents(eeKT_FCCAnalysesJets_ee_ES)")
               .Define("eeKT_jets_ee_ES_e",         "JetClusteringUtils::get_e(eeKT_jets_ee_ES)")
               .Define("eeKT_jets_ee_ES_px",        "JetClusteringUtils::get_px(eeKT_jets_ee_ES)")
               .Define("eeKT_jets_ee_ES_py",        "JetClusteringUtils::get_py(eeKT_jets_ee_ES)")
               .Define("eeKT_jets_ee_ES_pz",        "JetClusteringUtils::get_pz(eeKT_jets_ee_ES)")
               .Define("eeKT_jets_ee_ES_theta",     "JetClusteringUtils::get_theta(eeKT_jets_ee_ES)")         
               .Define("eeKT_jets_ee_ES_flavour",   "JetTaggingUtils::get_flavour(eeKT_jets_ee_ES, Particle)")
               .Define("eeKT_jets_ee_ES_btag_true", "JetTaggingUtils::get_btag(eeKT_jets_ee_ES_flavour, 1.0)")
               .Define("eeKT_jets_ee_ES_btag",      "JetTaggingUtils::get_btag(eeKT_jets_ee_ES_flavour, 0.80)")
               .Define("eeKT_jets_ee_ES_ctag_true", "JetTaggingUtils::get_ctag(eeKT_jets_ee_ES_flavour, 1.0)")
               .Define("eeKT_jets_ee_ES_ctag",      "JetTaggingUtils::get_ctag(eeKT_jets_ee_ES_flavour, 0.10)")
              

#################--------------------------------------------------------------------------


               #############################################
               ##              Build the thrust           ##
               #############################################


               .Define("RP_px_all_objects",         "ReconstructedParticle::get_px(ReconstructedParticles)")
               .Define("RP_py_all_objects",         "ReconstructedParticle::get_py(ReconstructedParticles)")
               .Define("RP_pz_all_objects",         "ReconstructedParticle::get_pz(ReconstructedParticles)")
               .Define("RP_e_all_objects",          "ReconstructedParticle::get_e(ReconstructedParticles)")


               .Define('EVT_thrustNP',      'Algorithms::minimize_thrust("Minuit2","Migrad")(RP_px_all_objects, RP_py_all_objects, RP_pz_all_objects)')
               .Define('RP_thrustangleNP',  'Algorithms::getAxisCosTheta(EVT_thrustNP, RP_px_all_objects, RP_py_all_objects, RP_pz_all_objects)')               
               .Define("EVT_thrust",        'Algorithms::getThrustPointing(1)(RP_thrustangleNP, RP_e_all_objects, EVT_thrustNP)')
               .Define('RP_thrustangle',    'Algorithms::getAxisCosTheta(EVT_thrust, RP_px_all_objects, RP_py_all_objects, RP_pz_all_objects)')

               
               .Define('EVT_thrust_Mag',  'EVT_thrust.at(0)')
               .Define('EVT_thrust_x',    'EVT_thrust.at(1)')
               .Define('EVT_thrust_x_err','EVT_thrust.at(2)')
               .Define('EVT_thrust_y',    'EVT_thrust.at(3)')
               .Define('EVT_thrust_y_err','EVT_thrust.at(4)')  
               .Define('EVT_thrust_z',    'EVT_thrust.at(5)')
               .Define('EVT_thrust_z_err','EVT_thrust.at(6)')


#################--------------------------------------------------------------------------


               #############################################
               ##              Build the sphericity       ##
               #############################################
               
               
               .Define('EVT_sphericity',      'Algorithms::minimize_sphericity("Minuit2","Migrad")(RP_px_all_objects, RP_py_all_objects, RP_pz_all_objects)')
               .Define('EVT_sphericity_Mag',  'EVT_sphericity.at(0)')
               .Define('EVT_sphericity_x',    'EVT_sphericity.at(1)')
               .Define('EVT_sphericity_x_err','EVT_sphericity.at(2)')
               .Define('EVT_sphericity_y',    'EVT_sphericity.at(3)')
               .Define('EVT_sphericity_y_err','EVT_sphericity.at(4)')
               .Define('EVT_sphericity_z',    'EVT_sphericity.at(5)')
               .Define('EVT_sphericity_z_err','EVT_sphericity.at(6)')
               
               .Define('RP_sphericityangle',  'Algorithms::getAxisCosTheta(EVT_sphericity, RP_px_all_objects, RP_py_all_objects, RP_pz_all_objects)') 
               


               #====================================================================================================================
                    # EFlowTrack
               #======================================================================================================================              
     

               #############################################
               ##               Build MC Vertex           ##
               #############################################
               .Define("MCVertexObject", "myUtils::get_MCVertexObject(Particle, Particle0)")
               .Define("MC_Vertex_x",    "myUtils::get_MCVertex_x(MCVertexObject)")
               .Define("MC_Vertex_y",    "myUtils::get_MCVertex_y(MCVertexObject)")
               .Define("MC_Vertex_z",    "myUtils::get_MCVertex_z(MCVertexObject)")
               .Define("MC_Vertex_ind",  "myUtils::get_MCindMCVertex(MCVertexObject)")
               .Define("MC_Vertex_ntrk", "myUtils::get_NTracksMCVertex(MCVertexObject)")
               .Define("MC_Vertex_n",    "int(MC_Vertex_x.size())")
               .Define("MC_Vertex_PDG",  "myUtils::get_MCpdgMCVertex(MCVertexObject, Particle)")
               .Define("MC_Vertex_PDGmother",  "myUtils::get_MCpdgMotherMCVertex(MCVertexObject, Particle)")
               .Define("MC_Vertex_PDGgmother", "myUtils::get_MCpdgGMotherMCVertex(MCVertexObject, Particle)")


               #############################################
               ##              Build Reco Vertex          ##
               #############################################
               .Define("VertexObject", "myUtils::get_VertexObject(MCVertexObject,ReconstructedParticles,EFlowTrack_1,MCRecoAssociations0,MCRecoAssociations1)")


               #############################################
               ##          Build PV var and filter        ##
               #############################################
               .Define("EVT_hasPV",    "myUtils::hasPV(VertexObject)")
               .Define("EVT_NtracksPV", "float(myUtils::get_PV_ntracks(VertexObject))")
               .Define("EVT_NVertex",   "float(VertexObject.size())")
               .Filter("EVT_hasPV==1")


               #############################################
               ##          Build RECO P with PID          ##
               #############################################
               .Define("RecoPartPID" ,"myUtils::PID(ReconstructedParticles, MCRecoAssociations0,MCRecoAssociations1,Particle)")

               #############################################
               ##    Build RECO P with PID at vertex      ##
               #############################################
               .Define("RecoPartPIDAtVertex" ,"myUtils::get_RP_atVertex(RecoPartPID, VertexObject)")

               #############################################
               ##         Build vertex variables          ##
               #############################################
               .Define("Vertex_x",        "myUtils::get_Vertex_x(VertexObject)")
               .Define("Vertex_y",        "myUtils::get_Vertex_y(VertexObject)")
               .Define("Vertex_z",        "myUtils::get_Vertex_z(VertexObject)")
               .Define("Vertex_xErr",     "myUtils::get_Vertex_xErr(VertexObject)")
               .Define("Vertex_yErr",     "myUtils::get_Vertex_yErr(VertexObject)")
               .Define("Vertex_zErr",     "myUtils::get_Vertex_zErr(VertexObject)")

               .Define("Vertex_chi2",     "myUtils::get_Vertex_chi2(VertexObject)")
               .Define("Vertex_mcind",    "myUtils::get_Vertex_indMC(VertexObject)")
               .Define("Vertex_ind",      "myUtils::get_Vertex_ind(VertexObject)")
               .Define("Vertex_isPV",     "myUtils::get_Vertex_isPV(VertexObject)")
               .Define("Vertex_ntrk",     "myUtils::get_Vertex_ntracks(VertexObject)")
               .Define("Vertex_n",        "int(Vertex_x.size())")
               .Define("Vertex_mass",     "myUtils::get_Vertex_mass(VertexObject,RecoPartPIDAtVertex)")

     

#################--------------------------------------------------------------------------

        )
        return df2

    #__________________________________________________________
    #Mandatory: output function, please make sure you return the branchlist as a python list
    def output():
        branchList = [
                "MC_pdg", "MC_size", "MC_px", "MC_py", "MC_pz", "MC_p", "MC_e", 
                "MC_theta", "MC_charge", "MC_mass", "MC_status", 
                "MC_PDG_ID", "MC_n", "MC_M1", "MC_M2", "MC_D1", "MC_D2",
                "n_jets", "seljet_eta", "seljet_p", "seljet_theta", "seljet_charge", "seljet_e", 
                "seljet_px", "seljet_py", "seljet_pz",  
                "eeKT_jetconstituents_ee_ES", 
                "eeKT_jets_ee_ES_e",
                "eeKT_jets_ee_ES_px",
                "eeKT_jets_ee_ES_py",
                "eeKT_jets_ee_ES_pz",
                "eeKT_jets_ee_ES_theta",
                "eeKT_jets_ee_ES_flavour",
                "eeKT_jets_ee_ES_btag_true",
                "eeKT_jets_ee_ES_btag",
                "eeKT_jets_ee_ES_ctag_true",
                "eeKT_jets_ee_ES_ctag", 
                "MC_Vertex_x", "MC_Vertex_y", "MC_Vertex_z",
                "MC_Vertex_ntrk", "MC_Vertex_n",
                "MC_Vertex_PDG","MC_Vertex_PDGmother","MC_Vertex_PDGgmother",
                "Vertex_x", "Vertex_y", "Vertex_z",
                "Vertex_xErr", "Vertex_yErr", "Vertex_zErr",
                "Vertex_isPV", "Vertex_ntrk", "Vertex_chi2", "Vertex_n",
                "Vertex_mass"
        ]
        return branchList
