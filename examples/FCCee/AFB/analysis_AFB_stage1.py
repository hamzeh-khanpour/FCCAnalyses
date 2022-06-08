#Mandatory: List of processes
processList = {
    'p8_ee_Zbb_ecm91':{'fraction':0.01},#Run the full statistics in one output file named <outputDir>/p8_ee_Zbb_ecm91.root
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
               .Alias("Muon0", "Muon#0.index")
               .Define("electrons", "ReconstructedParticle::get(Electron0, ReconstructedParticles)")
               .Define("muons", "ReconstructedParticle::get(Muon0, ReconstructedParticles)")

               .Define("n_electrons", "ReconstructedParticle::get_n(electrons)")
               .Define("n_muons", "ReconstructedParticle::get_n(muons)")

               .Define("leptons", "ReconstructedParticle::merge(electrons, muons)")
               .Define("n_leptons", "ReconstructedParticle::get_n(leptons)")
               .Define("lepton_e", "ReconstructedParticle::get_e(leptons)")
               .Define("lepton_px", "ReconstructedParticle::get_px(leptons)")
               .Define("lepton_py", "ReconstructedParticle::get_py(leptons)")
               .Define("lepton_pz", "ReconstructedParticle::get_pz(leptons)")
               .Define("lepton_pt", "ReconstructedParticle::get_pt(leptons)")
               .Define("lepton_eta", "ReconstructedParticle::get_eta(leptons)")
               .Define("lepton_phi", "ReconstructedParticle::get_phi(leptons)")
               .Define("lepton_mass", "ReconstructedParticle::get_mass(leptons)")
               .Define("lepton_charge", "ReconstructedParticle::get_charge(leptons)")               
               
               .Define("electron_e", "ReconstructedParticle::get_e(electrons)")
               .Define("electron_px", "ReconstructedParticle::get_px(electrons)")
               .Define("electron_py", "ReconstructedParticle::get_py(electrons)")
               .Define("electron_pz", "ReconstructedParticle::get_pz(electrons)")
               .Define("electron_pt", "ReconstructedParticle::get_pt(electrons)")
               .Define("electron_eta", "ReconstructedParticle::get_eta(electrons)")
               .Define("electron_phi", "ReconstructedParticle::get_phi(electrons)")
               .Define("electron_mass", "ReconstructedParticle::get_mass(electrons)")
               .Define("electron_charge", "ReconstructedParticle::get_charge(electrons)")               
               
               .Define("muon_e", "ReconstructedParticle::get_e(muons)")
               .Define("muon_px", "ReconstructedParticle::get_px(muons)")
               .Define("muon_py", "ReconstructedParticle::get_py(muons)")
               .Define("muon_pz", "ReconstructedParticle::get_pz(muons)")
               .Define("muon_pt", "ReconstructedParticle::get_pt(muons)")
               .Define("muon_eta", "ReconstructedParticle::get_eta(muons)")
               .Define("muon_phi", "ReconstructedParticle::get_phi(muons)")
               .Define("muon_mass", "ReconstructedParticle::get_mass(muons)")
               .Define("muon_charge", "ReconstructedParticle::get_charge(muons)")               
        )
        return df2

    #__________________________________________________________
    #Mandatory: output function, please make sure you return the branchlist as a python list
    def output():
        branchList = [
                "leptons", "lepton_e", "lepton_pt", "lepton_mass", "lepton_charge", "lepton_eta", "lepton_phi",  
                "electrons", "electron_e", "electron_pt", "electron_mass", "electron_charge", "electron_eta", "electron_phi",               
                "muons", "muon_e", "muon_pt", "muon_mass", "muon_charge", "muon_eta", "muon_phi"
        ]
        return branchList
