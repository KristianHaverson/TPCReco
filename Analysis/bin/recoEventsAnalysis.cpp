#include <cstdlib>
#include <iostream>
#include <vector>

#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "TString.h"
#include "TTreeIndex.h"

#include <boost/program_options.hpp>

#include "GeometryTPC.h"
#include "Track3D.h"
#include "HIGGS_analysis.h"
#include "HIGS_trees_analysis.h"
#include "HIGS_trees_analysis_kris.h"

#include "colorText.h"

int analyzeRecoEvents(const  std::string & geometryFileName, 
		      const  std::string & dataFileName,
		      const  float & beamEnergy, // [MeV]
		      const  TVector3 & beamDir, // unit vector in LAB detector frame
		      const  double & pressure, // [mbar]
		      const  double & temperature, // [K]
		      const  bool & makeTreeFlag);

enum class BeamDirection{
  X,
  MINUS_X
};

std::istream& operator>>(std::istream& in, BeamDirection& direction){
  std::string token;
  in >> token;
  if (token == "x" || token == "X"){
     direction = BeamDirection::X;
  }
  else if (token == "-x" || token == "-X"){
    direction = BeamDirection::MINUS_X;
  } 
  else
    in.setstate(std::ios_base::failbit);
  return in;
  }

boost::program_options::variables_map parseCmdLineArgs(int argc, char **argv){

  //  bool optionTree=true;
  boost::program_options::options_description cmdLineOptDesc("Allowed options");
  cmdLineOptDesc.add_options()
    ("help", "produce help message")
    ("geometryFile",  boost::program_options::value<std::string>()->required(), "string - path to the geometry file")
    ("dataFile",  boost::program_options::value<std::string>()->required(), "string - path to data file")
    ("beamEnergy", boost::program_options::value<float>()->required(), "float - LAB gamma beam energy [MeV]")
    ("beamDir", boost::program_options::value<BeamDirection>()->required(), "string - LAB gamma beam direction [\"x\" xor \"-x\"]")
    ("pressure", boost::program_options::value<float>()->required(), "float - CO2 pressure [mbar]")
    ("temperature", boost::program_options::value<float>()->default_value(273.15+20), "(option) float - CO2 temperature [K], default=293.15K (20 C)")
    ("noTree", boost::program_options::bool_switch()->default_value(false), "(option) bool - flag to skip creating additional TTree for 1,2,3-prongs, default=false ");
  
  boost::program_options::variables_map varMap;  

  try {     
    boost::program_options::store(boost::program_options::parse_command_line(argc, argv, cmdLineOptDesc), varMap);
    if (varMap.count("help")) {
      std::cout << "recoEventAnalysis" << "\n\n";
      std::cout << cmdLineOptDesc << std::endl;
      exit(1);
    }
    boost::program_options::notify(varMap);
  } catch (const std::exception &e) {
    std::cerr << e.what() << '\n';
    std::cout << cmdLineOptDesc << std::endl;
    exit(1);
  }

  return varMap;
}
/////////////////////////////////////
/////////////////////////////////////
int main(int argc, char **argv){

  auto varMap = parseCmdLineArgs(argc, argv);
  auto geometryFileName = varMap["geometryFile"].as<std::string>();
  auto dataFileName = varMap["dataFile"].as<std::string>();
  auto beamEnergy = varMap["beamEnergy"].as<float>();
  auto pressure = varMap["pressure"].as<float>();
  auto temperature = varMap["temperature"].as<float>();
  auto makeTreeFlag = !varMap["noTree"].as<bool>();
  TVector3 beamDir;
  switch(varMap["beamDir"].as<BeamDirection>()){
    case BeamDirection::X : 
      beamDir = TVector3(1,0,0);
      break;
    case BeamDirection::MINUS_X : 
      beamDir = TVector3(-1,0,0);
      break;
    default:
      return 1;
  }
  analyzeRecoEvents(geometryFileName, dataFileName, beamEnergy, beamDir, pressure, temperature, makeTreeFlag);
  return 0;
}
/////////////////////////////
////////////////////////////
std::shared_ptr<GeometryTPC> loadGeometry(const std::string fileName){
  return std::make_shared<GeometryTPC>(fileName.c_str(), false);
}
////////////////////////////
////////////////////////////
int analyzeRecoEvents(const  std::string & geometryFileName,
		      const  std::string & dataFileName,
		      const  float & beamEnergy, // [MeV]
		      const  TVector3 & beamDir, // unit vector in LAB detector frame
		      const  double & pressure, // [mbar]
		      const  double & temperature, // [K]
		      const  bool & makeTreeFlag){

  std::cout << __FUNCTION__ << ": Input parameters:" << std::endl
	    << "* geometry file: " << geometryFileName << std::endl
	    << "* input ROOT file: " << dataFileName << std::endl
	    << "* LAB gamma beam energy: " << beamEnergy << " MeV" << std::endl
	    << "* LAB gamma beam direction in DET coordinates: [x=" << beamDir.X() << ", y=" << beamDir.Y() << ", z=" << beamDir.Z() << "]" << std::endl
	    << "* CO2 pressure (for ion ranges):  " << pressure <<" mbar" <<std::endl
	    << "* CO2 temperature (for ion ranges): " << temperature << " K, "<< (temperature-273.15) << " C" << std::endl;

  TFile *aFile = new TFile(dataFileName.c_str());
  if(!aFile || !aFile->IsOpen()){
    std::cout<<KRED<<"Input file: "<<RST<<dataFileName
	     <<KRED<<" not found!"<<RST
	     <<std::endl;
    return -1;
  }
  std::shared_ptr<GeometryTPC> aGeometry = loadGeometry(geometryFileName);
  if(!aGeometry){
    std::cout<<KRED<<"Geometry file: "<<RST<<geometryFileName
	     <<KRED<<" not found!"<<RST
	     <<std::endl;
    return -1;
  }
  if(beamEnergy<=0) {
    std::cout<<KRED<<"Wrong beam energy: "<<RST<<beamEnergy<<" MeV"
	     <<std::endl;
    return -1;
  }
  if(beamDir.Mag2()==0) {
    std::cout<<KRED<<"Wrong beam direction vector: "<<RST<<"["<<beamDir.X()<<", "<<beamDir.Y()<<", "<<beamDir.Z()<<"]"
	     <<std::endl;
    return -1;
  }
  if(pressure<50.0 || pressure>1100.0) {
    std::cout<<KRED<<"Wrong CO2 pressure: "<<RST<<pressure<<" mbar"
	     <<std::endl;
    return -1;
  }
  if(temperature<273.15 || temperature>273.15+50) {
    std::cout<<KRED<<"Wrong CO2 temperature: "<<RST<<temperature<<" K"
	     <<std::endl;
    return -1;
  }

  std::shared_ptr<HIGGS_analysis> myAnalysis(new HIGGS_analysis(aGeometry, beamEnergy, beamDir, pressure, temperature));
  std::shared_ptr<HIGS_trees_analysis> myTreesAnalysis(0);
  if(makeTreeFlag) myTreesAnalysis=std::make_shared<HIGS_trees_analysis>(aGeometry, beamEnergy, beamDir, pressure, temperature);
  //!
  //std::shared_ptr<HIGS_trees_analysis_kris> myTreesAnalysisK(0);
  //if(makeTreeFlag) myTreesAnalysisK=std::make_shared<HIGS_trees_analysis_kris>(aGeometry, beamEnergy, beamDir, pressure, temperature);
  //! 
  TTree *aTree = (TTree*)aFile->Get("TPCRecoData");
  if(!aTree) {
    std::cerr<<KRED<<"ERROR: Cannot find 'TPCRecoData' tree!"<<RST<<std::endl;
    return -1;
  }
  // NOTE: * branch RecoEvent - must always be present
  //       * branch EventInfo - is optional (eg. for Monte Carlo)
  Track3D *aTrack = new Track3D();
  TBranch *aBranch  = aTree->GetBranch("RecoEvent");
  if(!aBranch) {
    std::cerr<<KRED<<"ERROR: Cannot find 'RecoEvent' branch!"<<RST<<std::endl;
    return -1;
  }
  aBranch->SetAddress(&aTrack);
  
  eventraw::EventInfo *aEventInfo = 0;
  TBranch *aBranchInfo = aTree->GetBranch("EventInfo");
  if(!aBranchInfo) {
   std::cerr<<KRED<<"WARNING: "
	    <<"Cannot find 'EventInfo' branch!"<<RST<<std::endl;
  }
  else{
    aEventInfo = new eventraw::EventInfo();
    aBranchInfo->SetAddress(&aEventInfo);
  }

  unsigned int nEntries = aTree->GetEntries();
  
  // When (optional) "EventInfo" branch is present try to sort input tree in ascending order of {runID, eventID}
  if(!aBranchInfo) {
    
    std::cout << __FUNCTION__ << ": Starting to loop " << nEntries << " events without sorting" << std::endl;
    for(unsigned int iEntry=0;iEntry<nEntries;++iEntry){
      aBranch->GetEntry(iEntry); 
      for (auto & aSegment: aTrack->getSegments())  aSegment.setGeometry(aGeometry); // need TPC geometry for track projections
      myAnalysis->fillHistos(aTrack);
      if(makeTreeFlag) myTreesAnalysis->fillTrees(aTrack, aEventInfo);
  //    if(makeTreeFlag) myTreesAnalysisK->fillTrees(aTrack, aEventInfo);
    }
    
  } else {
    
    std::cout << __FUNCTION__ << ": Starting to loop " << nEntries << " events with sorting by {runId, eventId}" << std::endl;
    aTree->BuildIndex("runId", "eventId");
    TTreeIndex *I=(TTreeIndex*)aTree->GetTreeIndex(); // get the tree index
    Long64_t* index=I->GetIndex();

    for(unsigned int iEntry=0;iEntry<nEntries;++iEntry){
      aBranch->GetEntry(index[iEntry]);
      aBranchInfo->GetEntry(index[iEntry]);

      ////// DEBUG
      //      std::cout << __FUNCTION__ << ": run=" << aEventInfo->GetRunId() << ", event=" << aEventInfo->GetEventId() << std::endl;
      ////// DEBUG
      
      for (auto & aSegment: aTrack->getSegments())  aSegment.setGeometry(aGeometry); // need TPC geometry for track projections
      myAnalysis->fillHistos(aTrack);
      if(makeTreeFlag) myTreesAnalysis->fillTrees(aTrack, aEventInfo);
    //  if(makeTreeFlag) myTreesAnalysisK->fillTrees(aTrack, aEventInfo);
    }
    
  }
  /*
    for(unsigned int iEntry=0;iEntry<nEntries;++iEntry){
    aBranch->GetEntry(iEntry); 
    if(aBranchInfo) aBranchInfo->GetEntry(iEntry); // this branch is optional
    for (auto & aSegment: aTrack->getSegments())  aSegment.setGeometry(aGeometry); // need TPC geometry for track projections
    myAnalysis->fillHistos(aTrack);
    if(makeTreeFlag) myTreesAnalysis->fillTrees(aTrack, aEventInfo);
    }			      
  */
  return 0;
}
/////////////////////////////
////////////////////////////

