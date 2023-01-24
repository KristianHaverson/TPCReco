#include <cstdlib>
#include <iostream>
#include <vector>

#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "TString.h"
#include "TStopwatch.h"

#include <boost/program_options.hpp>

#include "IonRangeCalculator.h"
#include "dEdxFitter.h"
#include "TrackBuilder.h"
#include "EventSourceMC.h"

#include "RecoOutput.h"
#include "RunIdParser.h"
#include "InputFileHelper.h"
#include "MakeUniqueName.h"
#include "colorText.h"

#include "EventTPC.h"
/////////////////////////////////////
/////////////////////////////////////
std::string createROOTFileName(const  std::string & grawFileName){

  std::string rootFileName = grawFileName;
  std::string::size_type index = rootFileName.find(",");
  if(index!=std::string::npos){
    rootFileName = grawFileName.substr(0,index);
  }
  index = rootFileName.rfind("/");
  if(index!=std::string::npos){
    rootFileName = rootFileName.substr(index+1,-1);
  }
  if(rootFileName.find("CoBo_ALL_AsAd_ALL")!=std::string::npos){
    rootFileName = rootFileName.replace(0,std::string("CoBo_ALL_AsAd_ALL").size(),"TrackTree");
  }
  else if(rootFileName.find("CoBo0_AsAd")!=std::string::npos){
    rootFileName = rootFileName.replace(0,std::string("CoBo0_AsAd").size()+1,"TrackTree");
  }
  else if(rootFileName.find("EventTPC")!=std::string::npos){
    rootFileName = rootFileName.replace(0,std::string("EventTPC").size(),"TrackTree");
  }
  else{
    std::cout<<KRED<<"File format unknown: "<<RST<<rootFileName<<std::endl;
    exit(1);
  }
  index = rootFileName.rfind("graw");
  if(index!=std::string::npos){
    rootFileName = rootFileName.replace(index,-1,"root");
  }
  
  return rootFileName;
}
/////////////////////////////////////
/////////////////////////////////////
int makeTrackTree(const  std::string & geometryFileName,
		  const  std::string & dataFileName);
/////////////////////////////////////
/////////////////////////////////////
boost::program_options::variables_map parseCmdLineArgs(int argc, char **argv){

  boost::program_options::options_description cmdLineOptDesc("Allowed options");
  cmdLineOptDesc.add_options()
    ("help", "produce help message")
    ("geometryFile",  boost::program_options::value<std::string>()->required(), "string - path to the geometry file.")
    ("dataFile",  boost::program_options::value<std::string>()->required(), "string - path to data file.");
  
  boost::program_options::variables_map varMap;        
try {     
    boost::program_options::store(boost::program_options::parse_command_line(argc, argv, cmdLineOptDesc), varMap);
    if (varMap.count("help")) {
      std::cout << "makeTrackTree" << "\n\n";
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

  TStopwatch aStopwatch;
  aStopwatch.Start();

  std::string geometryFileName, dataFileName;
  boost::program_options::variables_map varMap = parseCmdLineArgs(argc, argv);
  boost::property_tree::ptree tree;
  if(argc<3){
    char text[] = "--help";
    char *argvTmp[] = {text, text};
    parseCmdLineArgs(2,argvTmp);
    return 1;
  }
  if (varMap.count("geometryFile")) {
    geometryFileName = varMap["geometryFile"].as<std::string>();
  }
  if (varMap.count("dataFile")) {
    dataFileName = varMap["dataFile"].as<std::string>();
  }

  int nEntriesProcessed = 0;
  if(dataFileName.size() && geometryFileName.size()){
    nEntriesProcessed = makeTrackTree(geometryFileName, dataFileName);
  }
  else{
    std::cout<<KRED<<"Configuration not complete: "<<RST
	     <<" geometryFile: "<<geometryFileName<<"\n"
	     <<" dataFile: "<<dataFileName
	     <<std::endl;
  }

  aStopwatch.Stop();
  std::cout<<KBLU<<"Real time:       "<<RST<<aStopwatch.RealTime()<<" s"<<std::endl;
  std::cout<<KBLU<<"CPU time:        "<<RST<<aStopwatch.CpuTime()<<" s"<<std::endl;
  std::cout<<KBLU<<"Processing rate: "<<RST<<nEntriesProcessed/aStopwatch.RealTime()<< " ev/s"<<std::endl;

  return 0;
}
/////////////////////////////
////////////////////////////
// Define some simple structures
typedef struct {Float_t eventId, frameId,
    eventTypeGen,
    alphaRangeGen,
    alphaEnergyGen,
    chargeGen,
    cosThetaGen, phiGen,
    ///
    eventTypeReco,
    alphaRangeReco,
    alphaEnergyReco,
    chargeReco,
    cosThetaReco, phiReco;
    } TrackData;
/////////////////////////
int makeTrackTree(const  std::string & geometryFileName, const  std::string & dataFileName) {

  // ** EVENT SOURCE ** //
  std::shared_ptr<EventSourceMC> myEventSource = std::make_shared<EventSourceMC>(geometryFileName);

  // ** TTREE SETUP ** //
  std::string rootFileName = createROOTFileName(dataFileName);
  TFile outputROOTFile(rootFileName.c_str(),"RECREATE");
  TTree *tree = new TTree("trackTree", "Track tree");
  TrackData track_data;
  std::string leafNames = "";
  leafNames += "eventId:frameId:";
  leafNames += "eventTypeGen:alphaRangeGen:alphaEnergyGen:chargeGen:cosThetaGen:phiGen:";  
  leafNames += "eventTypeReco:alphaRangeReco:alphaEnergyReco:chargeReco:cosThetaReco:phiReco";
  tree->Branch("track",&track_data,leafNames.c_str());

  // ** GEOMETRY ** //
  int index = geometryFileName.find("mbar");
  double pressure = stof(geometryFileName.substr(index-3, 3));
  TrackBuilder myTkBuilder;
  myTkBuilder.setGeometry(myEventSource->getGeometry());
  myTkBuilder.setPressure(pressure);
  IonRangeCalculator myRangeCalculator(gas_mixture_type::CO2,pressure,293.15);

  // ** OUTFILE ** //
  std::string fileName = InputFileHelper::tokenize(dataFileName)[0];
  std::size_t last_dot_position = fileName.find_last_of(".");
  std::size_t last_slash_position = fileName.find_last_of("//");
  std::string recoFileName = MakeUniqueName("Reco_"+fileName.substr(last_slash_position+1,last_dot_position-last_slash_position-1)+".root");
  std::shared_ptr<eventraw::EventInfo> myEventInfo = std::make_shared<eventraw::EventInfo>();
  RecoOutput myRecoOutput.open(recoFileName);
  
  // ** LOAD DATA ** //
  myEventSource->loadDataFile(dataFileName);
  std::cout<<KBLU<<"File with "<<RST<<myEventSource->numberOfEntries()<<" frames loaded."<<std::endl;

  // ** MAIN LOOP ** //
  unsigned int nEntries = myEventSource->numberOfEntries();
  nEntries = 10; //TEST

  for(unsigned int iEntry=0;iEntry<nEntries;++iEntry){

    // ** COUNTER **//
    std::cout<<KBLU<<"Processed: "<<iEntry<<"/"<<nEntries<<RST<<std::endl;
    
    // ** RECONSTRUCT ** //
    myEventSource->loadFileEntry(iEntry);    
    *myEventInfo = myEventSource->getCurrentEvent()->GetEventInfo();    
    myTkBuilder.setEvent(myEventSource->getCurrentEvent());
    myTkBuilder.reconstruct();

    // ** EXTRACT n FILL ** //
    int eventId = myEventSource->getCurrentEvent()->GetEventInfo().GetEventId();
    const Track3D & aTrack3DGen = myEventSource->getGeneratedTrack();
    const Track3D & aTrack3DReco = myTkBuilder.getTrack3D(0);

    track_data.frameId = iEntry;
    track_data.eventId = eventId;

    track_data.eventTypeGen = aTrack3DGen.getSegments().front().getPID() + aTrack3DGen.getSegments().back().getPID();    
    track_data.alphaRangeGen =  aTrack3DGen.getSegments().front().getLength();    
    track_data.alphaEnergyGen = track_data.alphaRangeGen>0 ? 1E3*myRangeCalculator.getIonEnergyMeV(pid_type::ALPHA, track_data.alphaRangeGen):0.0;
    track_data.chargeGen = 100*track_data.alphaEnergyGen;//aTrack3DGen.getIntegratedCharge(track_data.alphaRangeGen);
    const TVector3 & tangentGen = aTrack3DGen.getSegments().front().getTangent();
    track_data.cosThetaGen = -tangentGen.X();
    track_data.phiGen = atan2(-tangentGen.Z(), tangentGen.Y());

    track_data.eventTypeReco = aTrack3DReco.getSegments().front().getPID() + +aTrack3DReco.getSegments().back().getPID();    
    track_data.alphaRangeReco =  aTrack3DReco.getSegments().front().getLength();    
    track_data.alphaEnergyReco = track_data.alphaRangeReco>0 ? myRangeCalculator.getIonEnergyMeV(pid_type::ALPHA, track_data.alphaRangeReco):0.0;
    track_data.chargeReco = aTrack3DReco.getIntegratedCharge(track_data.alphaRangeReco);
    const TVector3 & tangentReco = aTrack3DReco.getSegments().front().getTangent();
    track_data.cosThetaReco = -tangentReco.X();
    track_data.phiReco = atan2(-tangentReco.Z(), tangentReco.Y());

    
    tree->Fill();    
  }
  outputROOTFile.Write();
  return nEntries;
}
/////////////////////////////
////////////////////////////

