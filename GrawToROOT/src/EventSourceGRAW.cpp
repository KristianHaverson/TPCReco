#include <cstdlib>
#include <iostream>
#include <ostream>
#include <iomanip>
#include <iterator>

#include "EventSourceGRAW.h"
#include "colorText.h"


#include "get/graw2dataframe.h"
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
EventSourceGRAW::EventSourceGRAW(const std::string & geometryFileName) {

  loadGeometry(geometryFileName);
  GRAW_EVENT_FRAGMENTS = myGeometryPtr->GetAsadNboards();
  myPedestalCalculator.SetGeometryAndInitialize(myGeometryPtr);

  std::string formatsFilePath = "./CoboFormats.xcfg";
  myFrameLoader.initialize(formatsFilePath);
 
  //minSignalCell = 2;//FIXME read from config
 // maxSignalCell = 500;//FIXME read from config
}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
EventSourceGRAW::~EventSourceGRAW() { }
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
void EventSourceGRAW::setRemovePedestal(bool aFlag){
  removePedestal = aFlag;
}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
std::shared_ptr<EventTPC> EventSourceGRAW::getNextEvent(){

  unsigned long int currentEventId = myCurrentEvent->GetEventId();
  auto it = myFramesMap.find(currentEventId);
  unsigned int lastEventFrame = *it->second.rbegin();
  if(lastEventFrame<nEntries-1) ++lastEventFrame;
  loadFileEntry(lastEventFrame);
  return myCurrentEvent;  
}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
std::shared_ptr<EventTPC> EventSourceGRAW::getPreviousEvent(){

  unsigned int currentEventId = myCurrentEvent->GetEventId();
  auto it = myFramesMap.find(currentEventId);
  unsigned int startingEventIndexFrame = *it->second.begin();
  if(startingEventIndexFrame>0) --startingEventIndexFrame; 
  loadFileEntry(startingEventIndexFrame);
  return myCurrentEvent;  
}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
void EventSourceGRAW::loadDataFile(const std::string & fileName){

  EventSourceBase::loadDataFile(fileName);

  myFile =  std::make_shared<TGrawFile>(fileName.c_str());
  if(!myFile){
    std::cerr<<KRED<<"Can not open file: "<<fileName<<"!"<<RST<<std::endl;
    exit(1);
  }
  nEntries = myFile->GetGrawFramesNumber();
  
  const int firstEventSize=10;
  if(fileName!=myFilePath || nEntries<firstEventSize){
    findStartingIndex(firstEventSize);
  }

  myFilePath = fileName;
  myNextFilePath = getNextFilePath();
  myFramesMap.clear();
  myASADMap.clear();
  myReadEntriesSet.clear();
  isFullFileScanned = false;
  /*
  for(unsigned int iEntry=0;iEntry<nEntries;++iEntry){
    loadGrawFrame(iEntry, false);
    int currentEventId = myDataFrame.fHeader.fEventIdx;
    int ASAD_idx = myDataFrame.fHeader.fAsadIdx;
    std::cout<<"iEntry: "<<iEntry
	     <<" currentEventId: "<<currentEventId
	     <<" ASAD: "<<ASAD_idx
	     <<std::endl;
  }
  std::cout<<KBLU<<"End of file"<<RST<<std::endl;  
  */ 
}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
bool EventSourceGRAW::loadGrawFrame(unsigned int iEntry, bool readFullEvent){

  std::string tmpFilePath = myFilePath;
  if(iEntry>=nEntries){
    tmpFilePath = myNextFilePath;
    iEntry -= nEntries;
  }
  std::cout.setstate(std::ios_base::failbit);
  bool dataFrameRead = myFrameLoader.getGrawFrame(tmpFilePath, iEntry+1, myDataFrame, readFullEvent);///FIXME getGrawFrame counts frames from 1 (WRRR!)
  std::cout.clear();
  
  if(!dataFrameRead){
    std::cerr <<KRED<<std::endl<<"ERROR: cannot read file entry: " <<RST<<iEntry<<std::endl
	      <<KRED<<"from file: "<<std::endl
	      <<RST<<tmpFilePath
	      << std::endl;
    std::cerr <<KRED
	      <<"Please check if you are running the application from the resources directory."
	      <<std::endl<<"or if the data file is missing."
	      <<RST
	      <<std::endl;
    exit(1);
  }
  return dataFrameRead;
}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
void EventSourceGRAW::loadEventId(unsigned long int eventId){

  std::cout<<KBLU
	   <<"Start looking for the event id: "<<eventId
	   <<RST<<std::endl;
  
  auto it = myFramesMap.find(eventId);

  if(!isFullFileScanned && it!=myFramesMap.end() &&
     it->second.size()<GRAW_EVENT_FRAGMENTS){
    unsigned int iEntry =  *it->second.rbegin();
    findEventFragments(eventId, iEntry);
  }
  else if(!isFullFileScanned && it==myFramesMap.end() ){
    findEventFragments(eventId, GRAW_EVENT_FRAGMENTS *( eventId- startingEventIndex));
  }
  collectEventFragments(eventId);

  std::cout<<KBLU
	   <<"Finished looking for the event id: "<<eventId
	   <<RST<<std::endl;
}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
void EventSourceGRAW::loadFileEntry(unsigned long int iEntry){

  myCurrentEntry = iEntry;
  loadGrawFrame(iEntry, false);
  unsigned long int eventId = myDataFrame.fHeader.fEventIdx;

  std::cout<<KBLU
	   <<"Looking for event fragments from file entry id: "<<RST<<iEntry
	   <<KBLU<<" event id: "<<RST<<eventId
	   <<std::endl;
  
  if(!isFullFileScanned &&
     (myFramesMap.find(eventId)==myFramesMap.end() ||
      myFramesMap.find(eventId)->second.size()<GRAW_EVENT_FRAGMENTS)){
    findEventFragments(eventId, iEntry);
  }
  collectEventFragments(eventId);
}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
void EventSourceGRAW::checkEntryForFragments(unsigned int iEntry){

  if(myReadEntriesSet.count(iEntry)) return;
  loadGrawFrame(iEntry, false);
  myReadEntriesSet.insert(iEntry);
  unsigned long int currentEventId = myDataFrame.fHeader.fEventIdx;
  int ASAD_idx = myDataFrame.fHeader.fAsadIdx;  
  if(!myASADMap[currentEventId].count(ASAD_idx)) {
    myFramesMap[currentEventId].insert(iEntry);
    myASADMap[currentEventId].insert(ASAD_idx);
  }
}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
void EventSourceGRAW::findEventFragments(unsigned long int eventId, unsigned int iInitialEntry){

  auto it=myFramesMap.find(eventId);
  unsigned int nFragments = 0;
  if(it!=myFramesMap.end()){
    nFragments = it->second.size();
  }

  //bool reachStartOfFile = false;
  //bool reachEndOfFile = false;
  
  unsigned int lowEndScanRange = std::max((unsigned int)0, iInitialEntry-frameLoadRange);
  unsigned int highEndScanRange = std::min((unsigned int)nEntries, iInitialEntry+frameLoadRange);

  for(unsigned int iEntry=iInitialEntry;
      iEntry>=lowEndScanRange && iEntry<nEntries && nFragments<GRAW_EVENT_FRAGMENTS;
      --iEntry){
    checkEntryForFragments(iEntry);
    nFragments =  myFramesMap[eventId].size();
    //reachStartOfFile = (iEntry==0);
    std::cout<<"\r going backward and reading file entry: "<<iEntry
	     <<" fragments found so far: "
	     <<nFragments
	     <<" expected: "<<GRAW_EVENT_FRAGMENTS;
  }
  for(unsigned int iEntry=iInitialEntry;iEntry<highEndScanRange && nFragments<GRAW_EVENT_FRAGMENTS;++iEntry){
    checkEntryForFragments(iEntry);
    nFragments =  myFramesMap[eventId].size();
    //reachEndOfFile = (iEntry==nEntries);
    std::cout<<"\r going forward and reading file entry: "<<iEntry
	     <<" fragments found so far: "
	     <<nFragments
	     <<" expected: "<<GRAW_EVENT_FRAGMENTS;        
  }
  for(unsigned int iEntry=0;iEntry<frameLoadRange && nFragments<GRAW_EVENT_FRAGMENTS;++iEntry){
    checkEntryForFragments(iEntry+nEntries);
    nFragments =  myFramesMap[eventId].size();
    std::cout<<"\r going to next run file and reading file entry: "<<iEntry
	     <<" fragments found so far: "
	     <<nFragments
	     <<" expected: "<<GRAW_EVENT_FRAGMENTS;        
  }
  
  std::cout<<std::endl;
  if(myFramesMap.size()>=nEntries) isFullFileScanned = true;
}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
void EventSourceGRAW::collectEventFragments(unsigned int eventId){

  auto it = myFramesMap.find(eventId);
  if(it==myFramesMap.end()) return;
  if(it->second.size()!=GRAW_EVENT_FRAGMENTS){
      std::cerr<<KRED<<"Fragment counts for eventId = "<<RST<<eventId
	       <<KRED<<" mismatch. Expected: "<<RST<<GRAW_EVENT_FRAGMENTS
	       <<KRED<<" found: "<<RST<<it->second.size()
	       <<RST<<std::endl;
  }
  //long int eventNumberInFile = std::distance(myFramesMap.begin(), it);
 
  myCurrentEvent->Clear();
  myCurrentEvent->SetEventId(eventId);
  myCurrentEvent->SetGeoPtr(myGeometryPtr);

  std::cout<<KYEL<<"Creating a new event with eventId: "<<eventId<<RST<<std::endl;

  for(auto aFragment: it->second){
    loadGrawFrame(aFragment, true);
    myCurrentEvent->SetEventTime(myDataFrame.fHeader.fEventTime);
    myCurrentEvent->SetEventNumber(aFragment/myGeometryPtr->GetAsadNboards());
    int  ASAD_idx = myDataFrame.fHeader.fAsadIdx;
    unsigned long int eventId_fromFrame = myDataFrame.fHeader.fEventIdx;
    if(eventId!=eventId_fromFrame){
      std::cerr<<KRED<<__FUNCTION__
	       <<" Event id mismath!: eventId = "<<eventId
	       <<" eventId_fromFrame: "<<eventId_fromFrame
	       <<RST<<std::endl;
      return;
    }     
    std::cout<<KBLU<<"Found a frame for eventId: "<<RST<<eventId;
    if(aFragment<nEntries) std::cout<<KBLU<<" in file entry: "<<RST<<aFragment<<RST;
    else std::cout<<KBLU<<" in next file entry: "<<RST<<aFragment-nEntries<<RST;
    std::cout<<KBLU<<" for  ASAD: "<<RST<<ASAD_idx<<RST<<std::endl;      
    fillEventFromFrame(myDataFrame);
  }
}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
void EventSourceGRAW::fillEventFromFrame(GET::GDataFrame & aGrawFrame){

  myPedestalCalculator.CalculateEventPedestals(aGrawFrame);

  int  COBO_idx = aGrawFrame.fHeader.fCoboIdx;
  int  ASAD_idx = aGrawFrame.fHeader.fAsadIdx;
  if(ASAD_idx >= myGeometryPtr->GetAsadNboards()){
    std::cout<<KRED<<__FUNCTION__
	     <<" Data format mismatch!. ASAD: "<<ASAD_idx
	     <<" number of ASAD boards in geometry: "<<myGeometryPtr->GetAsadNboards()
	     <<" Skipping the frame."
	     <<RST<<std::endl;
    return;
  }
  
  for (Int_t agetId = 0; agetId < myGeometryPtr->GetAgetNchips(); ++agetId){
    // loop over normal channels and update channel mask for clustering
    for (Int_t chanId = 0; chanId < myGeometryPtr->GetAgetNchannels(); ++chanId){
      int iChannelGlobal     = myGeometryPtr->Global_normal2normal(COBO_idx, ASAD_idx, agetId, chanId);// 0-255 (without FPN)
      GET::GDataChannel* channel = aGrawFrame.SearchChannel(agetId, myGeometryPtr->Aget_normal2raw(chanId));
      if (!channel) continue;
	  
      for (Int_t i = 0; i < channel->fNsamples; ++i){
	GET::GDataSample* sample = (GET::GDataSample*) channel->fSamples.At(i);
	// skip cells outside signal time-window
	Int_t icell = sample->fBuckIdx;
	if(icell<2 || icell>509 || icell<myPedestalCalculator.GetMinSignalCell() || icell>myPedestalCalculator.GetMaxSignalCell()) continue;
	    
	Double_t rawVal  = sample->fValue;
	Double_t corrVal = rawVal;
//	if(myGeometryPtr->GetAsadNboards()==1){
    // Beware HACK!!!
  //TProfile with pedestals is only 256 (max chans in frame) long, pedestals are calculated for each frame and reset
  //to fit into TProfile the global number of first chan in COBO/ASAD has to be substracted from global chanel
  if(removePedestal){
    int minChannelGlobal = myGeometryPtr->Global_normal2normal(COBO_idx, ASAD_idx, 0, 0);
    corrVal -= myPedestalCalculator.GetPedestalCorrection(iChannelGlobal-minChannelGlobal, agetId, icell);
  }
//	} 
	myCurrentEvent->AddValByAgetChannel(COBO_idx, ASAD_idx, agetId, chanId, icell, corrVal);
      }
    }
  }
}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
void EventSourceGRAW::findStartingIndex(unsigned long int size){
  if(nEntries==0){
    startingEventIndex=0;
  } else {
    auto preloadSize=std::min(nEntries,size);
    startingEventIndex=std::numeric_limits<UInt_t>::max();
    for(unsigned long int i=0; i<preloadSize; ++i){
      myFile->GetGrawFrame(myDataFrame,i);
      startingEventIndex=std::min(startingEventIndex, static_cast<unsigned long int>(myDataFrame.fHeader.fEventIdx));
    }
  }
}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
void EventSourceGRAW::configurePedestal(const boost::property_tree::ptree &config){
  auto parser=[this, &config](std::string &&parameter, void (PedestalCalculator::*setter)(int)){
    if(config.find(parameter)!=config.not_found()){
      (this->myPedestalCalculator.*setter)(config.get<int>(parameter));
    }
  };
  parser("minPedestalCell", &PedestalCalculator::SetMinPedestalCell);
  parser("maxPedestalCell", &PedestalCalculator::SetMaxPedestalCell);
  parser("minSignalCell", &PedestalCalculator::SetMinSignalCell);
  parser("maxSignalCell", &PedestalCalculator::SetMaxSignalCell);
}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
std::string EventSourceGRAW::getNextFilePath(){

  std::string token = "_";
  int index = myFilePath.rfind(token);
  int fieldLength = 4;
  std::string fileIndex = myFilePath.substr(index+token.size(), fieldLength);
  int nextFileIndex = std::stoi(fileIndex) + 1;  
  ///
  token = "_";
  index = myFilePath.rfind(token);
  std::string fileNamePrefix = myFilePath.substr(0,index+token.size());
  ///
  std::ostringstream ostr;
  ostr<<fileNamePrefix<<std::setfill('0') <<std::setw(4)<<nextFileIndex<<".graw";
  std::string nextFilePath = ostr.str();
  return nextFilePath;
}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
