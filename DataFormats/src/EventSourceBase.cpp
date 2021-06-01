#include <cstdlib>
#include <iostream>
#include <fstream>

#include "colorText.h"
#include "EventSourceBase.h"
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
EventSourceBase::EventSourceBase() {

  myCurrentEvent =  std::make_shared<EventTPC>();
  myCurrentEntry = 0;
  nEntries = 0;
}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
EventSourceBase::~EventSourceBase() { }
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
void EventSourceBase::loadDataFile(const std::string & fileName){

  if(!std::ifstream(fileName)){
    std::cout<<KRED<<"Input data file: "<<RST<<fileName<<KRED<<" not found!"<<RST<<std::endl;
    exit(1);
  }

  currentFilePath = fileName;
}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
void EventSourceBase::loadGeometry(const std::string & fileName){
  
  myGeometryPtr = std::make_shared<GeometryTPC>(fileName.c_str(), false);
  if(!myGeometryPtr){
    std::cerr<<"Geometry not loaded! Refuse to work anymore."<<std::endl;
    exit(0);
  }
}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
std::shared_ptr<EventTPC> EventSourceBase::getCurrentEvent() const{

  return myCurrentEvent;
}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
std::shared_ptr<EventTPC> EventSourceBase::getLastEvent(){

  if(nEntries>0){
    loadFileEntry(nEntries-1);
    myCurrentEntry = nEntries-1;
  }

  return myCurrentEvent;
}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
unsigned long int EventSourceBase::currentEntryNumber() const{

  return myCurrentEntry;
}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
std::shared_ptr<GeometryTPC> EventSourceBase::getGeometry() const{ return myGeometryPtr; }
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
unsigned long int EventSourceBase::numberOfEntries() const{ return nEntries; }
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
unsigned long int EventSourceBase::currentEventNumber() const{

  if(myCurrentEvent){
    return myCurrentEvent->GetEventId();
  }
  return 0;
}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
std::string EventSourceBase::getCurrentPath() const{

  return currentFilePath;
}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
std::shared_ptr<EventTPC> EventSourceBase::getNextEventLoop(){
  unsigned int currentEventIdx;
  do{
    currentEventIdx=myCurrentEvent->GetEventId();
    getNextEvent();
  }
  while(!eventFilter.pass(*myCurrentEvent) && currentEventIdx!=myCurrentEvent->GetEventId());
  return myCurrentEvent;
}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
std::shared_ptr<EventTPC> EventSourceBase::getPreviousEventLoop(){
  unsigned int currentEventIdx;
  do{
    currentEventIdx=myCurrentEvent->GetEventId();
    getPreviousEvent();
  }
  while(!eventFilter.pass(*myCurrentEvent) && currentEventIdx!=myCurrentEvent->GetEventId());
  return myCurrentEvent;
}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
