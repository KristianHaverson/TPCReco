#include<cstdlib>
#include <iostream>

#include "EventSourceBase.h"
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
EventSourceBase::EventSourceBase() {

  myCurrentEvent =  std::make_shared<EventTPC>();
  myCurrentEntry = 0;
  nEvents = 0;
}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
EventSourceBase::~EventSourceBase() { }
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
void EventSourceBase::loadDataFile(const std::string & fileName){

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
void EventSourceBase::loadEventId(unsigned long int iEvent){

  unsigned long int iEntry = 0;
  while(currentEventNumber()!=iEvent && iEntry<nEvents){  
    loadFileEntry(iEntry);
    ++iEntry;
  }
}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
std::shared_ptr<EventTPC> EventSourceBase::getCurrentEvent() const{

  return myCurrentEvent;
}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
std::shared_ptr<EventTPC> EventSourceBase::getNextEvent(){

  std::cout<<"nEvents: "<<nEvents
	   <<" myCurrentEntry: "<<myCurrentEntry
	   <<std::endl;

  if(nEvents>0 && myCurrentEntry<nEvents-1){
    loadFileEntry(++myCurrentEntry);
  }
  return myCurrentEvent;
}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
std::shared_ptr<EventTPC> EventSourceBase::getPreviousEvent(){

  if(myCurrentEntry>0 && nEvents>0){
    loadFileEntry(--myCurrentEntry);
  }
  return myCurrentEvent;
}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
std::shared_ptr<EventTPC> EventSourceBase::getLastEvent(){

  if(nEvents>0){
    loadFileEntry(nEvents-1);
    myCurrentEntry = nEvents-1;
  }

  return myCurrentEvent;
}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
std::shared_ptr<GeometryTPC> EventSourceBase::getGeometry() const{ return myGeometryPtr; }
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
unsigned long int EventSourceBase::numberOfEvents() const{ return nEvents; }
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
