#include <vector>
#include <string>
#include <iostream>
#include <algorithm>

#include "TFile.h"
#include "TTree.h"
#include "TVector3.h"
#include "TH2Poly.h"

#include "GeometryTPC.h"
#include "Track3D.h"
#include "TrackSegment3D.h"
#include "EventInfo.h"
#include "HIGS_trees_analysis_kris.h"
#include "UtilsEventInfo.h"

#include "colorText.h"

using std::chrono::duration;
using std::chrono::duration_cast;

///////////////////////////////
///////////////////////////////
HIGS_trees_analysis_kris::HIGS_trees_analysis_kris(std::shared_ptr<GeometryTPC> aGeometryPtr, // definition of LAB detector coordinates
					 float beamEnergy,   // nominal gamma beam energy [MeV] in detector LAB frame
					 TVector3 beamDir,   // nominal gamma beam direction in detector LAB frame
					 double pressure,    // CO2 pressure [mbar]
					 double temperature){// CO2 temperature [K]
  setGeometry(aGeometryPtr);
  setBeamProperties(beamEnergy, beamDir);
  setIonRangeCalculator(pressure, temperature);
  setCuts();
  open();
//  bookHistos();

}

void HIGS_trees_analysis_kris::setCuts(){
  xyAreaCut.initialize(myGeometryPtr, 5.0); // TODO - TO BE PARAMETERIZED
  zRangeCut.initialize(myGeometryPtr, 25, 5); // TODO - TO BE PARAMETERIZED
}






bool HIGS_trees_analysis_kris::eventFilter(Track3D *aTrack){

  //  return true;


  // print statistics on demand
  const auto printAccepted=false; // TODO - TO BE PARAMETERIZED !!!
  const auto printRejected=false; // TODO - TO BE PARAMETERIZED !!!

  static TrackSegment3DCollection list;
  static TVector3 vertexPos;
  bool result=true;

  // cut #1 :reject empty events
  if(!aTrack || (list=aTrack->getSegments()).size()==0) {
    result=false;
    if(printRejected) {
      std::cout<<KRED<<__FUNCTION__<<": REJECTED (empty event)"<<RST<<std::endl;
    }
  }

  // get sorted list of tracks (descending order by track length)
  std::sort(list.begin(), list.end(),
	    [](const TrackSegment3D& a, const TrackSegment3D& b) {
	      return a.getLength() > b.getLength();
	    });


  
  // cut #2 : XY plane : minimal distance to the border of UVW active area
  // - less strict than simple XY rectangular cut, allows to gain some statistics
  if(result) {
    for(auto &track: list) {
      if( !xyAreaCut.IsInside(track.getEnd().X(), track.getEnd().Y()) ||
	  !xyAreaCut.IsInside(track.getStart().X(), track.getStart().Y()) ) {
	result=false;
	if(printRejected) {
	  std::cout<<KRED<<__FUNCTION__<<": REJECTED (horizontal: track too close to UVW border)"<<RST<<std::endl;
	}
	break;
      }
    }
  }

    // cut #3 : XY plane : vertex position per event, corrected for beam tilt
  if(result) {
    vertexPos = list.front().getStart();
    if( fabs( vertexPos.Y() - (beam_offset+beam_slope*vertexPos.X()) ) > 0.5*beam_diameter ) {
      result=false;
      if(printRejected) {
	std::cout<<KRED<<__FUNCTION__<<": REJECTED (horizontal: vertex too far from beam axis)"<<RST<<std::endl;
      }
    }
  }



                                                                                                                                                 

  // cut #4 : global Z-span per event, verifies that:
  // - vertical projection length is below physical drift cage length
  // - tracks do not overlap with pedestal exclusion zone, begin of history buffer
  // - tracks not too close to end of history buffer
  if(result) {
    auto zmin=list.at(0).getStart().Z(), zmax=zmin;
    for(auto i=0u; i<list.size(); i++) {
      zmin=std::min(zmin, (std::min(list.at(i).getStart().Z(), list.at(i).getEnd().Z())));
      zmax=std::max(zmax, (std::max(list.at(i).getStart().Z(), list.at(i).getEnd().Z())));
    }
    if(!zRangeCut.IsInside(zmin, zmax)) {
      result=false;
      if(printRejected) {
	std::cout<<KRED<<__FUNCTION__<<": REJECTED (vertical: too close to electronics limits / too close to drift cage height)"<<RST<<std::endl;
      }
    }
  }

  // cut #5 : Z-span wrt vertex per track per event, verifies that:
  // - vertical distance of endpoint to vertex is less than half of drift cage height
  //   corrected for maximal vertical beam spread
  // - ensures that 2,3-prong events hit neither the GEM plane nor the cathode plane
  // NOTE: does not protect against 1-prong events (eg. background) originating
  //       from the GEM plane or the cathode plane
  if(result) {
    for(auto i=0u; i<list.size(); i++) {
      if(fabs(list.at(i).getEnd().Z()-vertexPos.Z())>0.5*(myGeometryPtr->GetDriftCageZmax()-myGeometryPtr->GetDriftCageZmin()-beam_diameter)) {
      result=false;
	if(printRejected) {
	  std::cout<<KRED<<__FUNCTION__<<": REJECTED (vertical: too long wrt vertex)"<<RST<<std::endl;
	}
	break;
      }
    }
  }

 

  if(printAccepted && result) {
    std::cout<<KGRN<<__FUNCTION__<<": ACCEPTED"<<RST<<std::endl;
  }
  return result;
}



///////////////////////////////
///////////////////////////////
HIGS_trees_analysis_kris::~HIGS_trees_analysis_kris(){
  Output3prongTreePtr->Write("", TObject::kOverwrite);
  close();
}
///////////////////////////////
///////////////////////////////
void HIGS_trees_analysis_kris::open(){

  std::string fileName = "Trees_3p.root";
  
  std::string tree3Name = "Tree_3prong_events";
  OutputFilePtr = std::make_shared<TFile>(fileName.c_str(),"RECREATE");
  if(!OutputFilePtr) {
    std::cout<<KRED<<"HIGS_trees_analysis::open: Cannot create new ROOT file: "<<RST<<fileName
	     <<std::endl;
    return;
  }
  Output3prongTreePtr = std::make_shared<TTree>(tree3Name.c_str(),"");
  
  Output1prongTreePtr->Branch("data", &event1prong_);
  Output2prongTreePtr->Branch("data", &event2prong_);
  Output3prongTreePtr->Branch("data", &event3prong_);
    
}
///////////////////////////////
///////////////////////////////
void HIGS_trees_analysis_kris::close(){
  if(!OutputFilePtr){
    std::cout<<KRED<<"HIGS_trees_analysis::close: "<<RST
	     <<" pointer to output file not set!"
	     <<std::endl;
    return;
  }
  Output2prongTreePtr->Write("", TObject::kOverwrite);
  OutputFilePtr->Close();
}
///////////////////////////////
///////////////////////////////


void HIGS_trees_analysis_kris::setGeometry(std::shared_ptr<GeometryTPC> aGeometryPtr){
  myGeometryPtr = aGeometryPtr;
  if(!myGeometryPtr) {
    std::cout<<KRED<<"HIGS_trees_analysis::setGeometry: "<<RST
	     <<" pointer to TPC geometry not set!"
	     <<std::endl;
    exit(-1);
  }
}


///////////////////////////////
///////////////////////////////
void HIGS_trees_analysis_kris::setBeamProperties(float beamEnergy,   // nominal gamma beam energy [MeV] in detector LAB frame
					    TVector3 beamDir) { // nominal gamma beam direction in detector LAB frame
  
  //!!
  photonEnergyInMeV_LAB = fabs(beamEnergy);
  photonUnitVec_DET_LAB = beamDir.Unit();
  beam_slope=tan(3.0e-3); // [rad], measured slope: Y_DET(X_DET)=offset+slope*X_DET
  beam_offset=-0.73; // [mm], measured offset: Y_DET of beam axis at X_DET=0
  beam_diameter=11.8; // [mm] // TODO - TO BE PARAMETERIZED !!!
  //!

}
//////////////////////////
//////////////////////////
void HIGS_trees_analysis_kris::setIonRangeCalculator(double pressure, double temperature){ // CO2 pressure [mbar] and temperature [K]

  // set current conditions: gas=CO2, arbitrary temperature [K] and pressure [mbar]
  myRangeCalculator.setGasConditions(/*IonRangeCalculator::*/CO2, fabs(pressure), fabs(temperature));
}
///////////////////////////////
///////////////////////////////
void HIGS_trees_analysis_kris::fillTrees(Track3D *aTrack, eventraw::EventInfo *aEventInfo){
  
  if( !eventFilter(aTrack) ) return;
  static bool isFirst_3prong=true;
  const int ntracks = aTrack->getSegments().size();
  if (ntracks==3) fillTrees3prong(aTrack, aEventInfo, isFirst_3prong);  
}

///////////////////////////////
///////////////////////////////
void HIGS_trees_analysis_kris::fillTrees3prong(Track3D *aTrack, eventraw::EventInfo *aEventInfo, bool & isFirst){

  if(!Output3prongTreePtr){
    std::cout<<KRED<<"HIGS_trees_analysis::fillTrees"<<RST
	     <<" pointer to 3 prong output tree not set!"
	     <<std::endl;
    return;
  }

  // get sorted list of tracks (descending order by track length)
  TrackSegment3DCollection list = aTrack->getSegments();
  std::sort(list.begin(), list.end(),
	    [](const TrackSegment3D& a, const TrackSegment3D& b) {
	      return a.getLength() > b.getLength();
	    });  
 
  event3prong_->runId=(aEventInfo ? aEventInfo->GetRunId() : -1);
  event3prong_->eventId=(aEventInfo ? aEventInfo->GetEventId() : -1);
 
  event3prong_->eventType = aEventInfo ? aEventInfo->GetEventType().to_ulong() : 0;
  event1prong_->unixTimeSec =
      (aEventInfo
           ? duration_cast<duration<long double>>(
                 tpcreco::eventAbsoluteTime(*aEventInfo).time_since_epoch())
                 .count()
           : -1); // absolute Unix time [s]
  static double last_timestamp = 0;
  if(isFirst) {
    last_timestamp=event3prong_->unixTimeSec;
    isFirst=false;
  }
  event1prong_->runTimeSec =
      (aEventInfo ? duration_cast<duration<long double>>(
                        tpcreco::eventRelativeTime(*aEventInfo))
                        .count()
                  : -1); // [s]
  event3prong_->deltaTimeSec=(aEventInfo ? event3prong_->unixTimeSec - last_timestamp : -1); // [s] time difference for rate measurements
  last_timestamp=event3prong_->unixTimeSec;

  event3prong_->vertexPos = list.front().getStart();
  auto track1=list.at(0);
  event3prong_->alpha1_endPos = track1.getEnd();
  event3prong_->alpha1_length = track1.getLength(); // longest alpha
  event3prong_->alpha1_phiDET = track1.getTangent().Phi();
  event3prong_->alpha1_cosThetaDET = cos(track1.getTangent().Theta());
  event3prong_->alpha1_phiBEAM = atan2(-track1.getTangent().Z(), track1.getTangent().Y()); // [rad], azimuthal angle from horizontal axis;
  event3prong_->alpha1_cosThetaBEAM = track1.getTangent()*photonUnitVec_DET_LAB; // cos polar angle wrt beam axis
  event3prong_->alpha1_chargeProfile = track1.getChargeProfile();

  auto track2=list.at(1);
  event3prong_->alpha2_endPos = track2.getEnd();
  event3prong_->alpha2_length = track2.getLength(); // middle alpha
  event3prong_->alpha2_phiDET = track2.getTangent().Phi();
  event3prong_->alpha2_cosThetaDET = cos(track2.getTangent().Theta());
  event3prong_->alpha2_phiBEAM = atan2(-track2.getTangent().Z(), track2.getTangent().Y()); // [rad], azimuthal angle from horizontal axis;
  event3prong_->alpha2_cosThetaBEAM = track2.getTangent()*photonUnitVec_DET_LAB; // cos polar angle wrt beam axis
  event3prong_->alpha2_chargeProfile = track2.getChargeProfile();

  auto track3=list.at(2);
  event3prong_->alpha3_endPos = track3.getEnd();
  event3prong_->alpha3_length = track3.getLength(); // shortest alpha
  event3prong_->alpha3_phiDET = track3.getTangent().Phi();
  event3prong_->alpha3_cosThetaDET = cos(track3.getTangent().Theta());
  event3prong_->alpha3_phiBEAM = atan2(-track3.getTangent().Z(), track3.getTangent().Y()); // [rad], azimuthal angle from horizontal axis;
  event3prong_->alpha3_cosThetaBEAM = track3.getTangent()*photonUnitVec_DET_LAB; // polar angle wrt beam axis
  event3prong_->alpha3_chargeProfile = track3.getChargeProfile();
  


  Output3prongTreePtr->Fill();
}
