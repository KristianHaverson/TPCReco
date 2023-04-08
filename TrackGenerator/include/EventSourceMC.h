#ifndef _EventSourceMC_H_
#define _EventSourceMC_H_

#include <memory>
#include <boost/property_tree/json_parser.hpp>

#include "TH1F.h"
#include "TRandom3.h"

#include "EventSourceBase.h"
#include "UVWprojector.h"
#include "Track3D.h"
#include "IonRangeCalculator.h"
#include "TLorentzVector.h"
#include "TLorentzRotation.h"
#include "IonRangeCalculator.h"
#include "TVector3.h"

// *********************************** //
// *********************************** //
// *********************************** //
// *********************************** //

class Prong3_TrackInfo {
private:
  double trackLength[3];
  double trackTheta[3];
  double trackPhi[3];
  double trackEnergy[3];
  TVector3 tangentA1;
  TVector3 tangentA2;
  TVector3 tangentA3;

public:
  // Setters for length, angle, and energy of each alpha particle

  void setTangent(int index, TVector3 tangent) {
    if(index==0)tangentA1=tangent;
    if(index==1)tangentA2=tangent;
    if(index==2)tangentA3=tangent;
  
  }

  TVector3 getTangent(int index) {
    if(index==0)return tangentA1;
    if(index==1)return tangentA2;
    if(index==2)return tangentA3;
    return tangentA1;
  }


  void setTrackLength(int index, double value) {
    trackLength[index] = value;
  }
  void setTheta(int index, double value) {
    trackTheta[index] = value;
  }
  void setPhi(int index, double value) {
    trackPhi[index] = value;
  }
  void setEnergy(int index, double value) {
    trackEnergy[index] = value;
  }

  // Getters for length, angle, and energy of each alpha particle
  double getTrackLength(int index) {
    return trackLength[index];
  }
  double getPhi(int index) {
    return trackPhi[index];
  }
    double getTheta(int index) {
    return trackTheta[index];
  }
  double getEnergy(int index) {
    return trackEnergy[index];
  }
};



// *********************************** //
// *********************************** //
// *********************************** //
// *********************************** //

class Prong2_TrackInfo {
private:
  double trackLength[2];
  double trackTheta[2];
  double trackPhi[2];
  double trackEnergy[2];
  TVector3 tangentI1;
  TVector3 tangentI2;


  
public:

  //TLorentzVector PRONG2_CARBON_LAB;
  //TLorentzVector PRONG2_ALPHA_LAB;

  //void set2pCarbon(TLorentzVector vec){ PRONG2_CARBON_LAB = vec; }
  //void set2pAlpha(TLorentzVector& vec) { PRONG2_ALPHA_LAB = vec; }
  // Setters for length, angle, and energy of each alpha particle

  void setTangent(int index, TVector3 tangent) {
    if(index==0)tangentI1=tangent;
    if(index==1)tangentI2=tangent;
  
  }

  TVector3 getTangent(int index) {
    if(index==0)return tangentI1;
    if(index==1)return tangentI2;
    return tangentI2;
  }


  void setTrackLength(int index, double value) {
    trackLength[index] = value;
  }
  void setTheta(int index, double value) {
    trackTheta[index] = value;
  }
  void setPhi(int index, double value) {
    trackPhi[index] = value;
  }
  void setEnergy(int index, double value) {
    trackEnergy[index] = value;
  }
  // Getters for length, angle, and energy of each alpha particle
  double getTrackLength(int index) {
    return trackLength[index];
  }
  double getPhi(int index) {
    return trackPhi[index];
  }
    double getTheta(int index) {
    return trackTheta[index];
  }
  double getEnergy(int index) {
    return trackEnergy[index];
  }
};

// *********************************** //
// *********************************** //
// *********************************** //
// *********************************** //




class EventSourceMC: public EventSourceBase {
  
public:

  EventSourceMC(){};
  EventSourceMC(const std::string & geometryFileName);
  ~EventSourceMC();

  void loadFileEntry(unsigned long int iEntry);
  void loadEventId(unsigned long int iEvent);
  void loadDataFile(const std::string & fileName);

  std::shared_ptr<EventTPC> getNextEvent();
  std::shared_ptr<EventTPC> getPreviousEvent();
  std::shared_ptr<EventTPC> getLastEvent();

  const Track3D & getGeneratedTrack(unsigned int index=0) const {return myTracks3D.at(index);}
  pid_type getGeneratedEventType() const {return myGenEventType;}
  unsigned long int numberOfEvents() const;
  void loadGeometry(const std::string & fileName);
  TLorentzVector getMomentumVec()const {return myPvec;}



  void setGammaEnergy(double energy){GAMMA_ENERGY=energy;}
  double getGammaEnergy(){return GAMMA_ENERGY;}

  void setChamberPressure(double pressure){CHAMBER_PRESSURE=pressure;}
  double getChamberPressure(){return CHAMBER_PRESSURE;}


  // ------------------------------------------------- //
  bool DEBUG = true;

  double GAMMA_ENERGY = 0.0;
  double CHAMBER_PRESSURE = 0.0;

  const double CONV_Kev2MeV = 1./1000.;

  const double AMU        = 931.49410242;

  const double MASS_8Be   = AMU  * 8.00530510000; 
	const double MASS_4He   = AMU  * 4.00260325413; 
  const double MASS_12C   = AMU  *12.00000000000; 
  const double MASS_16O = 15.99491461956 * AMU; 


  const double Q_VALUE_16O  = - 7.162; 
  const double Q_VALUE_12C  = - 7367.0 * CONV_Kev2MeV; 
  const double Q_VALUE_8BE  =   92.0   * CONV_Kev2MeV; 


  


  // ------------------------------------------------- //


  private:
  void generateEvent();
  TVector3 createVertex() const;

  private:
  void generateTwoProng();
    double SRIMNormtrackLength( pid_type ion_id)const;
    Track3D createTrack(const TVector3 & aVtx, pid_type ion_id,TVector3 &aTangentNonBoost, TLorentzVector &p4outPut,Prong2_TrackInfo Prong2_TrackInfo ) const;

    Track3D createTrack(const TVector3 & aVtx, pid_type ion_id, Prong2_TrackInfo &Prong2_TrackInfo, const int index ) const;
    TrackSegment3D create2pSegment(const TVector3 vertexPos, pid_type ion_id,Prong2_TrackInfo &Prong2_TrackInfo, const int index) const;
    double setInfoFromSRIM_2p( pid_type ion_id,Prong2_TrackInfo &Prong2_TrackInfo, const int index)const;

    TrackSegment3D Prong2_createSegmentBoost(const TVector3 vertexPos, pid_type ion_id,TVector3 &aTangentNonBoost, TLorentzVector &p4outPut,Prong2_TrackInfo Prong2_TrackInfo) const;  
    double Prong2_Boost2Lab(double const gammaEnergy,double const length,TVector3 &tangent, pid_type ion_id, TLorentzVector &output,Prong2_TrackInfo Prong2_TrackInfo)const;

  
  private:
  void generateThreeProng(const int state);
    Track3D createTrack(const TVector3 & aVtx, pid_type ion_id,Prong3_TrackInfo &Prong3_TrackInfo, const int index , const int groun_or_exc) const;
    TrackSegment3D create3pSegment(const TVector3 vertexPos, pid_type ion_id,Prong3_TrackInfo &Prong3_TrackInfo, const int index, const int groun_or_exc ) const;
    double setInfoFromSRIM_3p( pid_type ion_id,Prong3_TrackInfo &Prong3_TrackInfo, const int index)const;

  




  TLorentzVector myPvec;
  TGraph* braggGraph_alpha, *braggGraph_12C;
  double braggGraph_alpha_energy, braggGraph_12C_energy;
  double keVToChargeScale{1.0};
  mutable TRandom3 myRndm{0};
  std::shared_ptr<UVWprojector> myProjectorPtr;
  mutable IonRangeCalculator myRangeCalculator;
  TH3D my3DChargeCloud;
  std::vector<Track3D> myTracks3D;
  pid_type myGenEventType;
  void fill3DChargeCloud(const Track3D & aTrack);
  void fillPEventTPC(const TH3D & h3DChargeCloud, const Track3D & aTrack);
  void generateSingleProng(pid_type ion_id=pid_type::ALPHA);
  TrackSegment3D createSegment(const TVector3 vertexPos, pid_type ion_id) const;  
  TH1F createChargeProfile(double ion_range, pid_type ion_id) const;  
  Track3D createTrack(const TVector3 & aVtx, pid_type ion_id) const;
  



  
};
#endif

