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


class AlphaData {
private:
  double alphaLength[3];
  double alphaTheta[3];
  double alphaPhi[3];
  double alphaEnergy[3];
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


  void setAlphaLength(int index, double value) {
    alphaLength[index] = value;
  }
  void setAlphaTheta(int index, double value) {
    alphaTheta[index] = value;
  }
  void setAlphaPhi(int index, double value) {
    alphaPhi[index] = value;
  }
  void setAlphaEnergy(int index, double value) {
    alphaEnergy[index] = value;
  }
  // Getters for length, angle, and energy of each alpha particle
  double getAlphaLength(int index) {
    return alphaLength[index];
  }
  double getAlphaPhi(int index) {
    return alphaPhi[index];
  }
    double getAlphaTheta(int index) {
    return alphaTheta[index];
  }
  double getAlphaEnergy(int index) {
    return alphaEnergy[index];
  }
};



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

  //std::shared_ptr<EventTPC> getCurrentEvent() const;



  
 private:
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

  TVector3 createVertex() const;
  TrackSegment3D createLorentzSegment(const TVector3 vertexPos, pid_type ion_id,TVector3 &aTangentNonBoost, TLorentzVector &p4outPut) const;  
  TrackSegment3D createSegment(const TVector3 vertexPos, pid_type ion_id) const;  
  TH1F createChargeProfile(double ion_range, pid_type ion_id) const;  
  Track3D createTrack(const TVector3 & aVtx, pid_type ion_id) const;
  Track3D createTrack(const TVector3 & aVtx, pid_type ion_id,TVector3 &aTangentNonBoost, TLorentzVector &p4outPut) const;
  Track3D createTrack(const TVector3 & aVtx, pid_type ion_id,AlphaData &alphaData, const int index , const int groun_or_exc) const;
  
TrackSegment3D create3pSegment(const TVector3 vertexPos, pid_type ion_id,AlphaData &alphaData, const int index, const int groun_or_exc ) const;
    void generateSingleProng(pid_type ion_id=pid_type::ALPHA);
  void generateTwoProng();
  void generateThreeProng();

  void fill3DChargeCloud(const Track3D & aTrack);
  void fillPEventTPC(const TH3D & h3DChargeCloud, const Track3D & aTrack);
  void generateEvent();

  
  double Boost2Lab(double const gammaEnergy,double const length,TVector3 &tangent, pid_type ion_id, TLorentzVector &output)const;
  double boostTheta(double thetaCM);
  double boostPhi(double phiCM);
  void SRIMtrackLength(const double gammaEnergy, const double pressure, pid_type ion_id, double &minLen, double &maxLen)const;
  double SRIMNormtrackLength(const double gammaEnergy, const double pressure, pid_type ion_id)const;
 // double SRIMNormtrackLength(const double gammaEnergy, const double pressure, pid_type ion_id, const int index)const;
double SRIMNormtrackLength(const double gammaEnergy, const double pressure, pid_type ion_id,AlphaData &alphaData, const int index)const;
double BoostbyVel(double const gammaEnergy,AlphaData &alphaData,const int index,TVector3 &tangent)const;
double Boost2Lab3p(double const gammaEnergy,AlphaData &alphaData, int index)const;
double BoostByVel3p(double const gammaEnergy,AlphaData &alphaData, int index)const;
};
#endif

