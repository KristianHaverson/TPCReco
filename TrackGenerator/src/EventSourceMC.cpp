#include <iostream>

#include "colorText.h"
#include "EventSourceMC.h"
#include "TLorentzVector.h"
#include "TLorentzRotation.h"
#include "CommonDefinitions.h"
#include "IonRangeCalculator.h"

/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
EventSourceMC::EventSourceMC(const std::string & geometryFileName) {

  loadGeometry(geometryFileName);

  int NbinsX=50, NbinsY=50, NbinsZ=50;
  double xmin=-150,  ymin=-100, zmin=-100;
  double xmax=150,  ymax=100,  zmax=100;

  braggGraph_alpha = new TGraph("dEdx_corr_alpha_10MeV_CO2_250mbar.dat", "%lg %lg");
  braggGraph_12C = new TGraph("dEdx_corr_12C_5MeV_CO2_250mbar.dat", "%lg %lg");
  double nominalPressure = 250.0; //[mbar]
  braggGraph_alpha_energy = 10; // [MeV]
  braggGraph_12C_energy = 5; // [MeV]

  keVToChargeScale = 100.0; // 1 keV = 100 charge in arb. units
  
  myRangeCalculator.setGasPressure(nominalPressure);


  my3DChargeCloud = TH3D("h3DChargeCloud", "charge cloud [arb. units];X [mm];Y [mm];Z [mm]",
			 NbinsX, xmin, xmax,
			 NbinsY, ymin, ymax,
			 NbinsZ, zmin, zmax);
  
  nEntries = 9999;
}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
EventSourceMC::~EventSourceMC(){ }
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
void EventSourceMC::loadDataFile(const std::string & fileName){ }
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
void EventSourceMC::loadFileEntry(unsigned long int iEntry){

  generateEvent();
  myCurrentEntry = iEntry;
  
}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
unsigned long int EventSourceMC::numberOfEvents() const{ return nEntries; }
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
std::shared_ptr<EventTPC> EventSourceMC::getNextEvent(){

  generateEvent();
  return myCurrentEvent;
}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
std::shared_ptr<EventTPC> EventSourceMC::getPreviousEvent(){

  generateEvent();
  return myCurrentEvent;
}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
void EventSourceMC::loadEventId(unsigned long int iEvent){

  myCurrentEntry = iEvent - 1;
  generateEvent();
 
}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
void EventSourceMC::loadGeometry(const std::string & fileName){

  EventSourceBase::loadGeometry(fileName);
  myProjectorPtr.reset(new UVWprojector(myGeometryPtr));
}

/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
TrackSegment3D EventSourceMC::createSegment(const TVector3 vertexPos, pid_type ion_id) const{

  double length = 60;
  double theta = 0, phi = 0.0;
  double minCosTheta = -1, maxCosTheta = 1;
  double minPhi = 0, maxPhi = 2*M_PI;
  double minLength = 20, maxLength = 60;



  if(ion_id==pid_type::ALPHA){
    minLength = 20;
    maxLength = 60;    
  }
  else if(ion_id==pid_type::C_12 && myTracks3D.size()==1){
    minLength = 10;
    maxLength = 15;
    const TVector3 &aTangent = -myTracks3D.front().getSegments().front().getTangent();
    minPhi = aTangent.Phi()-0.1;
    maxPhi = aTangent.Phi()+0.1;
    minCosTheta = cos(aTangent.Theta())-0.1;
    maxCosTheta = cos(aTangent.Theta())+0.1;
    if(minCosTheta<-1) minCosTheta = -1.0;
    if(maxCosTheta>1) maxCosTheta = 1.0;
  }
  else if(ion_id==pid_type::THREE_ALPHA){
    minLength = 5;
    maxLength = 20;
  }
 
  if(false && myCurrentEntry==0){
    theta = 0.0;
    phi = 0.0;
  }
  else if(false && myCurrentEntry==1){
    theta = M_PI/2.0;
    phi = M_PI/4.0;
  }
  else if(false && myCurrentEntry==2){
    theta = M_PI/4.0;
    phi = M_PI/4.0;
  }
  else{
    theta = TMath::ACos(myRndm.Uniform(minCosTheta, maxCosTheta));
    phi = myRndm.Uniform(minPhi, maxPhi);
    length = myRndm.Uniform(minLength, maxLength);
  }

  TVector3 tangent;
  tangent.SetMagThetaPhi(1.0, theta, phi);

  TrackSegment3D aSegment;
  aSegment.setGeometry(myGeometryPtr);
  aSegment.setStartEnd(vertexPos, vertexPos + tangent*length);
  aSegment.setPID(ion_id);
  
  return aSegment;
}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////


//!------------------------------------------------------
TVector3 EventSourceMC::createVertex() const{

  double x = 0.0, y = 0.0, z = 0.0;

  std::random_device rd; // Seed random number
  std::mt19937 gen1(rd());
  std::mt19937 gen2(rd());
  std::mt19937 gen3(rd());
  double xmin=-150;
  double xmax=150;
  std::uniform_real_distribution<> dis(xmin, xmax);
  std::normal_distribution<> disGaus(0, 3);
  
  x = dis(gen1); 
  y = disGaus(gen2);
  z = disGaus(gen3);
  std::cout<<"-> "<<x<<" "<<y<<" "<<z<<std::endl;



  TVector3 aVertex(x,y,z);
  return aVertex;
}
//!------------------------------------------------------

/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////

//!------------------------------------------------------
double EventSourceMC::Boost2Lab(double const length,TVector3 &tangent, pid_type ion_id)const{
  
  // ** KE ** //
  IonRangeCalculator myRangeCalculator(gas_mixture_type::CO2,250,293.15);
  double alphaEnergyKe = myRangeCalculator.getIonEnergyMeV(pid_type::ALPHA,length);
  double carbonEnergyKe =myRangeCalculator.getIonEnergyMeV(pid_type::CARBON_12, length);
 
  // ** ION INFO ** //
  const double AMU = 931.49410242; // MeV/c^2
  const double alphaMass  =  4.00260325413  *AMU ;
  const double oxygenMass = 15.99491461956*AMU; 
  const double carbonMass = 12.000*AMU; 

  // ** BEAM INFO ** //
  TVector3 beamDir_LAB(-1,0, 0); // unit vector
  double beamEnergy = 13.1;
  TLorentzVector p4_beam_CMS(beamDir_LAB.Unit()*beamEnergy, beamEnergy); // [MeV/c, MeV/c, MeV/c, MeV/c^2]
  TLorentzVector p4_target_CMS(0, 0, 0, oxygenMass); // [MeV/c, MeV/c, MeV/c, MeV/c^2]
  TLorentzVector p4_total_CMS =  p4_beam_CMS + p4_target_CMS;

  // ** BEAM CHECK ** //
  double totalEnergy_CMS=(p4_target_CMS+p4_beam_CMS).E();
  double totalEnergy_CMS_check=sqrt( oxygenMass*(2*beamEnergy+oxygenMass) ); 

  std::cout<<"Total Energy CMS: "<<KRED<<totalEnergy_CMS<<" <> "<<totalEnergy_CMS_check<<RST<<std::endl;



  double p_alpha_CMS =  sqrt(alphaEnergyKe*(alphaEnergyKe+2*alphaMass));
  double p_carbon_CMS = sqrt(carbonEnergyKe*(carbonEnergyKe+2*carbonMass));
  double Etot_alpha_CMS  = alphaMass+alphaEnergyKe;
  double Etot_carbon_CMS  = carbonMass+carbonEnergyKe;

  //std::cout<<"Total Energy Alpha  CMS: "<<KRED<<Etot_alpha_CMS<<RST<<std::endl;
  //std::cout<<"Total Energy Carbon CMS: "<<KRED<<Etot_carbon_CMS<<RST<<std::endl;
  //std::cout<<"Total Energy CMS xcheck: "<<KRED<<Etot_carbon_CMS+Etot_alpha_CMS <<RST<<std::endl;


  
  TVector3 boost_CMS_to_LAB = p4_total_CMS.BoostVector(); 
  TLorentzRotation rot_CMS_to_LAB(boost_CMS_to_LAB);
  double energyPostBoost=0;

  if(ion_id==pid_type::ALPHA){
    TLorentzVector alpha_CMS(tangent*p_alpha_CMS,Etot_alpha_CMS);
    TLorentzVector alpha_LAB = rot_CMS_to_LAB*alpha_CMS;
    tangent = alpha_LAB.Vect();
    energyPostBoost=alpha_LAB.E()-alpha_LAB.M();
    std::cout<<"ALPHA E = "<<alpha_LAB.E()-alpha_LAB.M()<<std::endl;   
  } else if(ion_id==pid_type::C_12 ){
    TLorentzVector carbon_CMS(tangent*p_carbon_CMS,Etot_carbon_CMS);
    TLorentzVector carbon_LAB = rot_CMS_to_LAB*carbon_CMS;
    tangent = carbon_LAB.Vect();
    energyPostBoost=carbon_LAB.E()-carbon_LAB.M();
    std::cout<<"CARBON E = "<<carbon_LAB.E()-carbon_LAB.M()<<std::endl;   

  }
  tangent.SetMag(1.0);
  //////////
  //TLorentzVector p4_beam = rot_CMS_to_LAB*p4_beam_CMS;
  //TLorentzVector p4_target = rot_CMS_to_LAB*p4_target_CMS;
  //TLorentzVector p4_total = p4_beam+p4_target;
  //p4_total.Print();
 // std::cout<<"^^^^^^^^^^^^^^^^^"<<std::endl;
  //////////
  double rangePostBoost=0;

  if(ion_id==pid_type::ALPHA)rangePostBoost = myRangeCalculator.getIonRangeMM(pid_type::ALPHA,energyPostBoost);
  if(ion_id==pid_type::C_12) rangePostBoost  = myRangeCalculator.getIonRangeMM(pid_type::CARBON_12, energyPostBoost);
  std::cout<<"----------------"<<std::endl;
  std::cout<<"rangeB4Boost -> "<<length<<std::endl;
  std::cout<<"rangePostBoost -> "<<rangePostBoost<<std::endl;
  std::cout<<"----------------"<<std::endl;
  

  return rangePostBoost;
}

//!------------------------------------------------------

/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////


//!------------------------------------------------------
TrackSegment3D EventSourceMC::createLorentzSegment(const TVector3 vertexPos, pid_type ion_id,TVector3 &aTangentNonBoost) const{
  
  double length = 0.0;
  double minLength = 0.0, maxLength = 0.0;
  double theta = 0.0, phi = 0.0;
  double minCosTheta = -1, maxCosTheta = 1;
  double minPhi = -TMath::Pi(), maxPhi = TMath::Pi();

  

  if(ion_id==pid_type::ALPHA){
    minLength = 20;
    maxLength = 60;    
    theta = TMath::ACos(myRndm.Uniform(minCosTheta, maxCosTheta));
    phi = myRndm.Uniform(minPhi, maxPhi);
  } 
  else if(ion_id==pid_type::C_12 && myTracks3D.size()==1){
    minLength = 10;
    maxLength = 15;
    theta = aTangentNonBoost.Theta();
    phi = aTangentNonBoost.Phi();
  }
  
  length = myRndm.Uniform(minLength, maxLength);
  std::cout<<"theta -> "<<theta<<std::endl;
  std::cout<<"phi -> "<<phi<<std::endl;
  std::cout<<"length -> "<<length<<std::endl;
 
  TVector3 tangent;
  tangent.SetMagThetaPhi(1.0, theta, phi);
  if(ion_id==pid_type::ALPHA){aTangentNonBoost=-tangent;}


  std::cout<<"====================="<<std::endl;
  std::cout<<"Pre Boost -> ";
  tangent.Print();
 double newLength =  Boost2Lab(length,tangent,ion_id);
  std::cout<<"Post Boost -> ";
  tangent.Print();
  std::cout<<"====================="<<std::endl;



  TrackSegment3D aSegment;
  aSegment.setGeometry(myGeometryPtr);
  aSegment.setStartEnd(vertexPos, vertexPos + tangent*newLength);
  aSegment.setPID(ion_id);
  
  return aSegment;
}
//!------------------------------------------------------

/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////

TH1F EventSourceMC::createChargeProfile(double ion_range, pid_type ion_id) const{
  
  double dEdx_max_energy =0.0, graph_range = 0.0;
  TGraph *dEdx_graph=0;
  
  if(ion_id==pid_type::ALPHA){
    dEdx_max_energy = braggGraph_alpha_energy;
    dEdx_graph=braggGraph_alpha;
    graph_range = 299;
  }
    else if(ion_id==pid_type::CARBON_12){
      dEdx_max_energy = braggGraph_12C_energy;
      dEdx_graph = braggGraph_12C;
      graph_range = 23.5;
    }
  
  double max_ion_range = myRangeCalculator.getIonRangeMM(ion_id, dEdx_max_energy);
  double delta = max_ion_range - ion_range;

  double x = 0.0;
  double value = 0.0;
  TH1F aChargeProfile{"hChargeProfile",";x [mm];dE/dx [keV/bin width]", 1024, -0.2*ion_range, 1.2*ion_range};
  for(int iBinX=1;iBinX<=aChargeProfile.GetNbinsX();++iBinX){
    x = aChargeProfile.GetBinCenter(iBinX);
    if(x<0 || (x+delta>graph_range)) continue;
    value = dEdx_graph->Eval(x+delta)*aChargeProfile.GetBinWidth(iBinX);
    aChargeProfile.SetBinContent(iBinX, value);
  }
  
  std::cout<<KBLU<<"Charge profile sum [keV]: \t"<<RST<<aChargeProfile.Integral()<<std::endl;
  std::cout<<KBLU<<"Ion energy from range calculator [keV]: \t"<<RST<<myRangeCalculator.getIonEnergyMeV(ion_id, ion_range)*1E3<<std::endl;
  return aChargeProfile;
}

/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////

//!------------------------------------------------------

Track3D EventSourceMC::createTrack(const TVector3 & aVtx, pid_type ion_id,TVector3 &aTangentNonBoost) const{
//Track3D EventSourceMC::createTrack(const TVector3 & aVtx,const TVector3 boost, pid_type ion_id) const{

  Track3D aTrack;
  TrackSegment3D aSegment = createLorentzSegment(aVtx, ion_id,aTangentNonBoost);
  aSegment.setVertexPos(aVtx);
  TH1F hChargeProfile = createChargeProfile(aSegment.getLength(), ion_id);
  aTrack.addSegment(aSegment);
  aTrack.setChargeProfile(hChargeProfile);
  return aTrack;
}
//!------------------------------------------------------

/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
void EventSourceMC::fill3DChargeCloud(const Track3D & aTrack){

  my3DChargeCloud.Reset();  
  TVector3 depositPosition, smearedPosition;
  TH1F hChargeProfile = aTrack.getChargeProfile();
  double lambda = 0.0, value = 0.0;
  double sigma = 2.0;
  int nTries = 100;
  
  for(int iBin=0;iBin<hChargeProfile.GetNbinsX();++iBin){
    value = hChargeProfile.GetBinContent(iBin);
    if(!value) continue;
    lambda = hChargeProfile.GetBinCenter(iBin);
    depositPosition = aTrack.getSegments().front().getStart() + lambda*aTrack.getSegments().front().getTangent();

    for(int iTry=0;iTry<nTries;++iTry){
      smearedPosition = TVector3(myRndm.Gaus(depositPosition.X(), sigma),
				 myRndm.Gaus(depositPosition.Y(), sigma),
				 myRndm.Gaus(depositPosition.Z(), sigma));
      my3DChargeCloud.Fill(smearedPosition.X(), smearedPosition.Y(), smearedPosition.Z(), value/nTries*keVToChargeScale);
    }
  }
  std::cout<<KBLU<<"Total charge cloud energy [100*keV]: "<<RST<<my3DChargeCloud.Integral()<<std::endl;
}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
void EventSourceMC::fillPEventTPC(const TH3D & h3DChargeCloud, const Track3D & aTrack){

  /*
  myProjectorPtr->SetEvent3D(my3DChargeCloud);
  myProjectorPtr->fillPEventTPC(myCurrentPEvent);
  */

  int iPolyBin = 0;
  double value = 0.0, totalCharge = 0.0;
  bool err_flag = false;  
  double sigma = 2.0;
  int nTries = 100;
  
  double lambda = 0.0;

  int iCell = 0;
  TVector3 depositPosition, smearedPosition;
  TH1F hChargeProfile = aTrack.getChargeProfile();

  for(int iBin=0;iBin<hChargeProfile.GetNbinsX();++iBin){
    value = hChargeProfile.GetBinContent(iBin);
    if(!value) continue;
    lambda = hChargeProfile.GetBinCenter(iBin);
    depositPosition = aTrack.getSegments().front().getStart() + lambda*aTrack.getSegments().front().getTangent();
    for(int iTry=0;iTry<nTries;++iTry){
      smearedPosition = TVector3(myRndm.Gaus(depositPosition.X(), sigma),
				 myRndm.Gaus(depositPosition.Y(), sigma),
				 myRndm.Gaus(depositPosition.Z(), sigma));
      iPolyBin = myGeometryPtr->GetTH2Poly()->FindBin(smearedPosition.X(), smearedPosition.Y());
      iCell = myGeometryPtr->Pos2timecell(smearedPosition.Z(), err_flag);
      std::shared_ptr<StripTPC> aStrip = myGeometryPtr->GetTH2PolyStrip(iPolyBin);
      if(aStrip){
	  myCurrentPEvent->AddValByStrip(aStrip, iCell, value/nTries*keVToChargeScale);
	  totalCharge+=value/nTries*keVToChargeScale;
      }
    }
  }
  std::cout<<KBLU<<"Total charge generated [100*keV]: "<<RST<<totalCharge<<std::endl;
}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
void EventSourceMC::generateSingleProng(pid_type ion_id){

  myGenEventType = ion_id;

  TVector3 aVtx = createVertex();
  myTracks3D.push_back(createTrack(aVtx, ion_id));
  std::cout<<KBLU<<"Generated track: "<<RST<<std::endl;
  std::cout<<myTracks3D.back()<<std::endl;
}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////


//!------------------------------------------------------
void EventSourceMC::generateTwoProng(){

  myGenEventType = CARBON_12;


  pid_type ion_id_first = pid_type::ALPHA;
  pid_type ion_id_second = pid_type::CARBON_12;
  TVector3 aVtx = createVertex();
  
 // TLorentzVector LorzentBoost = getBoost();
  TVector3 aTangentNonBoost;

  myTracks3D.push_back(createTrack(aVtx, ion_id_first,aTangentNonBoost));
  std::cout<<KBLU<<"First generated track: "<<RST<<std::endl;
  std::cout<<myTracks3D.back()<<std::endl;

  myTracks3D.push_back(createTrack(aVtx, ion_id_second,aTangentNonBoost));
  std::cout<<KBLU<<"Second generated track: "<<RST<<std::endl;
  std::cout<<myTracks3D.back()<<std::endl;


}
//!------------------------------------------------------

/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
void EventSourceMC::generateThreeProng(){

  myGenEventType = pid_type::THREE_ALPHA;

  pid_type ion_id = pid_type::ALPHA;
  TVector3 aVtx = createVertex();

  int nParts = 3;
  for(int iPart=0;iPart<nParts;++iPart){
    myTracks3D.push_back(createTrack(aVtx, ion_id));
    std::cout<<KBLU<<"Generated track number: "<<std::to_string(iPart)<<RST<<std::endl;
    std::cout<<myTracks3D.back()<<std::endl;
  }
}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
void EventSourceMC::generateEvent(){

  myCurrentEventInfo.SetEventId(myCurrentEntry);
  myCurrentEventInfo.SetRunId(0);
  myCurrentEventInfo.SetEventTimestamp(0);
  myCurrentEventInfo.SetPedestalSubtracted(true);
  
  myCurrentPEvent->Clear();
  myCurrentPEvent->SetEventInfo(myCurrentEventInfo);

  myTracks3D.clear();
  
  generateTwoProng();
  
  //double aRndm = myRndm.Uniform(0,1);
  //if(aRndm<0.33) generateSingleProng();
  //else if (aRndm<2*0.33) generateTwoProng();
  //else generateThreeProng();
  
  for(const auto & aTrack: myTracks3D) fillPEventTPC(my3DChargeCloud, aTrack);
  fillEventTPC();
  ++myCurrentEntry;
}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////

Track3D EventSourceMC::createTrack(const TVector3 & aVtx, pid_type ion_id) const{
//Track3D EventSourceMC::createTrack(const TVector3 & aVtx,const TVector3 boost, pid_type ion_id) const{

  Track3D aTrack;
  TrackSegment3D aSegment = createSegment(aVtx, ion_id);
  aSegment.setVertexPos(aVtx);
  TH1F hChargeProfile = createChargeProfile(aSegment.getLength(), ion_id);
  aTrack.addSegment(aSegment);
  aTrack.setChargeProfile(hChargeProfile);
  return aTrack;
}





double EventSourceMC::SRIMtrackLength(const double gammaEnergy, const double pressure, pid_type ion_id)const{
  IonRangeCalculator myRangeCalculator(gas_mixture_type::CO2,pressure,293.15);
  double energy=0; //calculate the energy given 
  double range=myRangeCalculator.getIonRangeMM(ion_id, energy);
  std::cout<<range<<std::endl;
  /*
  if(ion_id)std:cout<<pressure<<<std:endl;
  
  const double AMU = 931.49410242; // MeV/c^2
  const double alphaMass  =  4.00260325413  *AMU ;
  const double oxygenMass = 15.99491461956*AMU; 
  const double carbonMass = 12.000*AMU; 


  // Pre split
  TVector3 beamDir_LAB(-1,0, 0); // unit vector
  TLorentzVector p4_beam_CMS(beamDir_LAB.Unit()*gammaEnergy, gammaEnergy); // [MeV/c, MeV/c, MeV/c, MeV/c^2]
  TLorentzVector p4_target_CMS(0, 0, 0, oxygenMass); // [MeV/c, MeV/c, MeV/c, MeV/c^2]
  TLorentzVector p4_total_CMS =  p4_beam_CMS + p4_target_CMS;

  //T
  Vector3 tangent_alphaMin;
  tangent.SetMagThetaPhi(1.0, 180, 0);
  //TODO
  double totalEnergy_CMS=(p4_target_CMS+p4_beam_CMS).E();
  double totalEnergy_CMS_check=sqrt( oxygenMass*(2*beamEnergy+oxygenMass) ); 
  
  std::cout<<"totalE: "<<totalEnergy_CMS<<" <> "<<totalEnergy_CMS_check<<std::endl;
  */
  
  return 0.;
}