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
  double xmin=-165+5;
  double xmax=165+5;
  std::uniform_real_distribution<> dis(xmin, xmax);
  std::normal_distribution<> disGaus(0, 3);
  
  x = dis(gen1); 
  y = disGaus(gen2);
  z = disGaus(gen3);
  //std::cout<<"-> "<<x<<" "<<y<<" "<<z<<std::endl;



  TVector3 aVertex(x,y,z);
  return aVertex;
}
//!------------------------------------------------------

/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////

//!------------------------------------------------------
double EventSourceMC::Boost2Lab(double const gammaEnergy,double const length,TVector3 &tangent, pid_type ion_id, TLorentzVector &output)const{
  bool checkSwitch= false;

  // ** KE ** //
  IonRangeCalculator myRangeCalculator(gas_mixture_type::CO2,250,293.15);
  double alphaEnergyKe = myRangeCalculator.getIonEnergyMeV(pid_type::ALPHA,length);
  double carbonEnergyKe =myRangeCalculator.getIonEnergyMeV(pid_type::CARBON_12, length);
 
  // ** ION INFO ** //
  const double AMU        = 931.49410242; // MeV/c^2
  const double alphaMass  = 4.00260325413 * AMU ;
  const double oxygenMass = 15.99491461956 * AMU; 
  const double carbonMass = 12.000 * AMU; 

  // ** BEAM/TARGE LAB ** //
  TVector3 beamDir_LAB(-1,0, 0); // unit vector
  double beamEnergy_LAB = gammaEnergy;
  TLorentzVector p4_beam_LAB(beamDir_LAB.Unit()*beamEnergy_LAB, beamEnergy_LAB); // [MeV/c, MeV/c, MeV/c, MeV/c^2]
  TLorentzVector p4_target_LAB(0, 0, 0, oxygenMass); // [MeV/c, MeV/c, MeV/c, MeV/c^2]
  TLorentzVector p4_total_LAB =  p4_beam_LAB + p4_target_LAB;

  // ** DEFINE BOOSTS **//
  TVector3 boost_CMS_to_LAB = p4_total_LAB.BoostVector(); 
  TLorentzRotation rot_CMS_to_LAB(boost_CMS_to_LAB);

  TVector3 boost_LAB_to_CMS = p4_total_LAB.BoostVector(); 
  TLorentzRotation rot_LAB_to_CMS(-boost_LAB_to_CMS);

  // ** BEAM/TARGET CMS ** //
  TLorentzVector p4_beam_CMS   =  rot_LAB_to_CMS * p4_beam_LAB;
  TLorentzVector p4_target_CMS =  rot_LAB_to_CMS * p4_target_LAB;
  TLorentzVector p4_total_CMS  =  p4_beam_CMS + p4_target_CMS;

  // **  CHECKS ** //
  double p_alpha_CMS      = sqrt(alphaEnergyKe*(alphaEnergyKe+2*alphaMass));
  double p_carbon_CMS     = sqrt(carbonEnergyKe*(carbonEnergyKe+2*carbonMass));
  double Etot_alpha_CMS   = alphaMass+alphaEnergyKe;
  double Etot_carbon_CMS  = carbonMass+carbonEnergyKe;


  if(checkSwitch==true){
    if(ion_id==pid_type::C_12 )std::cout<<KGRN<<"__CARBON__"<<RST<<std::endl;
    if(ion_id==pid_type::ALPHA)std::cout<<KGRN<<"____ALPHA____"<<RST<<std::endl;
    std::cout<<KRED<<"********************************"<<RST<<std::endl;
    std::cout<<"Check that p4_total_LAB - > p4_total_CMS goes to zero under inverse boost"<<std::endl;  
    p4_total_LAB.Print();
    p4_total_CMS.Print();
    std::cout<<KRED<<"********************************"<<RST<<std::endl;
    std::cout<<"OXYGEN E in LAB: "<< p4_target_LAB.E()<<std::endl;
    std::cout<<"OXYGEN P in LAB: "<<p4_target_LAB.P()<<std::endl;
    std::cout<<"OXYGEN E in CMS: "<< p4_target_CMS.E()<<std::endl;
    std::cout<<"OXYGEN P in CMS: "<<p4_target_CMS.P()<<std::endl;
    std::cout<<KRED<<"********************************"<<RST<<std::endl;
    std::cout<<"GAMMA  E in LAB: "<<beamEnergy_LAB<<std::endl;
    std::cout<<"GAMMA  P in LAB: "<<p4_beam_LAB.P()<<std::endl;
    std::cout<<"GAMMA  E in CMS: "<< p4_beam_CMS.E()<<std::endl;
    std::cout<<"GAMMA  P in CMS: "<<p4_beam_CMS.P()<<std::endl;
    std::cout<<KRED<<"********************************"<<RST<<std::endl;
    std::cout<<"TOTAL  E in LAB: "<< p4_total_LAB.E()<<std::endl;
    std::cout<<"TOTAL  P in LAB: "<<p4_total_LAB.P()<<std::endl;
    std::cout<<KRED<<"********************************"<<RST<<std::endl;
    std::cout<<"TOTAL  E in CMS: "<<p4_total_CMS.E()<<std::endl;
    std::cout<<"TOTAL  P in CMS: "<<p4_total_CMS.P()<<std::endl;
    std::cout<<KRED<<"********************************"<<RST<<std::endl;
    std::cout<<"Alpha  KE in CMS: "<<alphaEnergyKe<<std::endl;
    std::cout<<"Alpha  E in CMS: "<<Etot_alpha_CMS<<std::endl;
    std::cout<<KRED<<"********************************"<<RST<<std::endl;
    std::cout<<"CARBON KE in CMS: "<<carbonEnergyKe<<std::endl;
    std::cout<<"CARBON E in CMS: "<<Etot_carbon_CMS<<std::endl;
    std::cout<<KRED<<"********************************"<<RST<<std::endl;
  }


  
  
  double energyPostBoost=0;
  if(ion_id==pid_type::ALPHA){

    // ** BOOST ** //
    TLorentzVector alpha_CMS(tangent * p_alpha_CMS, Etot_alpha_CMS);
    TLorentzVector alpha_LAB = rot_CMS_to_LAB * alpha_CMS;

    // ** THINGS TO SAVE ** //
    tangent = alpha_LAB.Vect();
    energyPostBoost=alpha_LAB.E()-alpha_LAB.M();
    output = alpha_LAB;


    if(checkSwitch ==true){
      alpha_LAB.Print();
      std::cout<<"ALPHA  KE in LAB: "<<alpha_LAB.E()-alpha_LAB.M()<<std::endl; 
      std::cout<<"ALPHA  KE in CMS (xcheck): "<<KBLU<<alpha_CMS.E()-alpha_CMS.M()<<RST<<std::endl; 
      std::cout<<"ALPHA  P in LAB: "<<alpha_LAB.P()<<std::endl;
      std::cout<<"ALPHA  P in CMS: "<<alpha_CMS.P()<<std::endl;
    }
  
  }else if(ion_id==pid_type::C_12 ){
    
    // ** BOOST ** //
    TLorentzVector carbon_CMS(tangent * p_carbon_CMS, Etot_carbon_CMS);
    TLorentzVector carbon_LAB = rot_CMS_to_LAB*carbon_CMS;

    output = carbon_LAB;
    carbon_LAB.Print();

    // ** THINGS TO SAVE ** //
    tangent = carbon_LAB.Vect();
    energyPostBoost=carbon_LAB.E()-carbon_LAB.M();

    if(checkSwitch ==true){
      std::cout<<"CARBON  KE in LAB : "<<carbon_LAB.E()-carbon_LAB.M()<<std::endl; 
      std::cout<<"CARBON  KE in CMS (xcheck): "<<KBLU<<carbon_CMS.E()-carbon_CMS.M()<<RST<<std::endl; 
      std::cout<<"CARBON  P in LAB: "<<carbon_LAB.P()<<std::endl;
      std::cout<<"CARBON  P in CMS: "<<carbon_CMS.P()<<std::endl;
    }
  }
  tangent.SetMag(1.0);

  double rangePostBoost=0;

  if(ion_id==pid_type::ALPHA)rangePostBoost  = myRangeCalculator.getIonRangeMM(pid_type::ALPHA,energyPostBoost);
  if(ion_id==pid_type::C_12 )rangePostBoost  = myRangeCalculator.getIonRangeMM(pid_type::CARBON_12, energyPostBoost);

  if(checkSwitch ==true){
    std::cout<<"----------------"<<std::endl;
    std::cout<<"range before Boost -> "<<length<<std::endl;
    std::cout<<"range after Boost -> "<<rangePostBoost<<std::endl;
    std::cout<<"----------------"<<std::endl;
  }

  return rangePostBoost;
}

//!------------------------------------------------------

/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////


//!------------------------------------------------------
TrackSegment3D EventSourceMC::createLorentzSegment(const TVector3 vertexPos, pid_type ion_id,TVector3 &aTangentNonBoost, TLorentzVector &p4outPut) const{
  
  double length = 0.0;
  double theta = 0.0, phi = 0.0;
  double minCosTheta = -1, maxCosTheta = 1;
  double minPhi = -TMath::Pi(), maxPhi = TMath::Pi();
  double gammaEnergy = 13.1;

  TLorentzVector p4outPut_alpha;



  if(ion_id==pid_type::ALPHA){
    length= SRIMNormtrackLength(gammaEnergy,250,ion_id); 
    theta = TMath::ACos(myRndm.Uniform(minCosTheta, maxCosTheta));
    phi = myRndm.Uniform(minPhi, maxPhi);
  } 
  else if(ion_id==pid_type::C_12 && myTracks3D.size()==1){
    p4outPut_alpha=p4outPut;
    length =     SRIMNormtrackLength(gammaEnergy,250,ion_id);
    theta = aTangentNonBoost.Theta();
    phi = aTangentNonBoost.Phi();
  }
  

  //std::cout<<"theta b4 boost-> "<<theta<<std::endl;
  //std::cout<<"phi b4 boost-> "<<phi<<std::endl;
  //std::cout<<"length b4 boost -> "<<length<<std::endl;
 
  TVector3 tangent;
  tangent.SetMagThetaPhi(1.0, theta, phi);
  if(ion_id==pid_type::ALPHA){aTangentNonBoost=-tangent;}
  TVector3 tangentCopy =tangent;

 
  //std::cout<<"====================="<<std::endl;
  double newLength =  Boost2Lab(gammaEnergy,length,tangent,ion_id,p4outPut);
  std::cout<<"Pre Boost -> ";
  tangentCopy.Print();
  std::cout<<"Post Boost -> ";
  tangent.Print();
  //std::cout<<"====================="<<std::endl;
  //p4outPut.Print();
  p4outPut.Print();
  TLorentzVector p4outPut_carbon;


 if(ion_id==pid_type::C_12){
    p4outPut_carbon = p4outPut;

    TLorentzVector totalP4 = p4outPut_alpha+p4outPut_carbon;
    //std::cout<<"=======££££=========="<<std::endl;
    //p4outPut_alpha.Print();
    //p4outPut_carbon.Print();
    std::cout<<"TOTAL P -> ";
    totalP4.Print();
    //std::cout<<totalP4.X() << std::endl;
    //std::cout<<totalP4.Y() << std::endl;
    //std::cout<<totalP4.Z() << std::endl;
    //std::cout<<"=======££££=========="<<std::endl;


    p4outPut = totalP4;

   }



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

Track3D EventSourceMC::createTrack(const TVector3 & aVtx, pid_type ion_id,TVector3 &aTangentNonBoost, TLorentzVector &p4outPut) const{
//Track3D EventSourceMC::createTrack(const TVector3 & aVtx,const TVector3 boost, pid_type ion_id) const{

  Track3D aTrack;
//  TLorentzVector p4outPut;
  TrackSegment3D aSegment = createLorentzSegment(aVtx, ion_id,aTangentNonBoost,p4outPut);
  aSegment.setVertexPos(aVtx);
   
   /*
  if(ion_id==pid_type::C_12){
    std::cout <<"HERE"<<p4outPut.X() <<std::endl;
    std::cout <<"HERE"<<p4outPut.Y() <<std::endl;
    std::cout <<"HERE"<<p4outPut.Z() <<std::endl;
    aSegment.setMomentumVec(p4outPut);}
    */
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
  TLorentzVector p4outPut;


  myTracks3D.push_back(createTrack(aVtx, ion_id_first,aTangentNonBoost,p4outPut));
  std::cout<<KBLU<<"First generated track: "<<RST<<std::endl;
  std::cout<<myTracks3D.back()<<std::endl;

  myTracks3D.push_back(createTrack(aVtx, ion_id_second,aTangentNonBoost,p4outPut));
  std::cout<<KBLU<<"Second generated track: "<<RST<<std::endl;
  std::cout<<myTracks3D.back()<<std::endl;

  myPvec = p4outPut;


 // return p4outPut;
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
  //generateThreeProng();
  
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

double EventSourceMC::SRIMNormtrackLength(const double gammaEnergy, const double pressure, pid_type ion_id)const{
  IonRangeCalculator myRangeCalculator(gas_mixture_type::CO2,pressure,293.15);
  const double AMU        = 931.49410242; // MeV/c^2
  const double alphaMass  = 4.00260325413 * AMU ;
  const double carbonMass = 12.000 * AMU; 
  double ke = gammaEnergy - 7.162; // -Q
  if(ion_id==pid_type::ALPHA){
    ke = ke*(carbonMass/(alphaMass+carbonMass));

  } else if(ion_id==pid_type::C_12){
    ke = ke * (alphaMass/(alphaMass+carbonMass));

  }



  double length = myRangeCalculator.getIonRangeMM(ion_id, ke);  

  return length;





} 


void EventSourceMC::SRIMtrackLength(const double gammaEnergy, const double pressure, pid_type ion_id, double &minLen ,double &maxLen)const{
  
  IonRangeCalculator myRangeCalculator(gas_mixture_type::CO2,pressure,293.15);
  
  
  
  const double AMU        = 931.49410242; // MeV/c^2
  const double alphaMass  = 4.00260325413 * AMU ;
  const double oxygenMass = 15.99491461956 * AMU; 
  const double carbonMass = 12.000 * AMU; 

  // BOOST: CMS->LAB
  TVector3 beamDir_LAB(-1,0, 0); // unit vector
  TLorentzVector p4_beam_LAB(beamDir_LAB.Unit()*gammaEnergy, gammaEnergy); // [MeV/c, MeV/c, MeV/c, MeV/c^2]
  TLorentzVector p4_target_LAB(0, 0, 0, oxygenMass); // [MeV/c, MeV/c, MeV/c, MeV/c^2]
  TLorentzVector p4_total_LAB =  p4_target_LAB + p4_beam_LAB;
  TVector3 boost_CMS_to_LAB = p4_total_LAB.BoostVector(); 
  TLorentzRotation rot(-boost_CMS_to_LAB);


  TVector3 tangent1; // max
  tangent1.SetMagThetaPhi(1.0, 360, 0);
 
  TVector3 tangent2; // min
  tangent2.SetMagThetaPhi(1.0, 180, 0);
  
  double p = 0;
  double ke = gammaEnergy - 7.162;
  double addMass=0;


  if(ion_id==pid_type::ALPHA){
    ke = ke*(carbonMass/(alphaMass+carbonMass));
    p = sqrt(ke*(ke+2*alphaMass));
     addMass= alphaMass;

  } else if(ion_id==pid_type::C_12){
    ke = ke * (alphaMass/(alphaMass+carbonMass));
    p = sqrt(ke*(ke+2*carbonMass));
     addMass= carbonMass;

  }




  TLorentzVector boost1(tangent1 * p, ke+addMass);
  TLorentzVector boosted1 = rot * boost1;
  
  
  TLorentzVector boost2(tangent2 * p, ke+addMass);
  TLorentzVector boosted2 = rot * boost2;


  double maxE = boosted1.E()-boosted1.M();
  double minE = boosted2.E()-boosted2.M();

  
  maxLen = myRangeCalculator.getIonRangeMM(ion_id, maxE);  
  minLen = myRangeCalculator.getIonRangeMM(ion_id, minE);

  


  return;


}