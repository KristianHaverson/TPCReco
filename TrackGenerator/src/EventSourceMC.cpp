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
  double minCosTheta = 0.9, maxCosTheta = 0.9;
  double minPhi = 0, maxPhi = 2*M_PI;
  double minLength = 20, maxLength = 20;



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
    minCosTheta = cos(- aTangent.Theta());//-0.1;
    maxCosTheta = cos(- aTangent.Theta());//+0.1;
    //if(minCosTheta<-1) minCosTheta = -1.0;
    //if(maxCosTheta>1) maxCosTheta = 1.0;
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
  //tangent.Print(0);

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

 // (-45.5864, 5.81417, -59.9457)

  // x = -45.5864; 
  // y = 5.81417;
  // z = -59.9457;
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
  const double minCosTheta = -1.0;
	const double maxCosTheta = 1.0;
  double minPhi =-TMath::Pi(), maxPhi = TMath::Pi();
  double gammaEnergy = 13.1;

  TLorentzVector p4outPut_alpha;



  if(ion_id==pid_type::ALPHA){
    length= SRIMNormtrackLength(gammaEnergy,250,ion_id); 
    theta = TMath::ACos(myRndm.Uniform(minCosTheta, maxCosTheta));
    phi = myRndm.Uniform(minPhi, maxPhi);

    //length = 79.9807;
    //theta=0.761485;
    //phi = 2.48371;


  } 
  else if(ion_id==pid_type::C_12 && myTracks3D.size()==1){
    p4outPut_alpha=p4outPut;
    length =     SRIMNormtrackLength(gammaEnergy,250,ion_id);
    theta = aTangentNonBoost.Theta();
    phi = aTangentNonBoost.Phi();
    
  //  length = 14.5987;
  //  theta  = 2.4118;
  //  phi = -0.734923;

  }
  

//  std::cout<<"theta b4 boost-> "<<theta*180/TMath::Pi()<<std::endl;
//  std::cout<<"phi b4 boost-> "<<phi*180/TMath::Pi()<<std::endl;
//  std::cout<<"length b4 boost -> "<<length<<std::endl;
 
  TVector3 tangent;
  tangent.SetMagThetaPhi(1.0, theta, phi);
  if(ion_id==pid_type::ALPHA){aTangentNonBoost=-tangent;}
  TVector3 tangentCopy =tangent;

 
  //std::cout<<"====================="<<std::endl;
  double newLength = length;//Boost2Lab(gammaEnergy,length,tangent,ion_id,p4outPut);
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


/////////////////////////////////////////////////////////

double EventSourceMC::BoostbyVel(double const gammaEnergy,AlphaData &alphaData,const int index,TVector3 &tangent)const{


  std::cout<<gammaEnergy<<std::endl;
  IonRangeCalculator myRangeCalculator(gas_mixture_type::CO2,250,293.15);

  const double AMU        = 931.49410242; // MeV/c^2
  const double alphaMass  = 4.00260325413 * AMU ;

  // ** BEAM/TARGE LAB ** //
  TVector3 beamDir_LAB = -alphaData.getTangent(0);
  double beamEnergy_LAB =  2*alphaData.getAlphaEnergy(1);
//  std::cout<<KRED<<beamEnergy_LAB<<RST<<std::endl;


  TLorentzVector p4_beam_LAB(beamDir_LAB.Unit()*beamEnergy_LAB, beamEnergy_LAB); // [MeV/c, MeV/c, MeV/c, MeV/c^2]
  TLorentzVector p4_target1_LAB(0, 0, 0, 2*alphaMass); // [MeV/c, MeV/c, MeV/c, MeV/c^2]
  TLorentzVector p4_total_LAB =  p4_beam_LAB+ p4_target1_LAB;//+p4_target2_LAB;
  
  // ** DEFINE BOOSTS **//
  TVector3 boost_CMS_to_LAB = p4_total_LAB.BoostVector(); 
  TLorentzRotation rot_CMS_to_LAB(boost_CMS_to_LAB);

  TVector3 boost_LAB_to_CMS = p4_total_LAB.BoostVector(); 
  TLorentzRotation rot_LAB_to_CMS(-boost_LAB_to_CMS);

//  std::cout<<KRED<<beamEnergy_LAB<<RST<<std::endl;
  // ** BEAM/TARGET CMS ** //
  TLorentzVector p4_beam_CMS   =  rot_LAB_to_CMS * p4_beam_LAB;
  TLorentzVector p4_target1_CMS =  rot_LAB_to_CMS * p4_target1_LAB;
  TLorentzVector p4_total_CMS  =  p4_beam_CMS + p4_target1_CMS;//+p4_target2_CMS;

    // **  CHECKS ** //
  double alphaEnergyKe = alphaData.getAlphaEnergy(index);
  double p_alpha1_CMS      = sqrt(alphaEnergyKe*(alphaEnergyKe+2*alphaMass));
  double p_alpha2_CMS      = sqrt(alphaEnergyKe*(alphaEnergyKe+2*alphaMass));
  double Etot_alpha1_CMS   = alphaMass+alphaEnergyKe;
  double Etot_alpha2_CMS   = alphaMass+alphaEnergyKe;


  std::cout<<KRED<<"********************************"<<RST<<std::endl;
  std::cout<<"Check that p4_total_LAB - > p4_total_CMS goes to zero under inverse boost"<<std::endl;  
  p4_total_LAB.Print();
  p4_total_CMS.Print();
  std::cout<<KRED<<"********************************"<<RST<<std::endl;
  std::cout<<"TOTAL  E in LAB: "<< p4_total_LAB.E()<<std::endl;
  std::cout<<"TOTAL  P in LAB: "<<p4_total_LAB.P()<<std::endl;
  std::cout<<KRED<<"********************************"<<RST<<std::endl;
  std::cout<<"TOTAL  E in CMS: "<<p4_total_CMS.E()<<std::endl;
  std::cout<<"TOTAL  P in CMS: "<<p4_total_CMS.P()<<std::endl;
  std::cout<<KRED<<"********************************"<<RST<<std::endl;
  std::cout<<"Alpha1  KE in CMS: "<<alphaEnergyKe<<std::endl;
  std::cout<<"Alpha1  E in CMS: "<<Etot_alpha1_CMS<<std::endl;
  std::cout<<"Alpha1  P in CMS: "<<p_alpha1_CMS<<std::endl;
  std::cout<<KRED<<"********************************"<<RST<<std::endl;
  std::cout<<"Alpha2  KE in CMS: "<<alphaEnergyKe<<std::endl;
  std::cout<<"Alpha2  E in CMS: "<<Etot_alpha2_CMS<<std::endl;
  std::cout<<"Alpha1  P in CMS: "<<p_alpha2_CMS<<std::endl;
  std::cout<<KRED<<"********************************"<<RST<<std::endl;

  //double p_alpha_CMS      = sqrt(alphaData.getAlphaEnergy(index)*(alphaData.getAlphaEnergy(index)+2*alphaMass));
  //double Etot_alpha_CMS   = alphaMass+alphaData.getAlphaEnergy(index);

  //TVector3 tangent =  alphaData.getTangent(index);
  TLorentzVector alpha_CMS(tangent * p_alpha1_CMS, Etot_alpha1_CMS);
  TLorentzVector alpha_LAB = rot_CMS_to_LAB * alpha_CMS;
  alpha_LAB.Print();
  std::cout<<"ALPHA  KE in LAB: "<<alpha_LAB.E()-alpha_LAB.M()<<std::endl; 
  std::cout<<"ALPHA  M in LAB: "<<alpha_LAB.M()<<std::endl; 
  std::cout<<"ALPHA  KE in CMS (xcheck): "<<KBLU<<alpha_CMS.E()-alpha_CMS.M()<<RST<<std::endl; 
  std::cout<<"ALPHA  M in CMS: "<<KBLU<<alpha_CMS.M()<<RST<<std::endl; 
  std::cout<<"ALPHA  P in LAB: "<<alpha_LAB.P()<<std::endl;
  std::cout<<"ALPHA  P in CMS: "<<alpha_CMS.P()<<std::endl;

  tangent.SetMag(1.0);

  alphaData.setTangent(index,alpha_LAB.Vect());
  double energyPostBoost=alpha_LAB.E()-alpha_LAB.M();
  double rangePostBoost  = myRangeCalculator.getIonRangeMM(pid_type::ALPHA,energyPostBoost);
  std::cout<<"energy pre -> "<<alphaEnergyKe<<std::endl;
  std::cout<<"energy post -> "<<energyPostBoost<<std::endl;
  std::cout<<"range before Boost -> "<<alphaData.getAlphaLength(index)<<std::endl;
  std::cout<<"range after Boost -> "<<rangePostBoost<<std::endl;

  alphaData.setAlphaLength(index,rangePostBoost);
  alphaData.setAlphaEnergy(index,energyPostBoost);



  
return 0;
}




TrackSegment3D EventSourceMC::create3pSegment(const TVector3 vertexPos, pid_type ion_id,AlphaData &alphaData, const int index, const int groun_or_exc) const{
  

  double gammaE = 13.1;

  TVector3 tangent;
  double length=0;

 


  // angles etc in COM
  if((ion_id==pid_type::ALPHA) && (index==0)){
    // Sets inital lengths & energies
    SRIMNormtrackLength(gammaE,250,ion_id,alphaData, 0); 
    SRIMNormtrackLength(gammaE,250,ion_id,alphaData, 1); 
    SRIMNormtrackLength(gammaE,250,ion_id,alphaData, 2); 
  }
  if(index==0){
    
	double Kev_to_MeV = 1./1000.;

	//========================================================================//
	// ** MASSES ** //
	const double AMU        = 931.49410242;
	const double mass_c12   = AMU  *12.00000000000; 
	const double mass_8be   = AMU  * 8.00530510000; 
	const double mass_alpha = AMU  * 4.00260325413; 

	//========================================================================//
	// ** ENERGY ** //
	const double Etot_gamma =gammaE;

	const double Q_VALUE_12C  = - 7367.0 * Kev_to_MeV; 
	const double Q_VALUE_8BE  = 92.0   * Kev_to_MeV; 

	// ** KE ** //
	double ECM1 = Etot_gamma + Q_VALUE_12C; 


  double exc_min = 0.001;
  double exc_xmax = ECM1-0.01;
  //std::cout<<exc_min<<std::endl;
  //std::cout<<exc_min<<std::endl;
  //std::cout<<exc_xmax<<std::endl;
  //std::cout<<exc_xmax<<std::endl;
  std::uniform_real_distribution<> dis(exc_min,exc_xmax);

  std::random_device rd; // Seed random number
  std::mt19937 gen1(rd());

  const double exct = dis(gen1);
//  std::cout<<"<>HERE<>"<<std::endl;
//  std::cout<<exct<<std::endl;
//  std::cout<<exct<<std::endl;
	
  if(groun_or_exc==1)ECM1=ECM1-exct;

  const double KE_alpha0  = ECM1 * (mass_8be  /(mass_alpha+mass_8be));
	const double KE_8be     = ECM1 * (mass_alpha/(mass_alpha+mass_8be));

	double ECM2 = Q_VALUE_8BE; 
  if(groun_or_exc==1)ECM2=ECM2+exct;
	const double KE_alpha12 = ECM2 * (mass_alpha/(mass_alpha+mass_alpha));

	// ** TOTAL ENERGY ** //
	const double Etot_alpha0  = mass_alpha + KE_alpha0;
	const double Etot_8be	  = mass_8be   + KE_8be;
	const double Etot_alpha12 = mass_alpha + KE_alpha12;

	//========================================================================//
	// ** ANGLES ** //
	double theta = 0.0; 
	double phi = 0.0;
 //double length;
	const double minCosTheta = -1.0;
	const double maxCosTheta = 1.0;
	const double minPhi = -TMath::Pi();
	const double maxPhi = TMath::Pi();

	//========================================================================//
	// ** MOMENTA ** //
	const double p_alpha0   = sqrt(KE_alpha0 *(KE_alpha0 +2*mass_alpha));
	const double p_8be	    = sqrt(KE_8be    *(KE_8be    +2*mass_8be));
	const double p_alpha12  = sqrt(KE_alpha12*(KE_alpha12+2*mass_alpha));


	//========================================================================//
	// ** CHECK KE ** //
	std::cout<<"*******************************************************"<<std::endl;   
	std::cout<<" -- CHECK TOTAL KE -- "<<std::endl;
	std::cout<<"    ECM1(abCM)-> " <<Etot_gamma<<" + "<<Q_VALUE_12C <<"   = "<<ECM1 <<std::endl;
	std::cout<<"    Total Ke  -> 8be             = "<<KE_8be<<std::endl; 
	std::cout<<"    Total Ke  -> a0              = "<<KE_alpha0<<std::endl; 
	std::cout<<"    Total Ke  -> 8be + a0        = "<<KE_8be+KE_alpha0<<std::endl;
	std::cout<<"*******************************************************"<<std::endl;   
	std::cout<<" -- CHECK TOTAL KE -- "<<std::endl;
	std::cout<<"    ECM2(bCM) -> " <<ECM2 <<std::endl;
	std::cout<<"    Total Ke  -> a1 + a2         = "<<KE_alpha12+KE_alpha12<<std::endl; 
	std::cout<<"*******************************************************"<<std::endl;   



	//========================================================================//
	// ** Alpha-0 and (8)Berlyllium frame of reference  -> abCM ** //

	// Alpha-0 Vector
	TVector3 tangent_alpha0_abCM;
	theta = TMath::ACos(myRndm.Uniform(minCosTheta, maxCosTheta));
	phi   = myRndm.Uniform(minPhi, maxPhi);
	tangent_alpha0_abCM.SetMagThetaPhi(1.0, theta, phi);
	TLorentzVector p4_alpha0_abCM(tangent_alpha0_abCM*p_alpha0 ,Etot_alpha0);

	// (8)Berlyllium Vector (B2B in abCM)
	TVector3 tangent_be8_abCM = -tangent_alpha0_abCM;
	TLorentzVector p4_be8_abCM(tangent_be8_abCM*p_8be , Etot_8be);


	TLorentzVector p4_total_abCM =  p4_alpha0_abCM + p4_be8_abCM;

	TVector3 boost_total_abCM  = p4_total_abCM.BoostVector(); 
	TLorentzRotation rot_to_abCM(-boost_total_abCM);

	p4_be8_abCM    = rot_to_abCM * p4_be8_abCM;
	p4_alpha0_abCM = rot_to_abCM * p4_alpha0_abCM;
	p4_total_abCM  = p4_be8_abCM + p4_alpha0_abCM;

	std::cout<<" -- abCM Frame -- "<<std::endl;
	std::cout<<"    a0 + 8be in abCM p4 ->";p4_total_abCM.Print();
	std::cout<<"    KE of a0     = "<<p4_alpha0_abCM.E()-p4_alpha0_abCM.M()<<std::endl;
	std::cout<<"*******************************************************"<<std::endl;   



	//========================================================================//
	// ** Alpha-1 and Alpha-2 in (8)Berlyllium frame of reference -> bCM ** //

	// Alpha-1 Vector
	TVector3 tangent_alpha1_bCM;
	theta = TMath::ACos(myRndm.Uniform(minCosTheta, maxCosTheta));
	phi   = myRndm.Uniform(minPhi, maxPhi);
	tangent_alpha1_bCM.SetMagThetaPhi(1.0, theta, phi);
	TLorentzVector p4_alpha1_bCM(tangent_alpha1_bCM*p_alpha12 ,Etot_alpha12);

	// Alpha-2 Vector (B2B in bCM)
	TVector3 tangent_alpha2_bCM =   -tangent_alpha1_bCM;
	TLorentzVector p4_alpha2_bCM(tangent_alpha2_bCM*p_alpha12 ,Etot_alpha12);

	//* total four momentum of alpha1, and alpah2 in bCM = beryllium energy in its CM)
	TLorentzVector p4_a12_bCM = p4_alpha1_bCM+ p4_alpha2_bCM;



	std::cout<<" -- bCM Frame -- "<<std::endl;
	std::cout<<"    a1 + a2 in bCM p4 ->";p4_a12_bCM.Print();
	std::cout<<"    KE of a1+a2  = "<<(p4_alpha1_bCM.E()-p4_alpha1_bCM.M())+(p4_alpha2_bCM.E()-p4_alpha2_bCM.M()) <<std::endl;


	//========================================================================//
	// ** Boost Alpha-1 & Alpha-2 from bCM -> abCM ** //

	TVector3 boost_8be = p4_be8_abCM.BoostVector(); 
	TLorentzRotation rot_to_abCM_by_8be(boost_8be);
	//=== TEST === //
	//* total four momentum of alpha2, and alpah3 in beryllium CM is correct (= beryllium energy in its CM)
	TLorentzRotation rot_from_abCM_by_8be(-boost_8be);
	TLorentzVector be_in_bCM   =  rot_from_abCM_by_8be * p4_be8_abCM;

	std::cout<<"    Energy 8Be in bCM    = "<<be_in_bCM.E()<<std::endl;
	std::cout<<"    Energy a1+a2 in bCM  = "<<p4_a12_bCM.E()<<std::endl;
	std::cout<<"*******************************************************"<<std::endl;   

	//==========================================


	// ** alpha1 and alpha2 in abCMS ** //
	TLorentzVector p4_alpha1_abCM   =  rot_to_abCM_by_8be * p4_alpha1_bCM;
	TLorentzVector p4_alpha2_abCM 	=  rot_to_abCM_by_8be * p4_alpha2_bCM;

	//========================================================================//
	// ** Boost abCM -> LAB ** //

	// ** BEAM/TARGET LAB ** //
	TVector3 beamDir_LAB(-1,0, 0); // unit vector
	double beamEnergy_LAB = Etot_gamma;
	TLorentzVector p4_beam_LAB(beamDir_LAB.Unit()*beamEnergy_LAB, beamEnergy_LAB); // [MeV/c, MeV/c, MeV/c, MeV/c^2]
	TLorentzVector p4_target_LAB(0, 0, 0, mass_c12); // [MeV/c, MeV/c, MeV/c, MeV/c^2]
	TLorentzVector p4_total_LAB =  p4_beam_LAB + p4_target_LAB;


	// ** DEFINE BOOSTS **//
	TVector3 boost_abCM_to_LAB = p4_total_LAB.BoostVector(); 
	TLorentzRotation rot_abCM_to_LAB(boost_abCM_to_LAB);

	// ** BOOST 2 LAB **//
	TLorentzVector p4_alpha0_LAB = rot_abCM_to_LAB * p4_alpha0_abCM; 
	TLorentzVector p4_alpha1_LAB = rot_abCM_to_LAB * p4_alpha1_abCM; 
	TLorentzVector p4_alpha2_LAB = rot_abCM_to_LAB * p4_alpha2_abCM; 

	TLorentzVector p4_total_a012_LAB  =  p4_alpha0_LAB + p4_alpha1_LAB+p4_alpha2_LAB;


	// Rotate Alpha-0 + 8(beryllium) from abCM to LAB
	TLorentzVector p4_a0_be8_LAB = rot_abCM_to_LAB * p4_total_abCM; 


	//========================================================================//


	std::cout<<" -- CHECK TOTAL MOMENTA -- "<<std::endl;
	std::cout<<"    g+12C    -> "<<p4_total_LAB.P() <<std::endl;
	std::cout<<"    a0+a1+a2 -> "<<p4_total_a012_LAB.P() <<std::endl;
	std::cout<<"    a0+8be   -> "<<p4_a0_be8_LAB.P() <<std::endl;
	std::cout<<"*******************************************************"<<std::endl;

	std::cout<<" -- CHECK TOTAL ENERGY -- "<<std::endl;
	std::cout<<"    g+12C    -> "<<p4_total_LAB.E()<<std::endl;
	std::cout<<"    a0+a1+a2 -> "<<p4_total_a012_LAB.E()  <<std::endl;
	std::cout<<"    a0+8be   -> "<<p4_a0_be8_LAB.E() <<"   <> "<< Etot_alpha0+Etot_8be<<std::endl;
	std::cout<<"*******************************************************"<<std::endl;

	//========================================================================//


	double energyPostBoost0 =  p4_alpha0_LAB.E()-p4_alpha0_LAB.M();
	double energyPostBoost1 =  p4_alpha1_LAB.E()-p4_alpha1_LAB.M();
	double energyPostBoost2 =  p4_alpha2_LAB.E()-p4_alpha2_LAB.M();

 // if(groun_or_exc==1)KE_8be=KE_8be+exct;


	std::cout<<"Total E -> "<<energyPostBoost0+energyPostBoost1+energyPostBoost2<<" (a0+a1+a2)"<<std::endl;
	if(groun_or_exc==0)std::cout<<"Total E -> "<<KE_alpha0+KE_8be<<" (a0+8be)"<<std::endl;
	if(groun_or_exc==1)std::cout<<"Total E -> "<<KE_alpha0+KE_8be+exct<<" (a0+8be)"<<std::endl;
	std::cout<<"Total E -> "<<energyPostBoost0+energyPostBoost1+energyPostBoost2-Q_VALUE_12C-Q_VALUE_8BE<<" (a0+a1+a2-Q1-Q2)"<<std::endl;
	if(groun_or_exc==0)std::cout<<"Total E -> "<<KE_alpha0+KE_8be-Q_VALUE_12C<<" (a0+8be-Q1)"<<std::endl;
	if(groun_or_exc==1)std::cout<<"Total E -> "<<KE_alpha0+KE_8be-Q_VALUE_12C+exct<<" (a0+8be-Q1)"<<std::endl;
	std::cout<<"*******************************************************"<<std::endl;
	std::cout<<"*******************************************************"<<std::endl;

	TLorentzVector totalLab = p4_alpha0_LAB+p4_alpha1_LAB+p4_alpha2_LAB;
	std::cout<<"Total P -> "<<totalLab.P()<<" (a0+a1+a2)"<<std::endl;
	std::cout<<"Total P -> "<<p4_total_LAB.P()<<" (g+12c)"<<std::endl;
	std::cout<<"*******************************************************"<<std::endl;
	std::cout<<"*******************************************************"<<std::endl;

	std::cout<<" -- CHECK INVAR MASS SQR -- "<<std::endl;
	std::cout<<"Total M2 -> "<<p4_total_LAB.Mag2()<<" (a0+a1+a2)"<<std::endl;
	std::cout<<"Total M2 -> "<<totalLab.Mag2()<<"  (g+12c)"<<std::endl;

	std::cout<<"*******************************************************"<<std::endl;
	std::cout<<"*******************************************************"<<std::endl;


  double rangePostBoost0 = myRangeCalculator.getIonRangeMM(pid_type::ALPHA,energyPostBoost0);
  double rangePostBoost1 = myRangeCalculator.getIonRangeMM(pid_type::ALPHA,energyPostBoost1);
  double rangePostBoost2 = myRangeCalculator.getIonRangeMM(pid_type::ALPHA,energyPostBoost2);

  TVector3 tangenta0 = p4_alpha0_LAB.Vect();
  TVector3 tangenta1 = p4_alpha1_LAB.Vect();
  TVector3 tangenta2 = p4_alpha2_LAB.Vect();

    tangenta0.SetMag(1.0);
      tangenta1.SetMag(1.0);
        tangenta2.SetMag(1.0);
  std::cout<<"EPB -> "<<energyPostBoost0<<" "<<energyPostBoost1<<" "<<energyPostBoost2<<std::endl;
  alphaData.setAlphaEnergy(0,energyPostBoost0);
  alphaData.setAlphaEnergy(1,energyPostBoost1);
  alphaData.setAlphaEnergy(2,energyPostBoost2);
  std::cout<<"RPB -> "<<rangePostBoost0<<" "<<rangePostBoost1<<" "<<rangePostBoost2<<std::endl;

  alphaData.setAlphaLength(0,rangePostBoost0);
  alphaData.setAlphaLength(1,rangePostBoost1);
  alphaData.setAlphaLength(2,rangePostBoost2);

  alphaData.setTangent(0,tangenta0);
  alphaData.setTangent(1,tangenta1);
  alphaData.setTangent(2,tangenta2);

  alphaData.setAlphaPhi(0,tangenta0.Phi());
  alphaData.setAlphaPhi(1,tangenta1.Phi());
  alphaData.setAlphaPhi(2,tangenta2.Phi());

  alphaData.setAlphaTheta(0,tangenta0.Theta());
  alphaData.setAlphaTheta(1,tangenta1.Theta());
  alphaData.setAlphaTheta(2,tangenta2.Theta());




  //tangent.SetMagThetaPhi(1.0, 0.803888, -0.34445);
  //length = 40.5754;



  tangent =  alphaData.getTangent(0);
  length =  alphaData.getAlphaLength(0);
  }
  if(index==1){

  //    tangent.SetMagThetaPhi(1.0, 1.76328, 3.11921);
  //length = 21.353;
    
  tangent =  alphaData.getTangent(1);
  length =  alphaData.getAlphaLength(1);
 

  } 
  if(index==2){
  //tangent.SetMagThetaPhi(1.0, 2.85048, 2.15383);
  //length = 17.4452;
  tangent =  alphaData.getTangent(2);
  length =  alphaData.getAlphaLength(2);
 

  } 


  TrackSegment3D aSegment;
  aSegment.setGeometry(myGeometryPtr);
  aSegment.setStartEnd(vertexPos, vertexPos + tangent*length);
  aSegment.setPID(ion_id);
  
  return aSegment;
}


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


Track3D EventSourceMC::createTrack(const TVector3 & aVtx, pid_type ion_id,AlphaData &alphaData, const int index, const int groun_or_exc ) const{

  Track3D aTrack;
  TrackSegment3D aSegment = create3pSegment(aVtx, ion_id,alphaData,index,groun_or_exc);
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
  double sigma =0.9;
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
  double sigma = 0.9;
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
  const int groun_or_exc = 1; // [0 - ground , 1- excited]

  TVector3 aVtx = createVertex();
  AlphaData alphaData;
 // int index = 0;

  myTracks3D.push_back(createTrack(aVtx, ion_id,alphaData,0,groun_or_exc));
  std::cout<<KBLU<<"Generated track number: "<<std::to_string(0)<<RST<<std::endl;
  std::cout<<myTracks3D.back()<<std::endl;
//  index++;

  myTracks3D.push_back(createTrack(aVtx, ion_id,alphaData,1,groun_or_exc));
  std::cout<<KBLU<<"Generated track number: "<<std::to_string(1)<<RST<<std::endl;
  std::cout<<myTracks3D.back()<<std::endl;
  //index++;

  myTracks3D.push_back(createTrack(aVtx, ion_id,alphaData,2,groun_or_exc));
  std::cout<<KBLU<<"Generated track number: "<<std::to_string(2)<<RST<<std::endl;
  std::cout<<myTracks3D.back()<<std::endl;

  std::cout<<"*****************"<<std::endl;
  std::cout<<"__ ALPHA 0 __ "<<std::endl;
  std::cout<<"length    = "<< alphaData.getAlphaLength(0)<<std::endl;
  std::cout<<"kinetic E = "<<alphaData.getAlphaEnergy(0) <<std::endl;
  std::cout<<"Theta     = "<<alphaData.getAlphaTheta(0) *180/TMath::Pi()<<std::endl;
  std::cout<<"Phi       = "<<alphaData.getAlphaPhi(0) <<std::endl;
  std::cout<<"__ ALPHA 1 __ "<<std::endl;
  std::cout<<"length    = "<<alphaData.getAlphaLength(1) <<std::endl;
  std::cout<<"kinetic E = "<<alphaData.getAlphaEnergy(1) <<std::endl;
  std::cout<<"Theta     = "<<alphaData.getAlphaTheta(1)*180/TMath::Pi() <<std::endl;
  std::cout<<"Phi       = "<<alphaData.getAlphaPhi(1) <<std::endl;
  std::cout<<"__ ALPHA 2 __ "<<std::endl;
  std::cout<<"length    = "<<alphaData.getAlphaLength(2) <<std::endl;
  std::cout<<"kinetic E = "<<alphaData.getAlphaEnergy(2) <<std::endl;
  std::cout<<"Theta     = "<<alphaData.getAlphaTheta(2)*180/TMath::Pi() <<std::endl;
  std::cout<<"Phi       = "<<alphaData.getAlphaPhi(2) <<std::endl;
  std::cout<<"*****************"<<std::endl;
  std::cout<<"Total E -> "<<alphaData.getAlphaEnergy(2)+alphaData.getAlphaEnergy(1)+alphaData.getAlphaEnergy(0)+7.367-92.0/1000 <<std::endl;

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
  
  //generateTwoProng();
  generateThreeProng();
  
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
double EventSourceMC::SRIMNormtrackLength(const double gammaEnergy, const double pressure, pid_type ion_id,AlphaData &alphaData, const int index)const{

  IonRangeCalculator myRangeCalculator(gas_mixture_type::CO2,pressure,293.15);
	
  double Kev_to_MeV = 1./1000.;
	// ** MASSES ** //
	const double AMU        = 931.49410242;
	//const double mass_c12   = AMU  *12.00000000000; 
	const double mass_8be   = AMU  * 8.00530510000; 
	const double mass_alpha = AMU  * 4.00260325413; 

  const double Q_VALUE_12C  = - 7367.0 * Kev_to_MeV; 
	const double Q_VALUE_8BE  = 92.0   * Kev_to_MeV; 

	// ** KE ** //
	const double ECM1 = gammaEnergy + Q_VALUE_12C; 
	const double KE_alpha0  = ECM1 * (mass_8be  /(mass_alpha+mass_8be));
	//const double KE_8be     = ECM1 * (mass_alpha/(mass_alpha+mass_8be));

	const double ECM2 = Q_VALUE_8BE; 
	const double KE_alpha12 = ECM2 * (mass_alpha/(mass_alpha+mass_alpha));







  double ke=0;

  if(ion_id==pid_type::ALPHA && index==0 ){
    ke = KE_alpha0;

  } else if(ion_id==pid_type::ALPHA && (index!=0)){
    ke = KE_alpha12;

  }


  double length = myRangeCalculator.getIonRangeMM(ion_id, ke);  
  alphaData.setAlphaEnergy(index, ke);
  alphaData.setAlphaLength(index, length);

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


double EventSourceMC::Boost2Lab3p(double const gammaEnergy,AlphaData &alphaData, int index)const{

  bool checkSwitch = true;

  const double AMU        = 931.49410242; // MeV/c^2
  const double alphaMass  = 4.00260325413 * AMU ;
  const double carbonMass = 12.000 * AMU;
  const double bery8Mass  =   8.00530510 * AMU ;
  double ke = gammaEnergy - 7.367; // Q value
  
  double alphaEnergyKe_0 = ke*(bery8Mass/(alphaMass+bery8Mass));
  double alphaEnergyKe_1_2 = (ke * (alphaMass/(alphaMass+bery8Mass)))/2;

  // ** BEAM/TARGE LAB ** //
  TVector3 beamDir_LAB(-1,0, 0); // unit vector
  double beamEnergy_LAB = gammaEnergy;
  TLorentzVector p4_beam_LAB(beamDir_LAB.Unit()*beamEnergy_LAB, beamEnergy_LAB); // [MeV/c, MeV/c, MeV/c, MeV/c^2]
  TLorentzVector p4_target_LAB(0, 0, 0, carbonMass); // [MeV/c, MeV/c, MeV/c, MeV/c^2]
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
  double p_alpha0_CMS      = sqrt(alphaEnergyKe_0*(alphaEnergyKe_0+2*alphaMass));
  double p_alpha1_CMS      = sqrt(alphaEnergyKe_1_2*(alphaEnergyKe_1_2+2*alphaMass));
  double p_alpha2_CMS      = sqrt(alphaEnergyKe_1_2*(alphaEnergyKe_1_2+2*alphaMass));
  double Etot_alpha0_CMS   = alphaMass+alphaEnergyKe_0;
  double Etot_alpha1_CMS   = alphaMass+alphaEnergyKe_1_2;
  double Etot_alpha2_CMS   = alphaMass+alphaEnergyKe_1_2;


  // ****************************************************** /
  // ****************************************************** /


  if(checkSwitch==true){
    std::cout<<KRED<<"********************************"<<RST<<std::endl;
    std::cout<<"Check that p4_total_LAB - > p4_total_CMS goes to zero under inverse boost"<<std::endl;  
    p4_total_LAB.Print();
    p4_total_CMS.Print();
    std::cout<<KRED<<"********************************"<<RST<<std::endl;
    std::cout<<"CARBON E in LAB: "<< p4_target_LAB.E()<<std::endl;
    std::cout<<"CARBON P in LAB: "<<p4_target_LAB.P()<<std::endl;
    std::cout<<"CARBON E in CMS: "<< p4_target_CMS.E()<<std::endl;
    std::cout<<"CARBON P in CMS: "<<p4_target_CMS.P()<<std::endl;
    std::cout<<KRED<<"********************************"<<RST<<std::endl;
    std::cout<<"GAMMA  E in LAB: "<<beamEnergy_LAB<<std::endl;
    std::cout<<"GAMMA  P in LAB: "<<p4_beam_LAB.P()<<std::endl;
    std::cout<<"GAMMA  E in CMS: "<< p4_beam_CMS.E()<<std::endl;
    std::cout<<"GAMMA  P in CMS: "<<p4_beam_CMS.P()<<std::endl;
    std::cout<<KRED<<"********************************"<<RST<<std::endl;
    std::cout<<"TOTAL  E in LAB: "<< p4_total_LAB.E()<<std::endl;
    std::cout<<"TOTAL  P in LAB: "<<p4_total_LAB.P()<<std::endl;
    std::cout<<"TOTAL  E in CMS: "<<p4_total_CMS.E()<<std::endl;
    std::cout<<"TOTAL  P in CMS: "<<p4_total_CMS.P()<<std::endl;
    std::cout<<KRED<<"********************************"<<RST<<std::endl;
    std::cout<<"Alpha0  KE in CMS: "<<alphaEnergyKe_0<<std::endl;
    std::cout<<"Alpha0  E in CMS: "<<Etot_alpha0_CMS<<std::endl;
    std::cout<<"Alpha0  P in CMS: "<<p_alpha0_CMS<<std::endl;
    std::cout<<KRED<<"********************************"<<RST<<std::endl;
    std::cout<<"Alpha1 KE in CMS: "<<alphaEnergyKe_1_2<<std::endl;
    std::cout<<"Alpha1 E in CMS: "<<Etot_alpha1_CMS<<std::endl;
    std::cout<<"Alpha1 P in CMS: "<<p_alpha1_CMS<<std::endl;
    std::cout<<KRED<<"********************************"<<RST<<std::endl;
    std::cout<<"Alpha2 KE in CMS: "<<alphaEnergyKe_1_2<<std::endl;
    std::cout<<"Alpha2 E in CMS: "<<Etot_alpha2_CMS<<std::endl;
    std::cout<<"Alpha2 P in CMS: "<<p_alpha2_CMS<<std::endl;
    std::cout<<KRED<<"********************************"<<RST<<std::endl;
  
  }
  TVector3 tangent = alphaData.getTangent(index);
  double energyPostBoost;
  double RangePostBoost;


  if(index==0){


    // ** BOOST ** //
    TLorentzVector alpha0_CMS(tangent * p_alpha0_CMS, Etot_alpha0_CMS);
    TLorentzVector alpha0_LAB = rot_CMS_to_LAB * alpha0_CMS;
    alpha0_LAB.Print();
    tangent = alpha0_LAB.Vect();
    energyPostBoost=alpha0_LAB.E()-alpha0_LAB.M();
    RangePostBoost  = myRangeCalculator.getIonRangeMM(pid_type::ALPHA,energyPostBoost);
    std::cout<<"E PRE -> "<<alphaData.getAlphaEnergy(index)<<std::endl;
    std::cout<<"E POST -> "<<energyPostBoost<<std::endl;
    std::cout<<"R PRE -> "<<alphaData.getAlphaLength(index)<<std::endl;
    std::cout<<"R POST -> "<<RangePostBoost<<std::endl;



  }else{
  // ** BOOST ** //

    // ** BEAM/TARGE LAB ** //
 // TVector3 b8_DIR = -alphaData.getTangent(0);
  //double b8_E_LAB =  2*alphaData.getAlphaEnergy(1) +bery8Mass;

  //TLorentzVector p4_beam_LAB2(b8_DIR.Unit()*b8_E_LAB, b8_E_LAB); // [MeV/c, MeV/c, MeV/c, MeV/c^2]
  //TLorentzVector p4_target_LAB2(0, 0, 0, carbonMass); // [MeV/c, MeV/c, MeV/c, MeV/c^2]
  //TLorentzVector p4_total_LAB2 =  p4_beam_LAB2 + p4_target_LAB2;
  //  TVector3 boost_CMS_to_LAB2 = p4_total_LAB2.BoostVector(); 
  //  TLorentzRotation rot_CMS_to_LAB2(boost_CMS_to_LAB2);

  TLorentzVector alpha1_CMS(tangent * p_alpha1_CMS, Etot_alpha1_CMS);
  TLorentzVector alpha1_LAB = rot_CMS_to_LAB * alpha1_CMS;
       //alpha1_LAB = rot_CMS_to_LAB2 * alpha1_LAB;

  alpha1_LAB.Print();
  tangent = alpha1_LAB.Vect();
  energyPostBoost=alpha1_LAB.E()-alpha1_LAB.M();
  RangePostBoost  = myRangeCalculator.getIonRangeMM(pid_type::ALPHA,energyPostBoost);
  std::cout<<"E PRE -> "<<alphaData.getAlphaEnergy(index)<<std::endl;
  std::cout<<"E POST -> "<<energyPostBoost<<std::endl;
  std::cout<<"R PRE -> "<<alphaData.getAlphaLength(index)<<std::endl;
  std::cout<<"R POST -> "<<RangePostBoost<<std::endl;





  }

  tangent.SetMag(1.0);
  alphaData.setTangent(index, tangent);
  alphaData.setAlphaEnergy(index, energyPostBoost);
  alphaData.setAlphaLength(index, RangePostBoost);
  return 0;
}




double EventSourceMC::BoostByVel3p(double const gammaEnergy,AlphaData &alphaData, int index)const{
  IonRangeCalculator myRangeCalculator(gas_mixture_type::CO2,250,293.15);
  const double AMU        = 931.49410242; // MeV/c^2
  const double alphaMass  = 4.00260325413 * AMU ;
  
  TVector3 b8_DIR = -alphaData.getTangent(0);
  double b8_E_LAB =  2*alphaData .getAlphaEnergy(1);// +bery8Mass;

  double be_p = sqrt(b8_E_LAB*(b8_E_LAB+2*(alphaMass+alphaMass)));
  double be_e = alphaMass+alphaMass+ b8_E_LAB;

  TLorentzVector p4_beam_LAB(b8_DIR.Unit()*(be_p), be_e); // [MeV/c, MeV/c, MeV/c, MeV/c^2]
  TLorentzVector p4_target_LAB1(0, 0, 0, alphaMass); // [MeV/c, MeV/c, MeV/c, MeV/c^2]
  TLorentzVector p4_target_LAB2(0, 0, 0, alphaMass); // [MeV/c, MeV/c, MeV/c, MeV/c^2]
  TLorentzVector p4_total_LAB =  p4_beam_LAB +p4_target_LAB1+ p4_target_LAB2;
  
  TVector3 boost_CMS_to_LAB2 = p4_total_LAB.BoostVector(); 
  TLorentzRotation rot_CMS_to_LAB(boost_CMS_to_LAB2);


  TVector3 boost_LAB_to_CMS = p4_total_LAB.BoostVector(); 
  TLorentzRotation rot_LAB_to_CMS(-boost_LAB_to_CMS);

  // ** BEAM/TARGET CMS ** //
  TLorentzVector p4_beam_CMS   =  rot_LAB_to_CMS * p4_beam_LAB;
  TLorentzVector p4_target1_CMS =  rot_LAB_to_CMS * p4_target_LAB1;
  TLorentzVector p4_target2_CMS =  rot_LAB_to_CMS * p4_target_LAB2;
  TLorentzVector p4_total_CMS  =  p4_beam_CMS + p4_target1_CMS+p4_target2_CMS;

 std::cout<<KRED<<"HERE_________________"<<std::endl;
  std::cout<<"Check that p4_total_LAB - > p4_total_CMS goes to zero under inverse boost"<<std::endl;  
  p4_total_LAB.Print();
  p4_total_CMS.Print();
  std::cout<<KRED<<"HERE_________________"<<std::endl;

  TVector3 tangent = alphaData.getTangent(index);
  double energyPostBoost;
  double RangePostBoost;
  TLorentzVector alpha1_CMS(tangent * (be_p/2), (be_e/2));
  TLorentzVector alpha1_LAB = rot_CMS_to_LAB * alpha1_CMS;
 alpha1_LAB.Print();
  tangent = alpha1_LAB.Vect();
  energyPostBoost=alpha1_LAB.E()-alpha1_LAB.M();
  RangePostBoost  = myRangeCalculator.getIonRangeMM(pid_type::ALPHA,energyPostBoost);
  std::cout<<"E PRE -> "<<alphaData.getAlphaEnergy(index)<<std::endl;
  std::cout<<"E POST -> "<<energyPostBoost<<std::endl;
  std::cout<<"R PRE -> "<<alphaData.getAlphaLength(index)<<std::endl;
  std::cout<<"R POST -> "<<RangePostBoost<<std::endl;

  tangent.SetMag(1.0);
  alphaData.setTangent(index, tangent);
  alphaData.setAlphaEnergy(index, energyPostBoost);
  alphaData.setAlphaLength(index, RangePostBoost);




return 0;  
}