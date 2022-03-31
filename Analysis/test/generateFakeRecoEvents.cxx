#ifndef __ROOTLOGON__
R__ADD_INCLUDE_PATH(../../DataFormats/include)
R__ADD_INCLUDE_PATH(../../Reconstruction/include)
R__ADD_LIBRARY_PATH(../lib)
#endif
#include "TROOT.h"
#include "Track3D.h"
#include "MultiKey.h"
#include "GeometryTPC.h"
#include "IonRangeCalculator.h"
//////////////////////////
//
// root
// root [0] .L generateFakeRecoEvents.cxx
// root [1] generateFakeRecoEvents("geometry.dat")
//
//
//
//////////////////////////
std::shared_ptr<IonRangeCalculator> loadRangeCalculator(){

  std::shared_ptr<IonRangeCalculator> db = std::make_shared<IonRangeCalculator>();

  // set current conditions: gas=CO2, pressure=250 mbar, temperature=20C
  db->setGasConditions(IonRangeCalculator::CO2, 250.0, 273.15+20);
  
  return db;
}
/////////////////////////
/////////////////////////
std::shared_ptr<GeometryTPC> loadGeometry(const std::string fileName){
  return std::make_shared<GeometryTPC>(fileName.c_str(), false);
}
//////////////////////////
//////////////////////////
bool isFullyContainedEvent(Track3D & aTrack){
  if(aTrack.getSegments().size()==0) return false;
  std::shared_ptr<GeometryTPC> aGeometry=aTrack.getSegments().front().getGeometry();
  double xmin, xmax, ymin, ymax, zmin, zmax;
  std::tie(xmin, xmax, ymin, ymax)=aGeometry->rangeXY();
  zmin=aGeometry->GetDriftCageZmin();
  zmax=aGeometry->GetDriftCageZmax();
  
  for(auto & aSegment: aTrack.getSegments()) {
    auto pack = aSegment.getStartEndXYZ(); // size = 6
    if(pack[0]<xmin || pack[0]>xmax || pack[0+3]<xmin || pack[0+3]>xmax || // x-coordinates
       pack[1]<ymin || pack[1]>ymax || pack[1+3]<ymin || pack[1+3]>ymax || // y-coordinates
       pack[2]<zmin || pack[2]>zmax || pack[2+3]<zmin || pack[2+3]>zmax)   // z-coordinates
      return false;
  }
  return true;
}
//////////////////////////
//////////////////////////
Track3D generateFake3AlphaEvent(std::shared_ptr<GeometryTPC> aGeometry, std::shared_ptr<IonRangeCalculator> aRangeCalculator, bool debug_flag=false){
  //Track3D generateFake3AlphaEvent(std::shared_ptr<GeometryTPC> aGeometry, IonRangeCalculator *aRangeCalculator, bool debug_flag=false){

  //  auto r=new Trandom3();
  auto r=gRandom; // use global ROOT pointer

  Track3D empty_result;

  ///Data and geometry loading.
  //  std::string fileName = "Reco_EventTPC_2021-06-22T12:01:56.568_0.root";
  //  Track3D *aTrack = loadRecoEvent(fileName);
  //  fileName = "/home/user1/scratch/akalinow/ELITPC/TPCReco/build/resources/geometry_ELITPC_250mbar_12.5MHz.dat";
  //  std::shared_ptr<GeometryTPC> aGeometry = loadGeometry(fileName);
  //  std::cout<<*aTrack<<std::endl;

  ///Fetching tracks segments from the full track.
  //  TrackSegment3D aSegment = aTrack->getSegments().front();
  //  aSegment.setGeometry(aGeometry);

  // define some constants
  double atomicMassUnit = 931.49410242; //MeV/c^2
  double alphaMass = 4*atomicMassUnit + 2.4249; //A*u + Delta MeV/c^2
  double carbonMass = 12*atomicMassUnit;

  // detector reference frame (DET=LAB)
  double beamEnergyResolution=0.05; // LAB beam energy smearing factor
  double beamEnergy_DET=13.0*r->Gaus(1, beamEnergyResolution); // smear beam energy by 5% // MeV
  TVector3 beamDir_DET(1, 0, 0); // unit vector
  TLorentzVector photonP4_DET(beamDir_DET.Unit()*beamEnergy_DET, beamEnergy_DET); // [MeV/c, MeV/c, MeV/c, MeV/c^2]
  TLorentzVector carbonP4_DET(0, 0, 0, carbonMass); // [MeV/c, MeV/c, MeV/c, MeV/c^2]

  // speed of photon-carbon CMS frame in DET frame
  double beta=beamEnergy_DET/(beamEnergy_DET+carbonMass);
  TVector3 beta_DET=beamDir_DET.Unit()*beta;

  // boosting P4 from DET/LAB frame to CMS frame (see TlorentzVector::Boost() convention!)
  TLorentzVector photonP4_CMS(photonP4_DET); 
  TLorentzVector carbonP4_CMS(carbonP4_DET); 
  photonP4_CMS.Boost(-beta_DET);
  carbonP4_CMS.Boost(-beta_DET);

  // check total energy in CMS frame
  double totalEnergy_CMS=(photonP4_CMS+carbonP4_CMS).E();
  double totalEnergy_CMS_xcheck=(carbonMass+beamEnergy_DET*(1-beta))/sqrt(1-beta*beta); // DEBUG
  double totalEnergy_CMS_xcheck2=sqrt( carbonMass*(2*beamEnergy_DET+carbonMass) ); // DEBUG

  // check beam energy in CMS frame
  double beamEnergy_CMS=photonP4_CMS.E();
  double beamEnergy_CMS_xcheck=beamEnergy_DET*sqrt(carbonMass/(carbonMass+2*beamEnergy_DET));

  // stationary excited carbon state
  double carbonMassExcited=totalEnergy_CMS;
  double carbonMassExcited_xcheck=beamEnergy_CMS_xcheck+sqrt(carbonMass*carbonMass+beamEnergy_CMS*beamEnergy_CMS); //
  double Qvalue_CMS=carbonMassExcited-3*alphaMass;
  if(Qvalue_CMS<0) {
    if(debug_flag) std::cout<<"Beam energy is too low to create 3 alpha particles!"<<std::endl;
    return empty_result;
  }
    
  // assign randomly kinetic energies for 3 alpha particles:
  //
  // - matrix element is constant
  // - Q = T1+T2+T3 = height of a regular triangle (Dalitz plot for same-particle 3-body decay)
  // - inscribed circle of radius R=Q/3 is uniformely populated (non-relativistic case, Ti<<M)
  // - T3=R+r0*sin(f0), 0<=r0<=R, 0<=f0<2*pi
  // - T2=R+r0*sin(f0+2*pi/3)
  // - T1=R+r0*sin(f0+4*pi/3)=Q-T2-T3
  //
  double phi0=r->Uniform(0, TMath::TwoPi());
  double r0=sqrt(r->Uniform(0, pow(Qvalue_CMS/3, 2))); // dS=dPhi*d(r^2)/2
  double T1_CMS=Qvalue_CMS/3+r0*sin(phi0);
  double T2_CMS=Qvalue_CMS/3+r0*sin(phi0+TMath::TwoPi()/3);
  double T3_CMS=Qvalue_CMS-T1_CMS-T2_CMS;
  double T3_CMS_xcheck=Qvalue_CMS/3+r0*sin(phi0+2*TMath::TwoPi()/3);
  double p1_CMS=sqrt(T1_CMS*(T1_CMS+2*alphaMass)); // sqrt(2*T1_CMS*alphaMass); // MeV
  double p2_CMS=sqrt(T2_CMS*(T2_CMS+2*alphaMass)); // sqrt(2*T2_CMS*alphaMass); // MeV
  double p3_CMS=sqrt(T3_CMS*(T3_CMS+2*alphaMass)); // sqrt(2*T3_CMS*alphaMass); // MeV

  // distribute isotropically 1st alpha particle momentum:
  double phi1=r->Uniform(0, TMath::TwoPi()); // rad
  double theta1=acos(r->Uniform(-1, 1)); // rad, dOmega=dPhi*d(cosTheta)
  TLorentzVector alpha1_P4_CMS(p1_CMS*TVector3(0, 0, 1), alphaMass+T1_CMS); // initial = along Z
  alpha1_P4_CMS.RotateY(theta1); // rotate by theta1 about Y axis
  alpha1_P4_CMS.RotateZ(phi1); // rotate by phi1 about Z axis

  // distribute isotropically 2nd alpha particle momentum
  // - cos(theta12) is constrained by momentum conservation: p3^2 = p1^2 + p2^2 + 2*p1*p2*cos(theta12)
  double phi12=r->Uniform(0, TMath::TwoPi()); // rad
  double theta12=acos( 0.5*(p3_CMS*p3_CMS - p1_CMS*p1_CMS - p2_CMS*p2_CMS)/p1_CMS/p2_CMS ); // rad
  TLorentzVector alpha2_P4_CMS(p2_CMS*alpha1_P4_CMS.Vect().Unit(), alphaMass+T2_CMS); // initial = along ALPHA1
  alpha2_P4_CMS.Rotate(theta12, alpha1_P4_CMS.Vect().Orthogonal()); // rotate by theta12 about axis orthogonal to ALPHA1
  alpha2_P4_CMS.Rotate(phi12, alpha1_P4_CMS.Vect()); // rotate by phi12 about ALPHA1 axis

  // calculate 3rd alpha particle momentum: P3=-(P1+P2)
  //  TLorentzVector alpha3_P4_CMS=-alpha1_P4_CMS-alpha2_P4_CMS; // sets correct momentum but not energy
  //  alpha3_P4_CMS.SetE(alphaMass+T3_CMS); // sets correct energy
  TLorentzVector alpha3_P4_CMS;
  alpha3_P4_CMS.SetVectM( -(alpha1_P4_CMS.Vect()+alpha2_P4_CMS.Vect()), alphaMass);

  // boosting P4 from CMS frame to DET/LAB frame (see TlorentzVector::Boost() convention!)
  TLorentzVector alpha1_P4_DET(alpha1_P4_CMS); 
  TLorentzVector alpha2_P4_DET(alpha2_P4_CMS); 
  TLorentzVector alpha3_P4_DET(alpha3_P4_CMS); 
  alpha1_P4_DET.Boost(beta_DET);
  alpha2_P4_DET.Boost(beta_DET);
  alpha3_P4_DET.Boost(beta_DET);

  // print debug info
  //  bool debug_flag=true;
  if(debug_flag) {
    TLorentzVector alpha3_P4_CMS_xcheck( -(alpha1_P4_CMS+alpha2_P4_CMS).Vect(), alphaMass+T3_CMS);
    TLorentzVector alpha3_P4_CMS_xcheck2( -(alpha1_P4_CMS.Vect()+alpha2_P4_CMS.Vect()).Unit()*p3_CMS, alphaMass+T3_CMS);
    std::cout<<"Alpha[1] kin.energy in CMS (2 methods): "<<T1_CMS<<", "<<alpha1_P4_CMS.E()-alpha1_P4_CMS.M()<<std::endl;
    std::cout<<"Alpha[1] momentum in CMS (2 methods): "<<p1_CMS<<", "<<alpha1_P4_CMS.P()<<std::endl;
    std::cout<<"Alpha[2] kin.energy in CMS (2 methods): "<<T2_CMS<<", "<<alpha2_P4_CMS.E()-alpha2_P4_CMS.M()<<std::endl;
    std::cout<<"Alpha[2] momentum in CMS (2 methods): "<<p2_CMS<<", "<<alpha2_P4_CMS.P()<<std::endl;
    std::cout<<"Alpha[3] kin.energy in CMS (3 methods): "<<T3_CMS<<", "<<alpha3_P4_CMS.E()-alpha3_P4_CMS.M()<<", "<<T3_CMS_xcheck<<std::endl;
    std::cout<<"Alpha[3] momentum in CMS (3 methods): "<<p3_CMS<<", "<<alpha3_P4_CMS.P()<<", "<<alpha3_P4_CMS_xcheck.P()<<std::endl;
    std::cout<<"Alpha[3] P4 in CMS (3 methods): "<<std::endl;
    alpha3_P4_CMS_xcheck2.Print();
    alpha3_P4_CMS.Print();
    alpha3_P4_CMS_xcheck.Print();
    std::cout<<"Alpha T1+T2+T3 in CMS (2 methods): "<<T1_CMS+T2_CMS+T3_CMS<<", "<<alpha1_P4_CMS.E()+alpha2_P4_CMS.E()+alpha3_P4_CMS.E()-alpha1_P4_CMS.M()-alpha2_P4_CMS.M()-alpha3_P4_CMS.M()<<std::endl;
    std::cout<<"Alpha p1+p2+p3 momentum in CMS (2 methods):"<<std::endl;
    (alpha1_P4_CMS.Vect()+alpha2_P4_CMS.Vect()+alpha3_P4_CMS_xcheck2.Vect()).Print();
    (alpha1_P4_CMS.Vect()+alpha2_P4_CMS.Vect()+alpha3_P4_CMS.Vect()).Print();
    (alpha1_P4_CMS.Vect()+alpha2_P4_CMS.Vect()+alpha3_P4_CMS_xcheck.Vect()).Print();
    std::cout<<"Alpha |p1+p2| momentum in CMS (2 methods):"<<p3_CMS<<", "<<(alpha1_P4_CMS+alpha2_P4_CMS).P()<<std::endl;
    double invariantMass1 = (photonP4_CMS+carbonP4_CMS).M();
    double invariantMass2 = (alpha1_P4_CMS+alpha2_P4_CMS+alpha3_P4_CMS).M();
    std::cout<<"Gamma beam energy in DET/LAB: "<<beamEnergy_DET<<std::endl;
    std::cout<<"Gamma beam energy in CMS (2 methods): "<<beamEnergy_CMS<<", "<<beamEnergy_CMS_xcheck<<std::endl;
    std::cout<<"Total energy in CMS (3 methods): "<<totalEnergy_CMS<<" "<<totalEnergy_CMS_xcheck<<" "<<totalEnergy_CMS_xcheck2<<std::endl; // DEBUG
    std::cout<<"Q-value in CMS: "<<Qvalue_CMS<<std::endl;
    std::cout<<"Alpha[i] kin.energy in CMS: "<<T1_CMS<<", "<<T2_CMS<<", "<<T3_CMS<<std::endl;
    std::cout<<"Alpha[i] mass in CMS: "<<alpha1_P4_CMS.M()<<" "<<alpha2_P4_CMS.M()<<" "<<alpha3_P4_CMS.M()<<std::endl; // DEBUG
    std::cout<<"Photon+Carbon system invariant mass = "<<invariantMass1<<std::endl;
    std::cout<<"Triple-alpha system invariant mass = "<<invariantMass2<<std::endl;
    std::cout<<"Carbon ground state mass = "<<carbonMass<<std::endl;
    std::cout<<"Carbon excited state mass (2 methods)= "<<carbonMassExcited<<", "<<carbonMassExcited_xcheck<<std::endl;
    std::cout<<"Alpha kin.energy in DET/LAB: "<<alpha1_P4_DET.E()-alpha1_P4_DET.M()<<", "<<alpha2_P4_DET.E()-alpha2_P4_DET.M()<<", "<<alpha3_P4_DET.E()-alpha3_P4_DET.M()<<std::endl;
    std::cout<<"Alpha mass in DET/LAB: "<<alpha1_P4_DET.M()<<" "<<alpha2_P4_DET.M()<<" "<<alpha3_P4_DET.M()<<std::endl; // DEBUG
  }

  // calculate alpha ranges [mm] in LAB/DET frame
  double range1_DET=aRangeCalculator->getIonRangeMM(IonRangeCalculator::ALPHA, alpha1_P4_DET.E()-alpha1_P4_DET.M()); // mm
  double range2_DET=aRangeCalculator->getIonRangeMM(IonRangeCalculator::ALPHA, alpha2_P4_DET.E()-alpha2_P4_DET.M()); // mm
  double range3_DET=aRangeCalculator->getIonRangeMM(IonRangeCalculator::ALPHA, alpha3_P4_DET.E()-alpha3_P4_DET.M()); // mm

  // create list of tracks
  std::vector<TVector3> list;
  list.push_back(alpha1_P4_DET.Vect().Unit()*range1_DET); // 1st track direction
  list.push_back(alpha2_P4_DET.Vect().Unit()*range2_DET); // 2nd track direction
  list.push_back(alpha3_P4_DET.Vect().Unit()*range3_DET); // 3rd track direction
  // sort tracks in ascending Z-coordinate order
  std::sort(list.begin(), list.end(),
	    [](const TVector3& a, const TVector3& b) {
      return a.Z() < b.Z();
    });

  // smear vertex position in XYZ
  // - assume gaussian spread of beam profile
  // - assume self-triggering mode
  // - assume trigger position at 10% of TIME scale
  double xmin, xmax, ymin, ymax, zmin, zmax;
  std::tie(xmin, xmax, ymin, ymax)=aGeometry->rangeXY();
  zmin=aGeometry->GetDriftCageZmin();
  zmax=aGeometry->GetDriftCageZmax();
  double beamSpread=5.0; // sigma in [mm]
  TVector3 vertex( r->Uniform(xmin, xmax), r->Gaus(0, beamSpread), r->Gaus(0.5*(zmin+zmax), beamSpread) );
  vertex.SetZ(zmin+0.1*(zmax-zmin)-list.front().Z()); // simulates 1st (triggering) track arrival at 10% of TIME scale
  
  // create TrackSegment3D collection
  Track3D aTrack;
  TrackSegment3D aSegment;
  for(auto & seg: list) {
    aSegment.setGeometry(aGeometry);
    aSegment.setStartEnd(vertex, vertex+seg);
    aTrack.addSegment(aSegment);
  }

  // check if all segments are contained
  if(!isFullyContainedEvent(aTrack)) {
    if(debug_flag) std::cout<<"Event is not fully contained inside the active volume!"<<std::endl;
    return empty_result;
  }	

  // accept event
  return aTrack;
}
//////////////////////////
//////////////////////////
void generateFakeRecoEvents(const std::string geometryName, bool debug_flag=false){

  if (!gROOT->GetClass("Track3D")){
    R__LOAD_LIBRARY(libDataFormats.so);
  }
  if (!gROOT->GetClass("IonRangeCalculator")){
    R__LOAD_LIBRARY(libReconstruction.so);
  }

  const std::string fileName="Reco_FakeEvents.root";
  auto myGeometry=loadGeometry(geometryName);
  if(!myGeometry->IsOK()) {
    std::cerr<<"ERROR: Wrong geometry file!"<<std::endl;
    return;
  }
  auto myRangeCalculator=loadRangeCalculator();
  if(!myRangeCalculator->IsOK()) {
    std::cerr<<"ERROR: Cannot initialize range/energy calculator!"<<std::endl;
    return;
  }
  
  TFile *aFile = new TFile(fileName.c_str(), "CREATE");
  if(!aFile || !aFile->IsOpen()) {
    std::cerr<<"ERROR: Cannot create output file!"<<std::endl;
    return;
  }
  aFile->cd();
  TTree *aTree = new TTree("TPCRecoData","");
  Track3D *aTrack = new Track3D();
  aTree->Branch("RecoEvent", "Track3D", &aTrack);

  for(auto i=0L; i<10000L; i++) {
    if(debug_flag || (i<1000L && (i%100L)==0L) || (i%1000L)==0L ) std::cout<<"Generating fake track i="<<i<<std::endl;
    *aTrack = generateFake3AlphaEvent(myGeometry, myRangeCalculator, debug_flag);

    // skip bad or parially contained events
    if(aTrack->getSegments().size()==0) {
      if(debug_flag) std::cout<<"Fake event i="<<i<<" rejected!"<<std::endl;
      continue; 
    }
    aTree->Fill();
  }
  aTree->Write("",TObject::kOverwrite);
  aFile->Close();
}
