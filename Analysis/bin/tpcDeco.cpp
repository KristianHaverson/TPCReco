#include <cstdlib>
#include <iostream>
#include <vector>

#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "TString.h"
#include "TStopwatch.h"

#include <boost/program_options.hpp>

#include "IonRangeCalculator.h"
#include "dEdxFitter.h"
#include "TrackBuilder.h"
#include "EventSourceMC.h"

#include "RecoOutput.h"
#include "RunIdParser.h"
#include "InputFileHelper.h"
#include "MakeUniqueName.h"
#include "colorText.h"
#include "TLorentzVector.h"
#include "TApplication.h"
#include "EventTPC.h"
#include "TSpectrum.h"

/*
 TH1D* SumHistograms(TH2D* h1, TH2D* h2, TH2D* h3) {

  TH1D* sum = new TH1D("","",290,10,280);
  int nbinsx = h1->GetNbinsX();

  int nbinsYu = h1->GetNbinsY();
  int nbinsYv = h2->GetNbinsY();
  int nbinsYw = h3->GetNbinsY();

  for (int i = 1; i <= nbinsx; ++i) {
	double sumVal=0;
	for (int y = 1; y <= nbinsYu; ++y){
		sumVal+= h1->GetBinContent(i, y);
	}
	for (int y = 1; y <= nbinsYv; ++y){
		sumVal+= h2->GetBinContent(i, y);
	}
	for (int y = 1; y <= nbinsYw; ++y){
		sumVal+= h3->GetBinContent(i, y);
	}


    sum->SetBinContent(i, sumVal);
  }


return sum;
}*/
/////////////////////////////////////
/////////////////////////////////////

/////////////////////////////////////
/////////////////////////////////////
int tpcDeco(const  std::string & geometryFileName,
		  const  std::string & dataFileName);
/////////////////////////////////////
/////////////////////////////////////
boost::program_options::variables_map parseCmdLineArgs(int argc, char **argv){

  boost::program_options::options_description cmdLineOptDesc("Allowed options");
  cmdLineOptDesc.add_options()
    ("help", "produce help message")
    ("geometryFile",  boost::program_options::value<std::string>()->required(), "string - path to the geometry file.")
    ("dataFile",  boost::program_options::value<std::string>()->required(), "string - path to data file.");
  
  boost::program_options::variables_map varMap;        
try {     
    boost::program_options::store(boost::program_options::parse_command_line(argc, argv, cmdLineOptDesc), varMap);
    if (varMap.count("help")) {
      std::cout << "makeTrackTree" << "\n\n";
      std::cout << cmdLineOptDesc << std::endl;
      exit(1);
    }
    boost::program_options::notify(varMap);
  } catch (const std::exception &e) {
    std::cerr << e.what() << '\n';
    std::cout << cmdLineOptDesc << std::endl;
    exit(1);
  }

  return varMap;
}
/////////////////////////////////////
/////////////////////////////////////
int main(int argc, char **argv){

  TStopwatch aStopwatch;
  aStopwatch.Start();

  std::string geometryFileName, dataFileName;
  boost::program_options::variables_map varMap = parseCmdLineArgs(argc, argv);
  boost::property_tree::ptree tree;
  if(argc<3){
    char text[] = "--help";
    char *argvTmp[] = {text, text};
    parseCmdLineArgs(2,argvTmp);
    return 1;
  }
  if (varMap.count("geometryFile")) {
    geometryFileName = varMap["geometryFile"].as<std::string>();
  }
  if (varMap.count("dataFile")) {
    dataFileName = varMap["dataFile"].as<std::string>();
  }

  int nEntriesProcessed = 0;
  if(dataFileName.size() && geometryFileName.size()){
    nEntriesProcessed = tpcDeco(geometryFileName, dataFileName);
  }
  else{
    std::cout<<KRED<<"Configuration not complete: "<<RST
	     <<" geometryFile: "<<geometryFileName<<"\n"
	     <<" dataFile: "<<dataFileName
	     <<std::endl;
  }

  aStopwatch.Stop();
  std::cout<<KBLU<<"Real time:       "<<RST<<aStopwatch.RealTime()<<" s"<<std::endl;
  std::cout<<KBLU<<"CPU time:        "<<RST<<aStopwatch.CpuTime()<<" s"<<std::endl;
  std::cout<<KBLU<<"Processing rate: "<<RST<<nEntriesProcessed/aStopwatch.RealTime()<< " ev/s"<<std::endl;

  return 0;
}

void readLineshape( std::vector<std::vector<float > >& data, std::string filename){


	std::ifstream myfile(filename, std::ios_base::in);
	if(myfile.fail()) std::cout<<"Lineshape img opening failed\n";
	if(myfile) std::cout<<"File opened \n";
	std::string a;

	std::ifstream in( filename );
	for ( std::string line; getline( in, line ); )
	{
		std::stringstream ss( line );
		std::vector<float> row;
		std::vector<std::vector<float>> indie;
		for ( float d; ss >> d; ) row.push_back( d );
		data.push_back( row );
	//	cout<<row<<endl;
	}



	myfile.close();
}

int  GetResHist(std::vector<std::vector<float> >& LineshapeTableImage,TH1D *RespHist ){

	
		for (unsigned int k=0; k < LineshapeTableImage[0].size(); k++) { // for each element
			RespHist->SetBinContent(k,LineshapeTableImage[0][k]);
	
	}



return 0;

}
#define MIN(X, Y) (((X) < (Y)) ? (X) : (Y))
#define MAX(X, Y) (((X) < (Y)) ? (Y) : (X))

double* convolve(double h[], double x[], int lenH, int lenX, int* lenY)
{
  int nconv = lenH+lenX-1;
  (*lenY) = nconv;
  int i,j,h_start,x_start,x_end;

double *y = (double*) calloc(nconv, sizeof(double));

  for (i=0; i<nconv; i++)
  {
    x_start = MAX(0,i-lenH+1);
    x_end   = MIN(i+1,lenX);
    h_start = MIN(i,lenH-1);
    for(j=x_start; j<x_end; j++)
    {
      y[i] += h[h_start--]*x[j];
    }
  }
  return y;
}

//using namespace RooFit;
/////////////////////////

/*
int tpcDeco(const  std::string & geometryFileName, const  std::string & dataFileName) {

  // ** EVENT SOURCE ** //
  std::shared_ptr<EventSourceMC> myEventSource = std::make_shared<EventSourceMC>(geometryFileName);
  
  // ** GEOMETRY ** //
  int index = geometryFileName.find("mbar");
  double pressure = stof(geometryFileName.substr(index-3, 3));
  IonRangeCalculator myRangeCalculator(gas_mixture_type::CO2,pressure,293.15);
  
	/////////////////////////////////////////////////////////////////////////////////////
	// ** READ RESPONSE /
	//std::vector<std::vector<float > > VRespo;
	//std::string VResponame ="/home/kris/Documents/ETPC/autoAnalysis/16_01_2023_/lineVRes.txt";
	//readLineshape(VRespo,VResponame);

	std::vector<std::vector<float > > URespo;
	std::string UResponame ="/home/kris/Documents/ETPC/autoAnalysis/16_01_2023_/lineURes.txt";
	readLineshape(URespo,UResponame);
  TH1D *URespHist = new TH1D("URespHist","URespHist",512,0,512);
  URespHist->GetXaxis()->SetRangeUser(0,20);
  GetResHist(URespo,URespHist);
  double ResponceU[512]; //URespHist
   for(int i=0;i<512;i++){
		  ResponceU[i] = URespHist->GetBinContent(i+1);
	  }

	std::vector<std::vector<float > > VRespo;
	std::string VResponame ="/home/kris/Documents/ETPC/autoAnalysis/16_01_2023_/lineVRes.txt";
	readLineshape(VRespo,VResponame);
  TH1D *VRespHist = new TH1D("VRespHist","VRespHist",512,0,512);
  VRespHist->GetXaxis()->SetRangeUser(0,20);
  GetResHist(VRespo,VRespHist);
  double ResponceV[512]; //URespHist
   for(int i=0;i<512;i++){
		  ResponceV[i] = VRespHist->GetBinContent(i+1);
	  }


	std::vector<std::vector<float > > WRespo;
	std::string WResponame ="/home/kris/Documents/ETPC/autoAnalysis/16_01_2023_/lineWRes.txt";
	readLineshape(WRespo,WResponame);
  TH1D *WRespHist = new TH1D("WRespHist","WRespHist",512,0,512);
  WRespHist->GetXaxis()->SetRangeUser(0,20);
  GetResHist(WRespo,WRespHist);
  double ResponceW[512]; //URespHist
   for(int i=0;i<512;i++){
		  ResponceW[i] = WRespHist->GetBinContent(i+1);
	  }


  // ** MAIN LOOP ** //
  unsigned int nEntries = myEventSource->numberOfEntries();
  nEntries = 400; //TEST


  TApplication application("",nullptr,nullptr);



  for(unsigned int iEntry=0;iEntry<nEntries;++iEntry){

    // ** COUNTER //
    std::cout<<KBLU<<"Processed: "<<iEntry<<"/"<<nEntries<<RST<<std::endl;
    
    // ** RECONSTRUCT ** //
    myEventSource->loadFileEntry(iEntry);    
    auto event = myEventSource->getCurrentEvent();
    //std::shared_ptr<TH2D> UHist2D = event->get2DProjection(projection_type::DIR_TIME_U, filter_type::none, scale_type::raw);
    //std::shared_ptr<TH2D> VHist2D = event->get2DProjection(projection_type::DIR_TIME_V, filter_type::none, scale_type::raw);
    //std::shared_ptr<TH2D> WHist2D = event->get2DProjection(projection_type::DIR_TIME_W, filter_type::none, scale_type::raw);

    std::shared_ptr<TH2D> UHist2D = event->get2DProjection(projection_type::DIR_TIME_U, filter_type::none, scale_type::raw);
	  TH2D *newUStrip = new TH2D("newUStrip","newUStrip",512,0,512, UHist2D->GetNbinsY(),0, UHist2D->GetNbinsY());	
	  TH2D *newUStripCON = new TH2D("newUStripCON","newUStripCON",512,0,512, UHist2D->GetNbinsY(),0, UHist2D->GetNbinsY());	

    std::shared_ptr<TH2D> VHist2D = event->get2DProjection(projection_type::DIR_TIME_V, filter_type::none, scale_type::raw);
	  TH2D *newVStrip = new TH2D("newVStrip","newVStrip",512,0,512, VHist2D->GetNbinsY(),0, VHist2D->GetNbinsY());	
	  TH2D *newVStripCON = new TH2D("newVStripCON","newVStripCON",512,0,512, VHist2D->GetNbinsY(),0, VHist2D->GetNbinsY());	

    std::shared_ptr<TH2D> WHist2D = event->get2DProjection(projection_type::DIR_TIME_W, filter_type::none, scale_type::raw);
	  TH2D *newWStrip = new TH2D("newWStrip","newWStrip",512,0,512, WHist2D->GetNbinsY(),0, WHist2D->GetNbinsY());	
	  TH2D *newWStripCON = new TH2D("newWStripCON","newWStripCON",512,0,512, WHist2D->GetNbinsY(),0, WHist2D->GetNbinsY());	



		for(int k=0; k<UHist2D->GetNbinsY();k++){ // For each strip
      TSpectrum *s = new TSpectrum();
      double Source[512]; 
			for(int i=0;i<512;i++) Source[i] = UHist2D->GetBinContent(i+1,k+1);			
      int lenY;
      double* conU =  convolve(ResponceU,Source, 512, 512, &lenY);
      for(int i=0; i<512;i++) newUStripCON->SetBinContent(i+1,k+1,conU[i]); 	 //fine
			s->Deconvolution(conU,ResponceU,512,100,2,2.);
			for(int i=0; i<512;i++) newUStrip->SetBinContent(i+1,k+1,conU[i]); 	 //fine
			for(int i=0; i<512;i++)if(newUStrip->GetBinContent(i+1,k+1)<=0)newUStrip->SetBinContent(i+1,k+1,0.1); 	 //
			for(int i=0; i<512;i++)if(newUStripCON->GetBinContent(i+1,k+1)<=0)newUStripCON->SetBinContent(i+1,k+1,0.1); 	 //
			for(int i=0; i<512;i++)if(UHist2D->GetBinContent(i+1,k+1)<=0)UHist2D->SetBinContent(i+1,k+1,0.1); 	 //
      delete s;
    }

    for(int k=0; k<VHist2D->GetNbinsY();k++){ // For each strip
      TSpectrum *s = new TSpectrum();
      double Source[512]; 
			for(int i=0;i<512;i++) Source[i] = VHist2D->GetBinContent(i+1,k+1);			
      int lenY;
      double* conV =  convolve(ResponceV,Source, 512, 512, &lenY);
      for(int i=0; i<512;i++) newVStripCON->SetBinContent(i+1,k+1,conV[i]); 	 //fine
			s->Deconvolution(conV,ResponceV,512,100,2,2.);
			for(int i=0; i<512;i++) newVStrip->SetBinContent(i+1,k+1,conV[i]); 	 //fine
			for(int i=0; i<512;i++)if(newVStrip->GetBinContent(i+1,k+1)<=0)newVStrip->SetBinContent(i+1,k+1,0.1); 	 //
			for(int i=0; i<512;i++)if(newVStripCON->GetBinContent(i+1,k+1)<=0)newVStripCON->SetBinContent(i+1,k+1,0.1); 	 //
			for(int i=0; i<512;i++)if(VHist2D->GetBinContent(i+1,k+1)<=0)VHist2D->SetBinContent(i+1,k+1,0.1); 	 //
      delete s;
    }


    for(int k=0; k<WHist2D->GetNbinsY();k++){ // For each strip
      TSpectrum *s = new TSpectrum();
      double Source[512]; 
			for(int i=0;i<512;i++) Source[i] = WHist2D->GetBinContent(i+1,k+1);			
      int lenY;
      double* conW =  convolve(ResponceW,Source, 512, 512, &lenY);
      for(int i=0; i<512;i++) newWStripCON->SetBinContent(i+1,k+1,conW[i]); 	 //fine
			s->Deconvolution(conW,ResponceW,512,100,2,2.);
			for(int i=0; i<512;i++) newWStrip->SetBinContent(i+1,k+1,conW[i]); 	 //fine
			for(int i=0; i<512;i++)if(newWStrip->GetBinContent(i+1,k+1)<=0)newWStrip->SetBinContent(i+1,k+1,0.1); 	 //
			for(int i=0; i<512;i++)if(newWStripCON->GetBinContent(i+1,k+1)<=0)newWStripCON->SetBinContent(i+1,k+1,0.1); 	 //
			for(int i=0; i<512;i++)if(WHist2D->GetBinContent(i+1,k+1)<=0)WHist2D->SetBinContent(i+1,k+1,0.1); 	 //
      delete s;
    }
    


//   TH1D* sum = SumHistograms(newUStripCON,newVStripCON,newWStripCON);


//
   // newUStripCON->GetXaxis()->SetRangeUser(10,280);
   // newUStripCON->GetYaxis()->SetRangeUser(28,79);
//
   // newVStripCON->GetXaxis()->SetRangeUser(10,280);
   // newVStripCON->GetYaxis()->SetRangeUser(34,105);
//
   // newWStripCON->GetXaxis()->SetRangeUser(10,280);
   // newWStripCON->GetYaxis()->SetRangeUser(63,105);



    TCanvas *canvas = new TCanvas("canavs","canvas",1400,600);
		canvas->SetWindowPosition(-100,-100);
		canvas->Divide(3,1);
		canvas->cd(1);
		newUStripCON->Draw("colz");
		canvas->cd(2);
		newVStripCON->Draw("colz");
		canvas->cd(3);
		newWStripCON->Draw("colz");
    newUStripCON->SetStats(0);
    //canvas->cd(4);
    //sum->Draw("HIST");

    newVStripCON->SetStats(0);
    newWStripCON->SetStats(0);

		canvas->Modified();
		canvas->Update();

    char step;
		std::cout<<"Step forward? y/n"<<std::endl;
		std::cin>>step;
		if(step=='n')return 0 ;

    //delete newUStrip;
    //delete newUStripCON;
    //delete newWStrip;
    //delete newWStripCON;
    //delete newVStrip;
    //delete newVStripCON;

  }
  application.Run();
  return nEntries;
}
/////////////////////////////
////////////////////////////

*/



typedef struct {Float_t eventId, frameId,
    eventTypeGen,
    alphaRangeGen,
    alphaEnergyGen,
    chargeGen,
    cosThetaGen, phiGen,
    ///
    genVertexX,genVertexY,genVertexZ,
    } TrackData;


int tpcDeco(const  std::string & geometryFileName, const  std::string & dataFileName) {




}