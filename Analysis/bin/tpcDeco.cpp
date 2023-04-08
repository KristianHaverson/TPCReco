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
#include "TStyle.h"
#include "TSystem.h"




 TH1D* SumHistograms( std::shared_ptr<TH2D> h1,  std::shared_ptr<TH2D> h2,  std::shared_ptr<TH2D>h3) {

  TH1D* sum = new TH1D("","",512,0,512);
  
  TH1D* temp1 = h1->ProjectionX();// new TH1D("temp1","temp1",512,0,512);
  TH1D* temp2 = h2->ProjectionX();// new TH1D("temp2","temp2",512,0,512);
  TH1D* temp3 = h3->ProjectionX();// new TH1D("temp3","temp3",512,0,512);

  for (int i = 1; i <= 512; ++i) {
      
    sum->SetBinContent(i,temp1->GetBinContent(i)+temp2->GetBinContent(i)+temp3->GetBinContent(i) );
  }

  /*
  int nbinsx = h1->GetNbinsX();

  int nbinsYu = h1->GetNbinsY();
  int nbinsYv = h2->GetNbinsY();
  int nbinsYw = h3->GetNbinsY();*/


/*
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
  }*/

  delete temp1;
  delete temp2;
  delete temp3;
return sum;
}
/////////////////////////////////////
/////////////////////////////////////

/////////////////////////////////////
/////////////////////////////////////
int tpcDeco(const  std::string & geometryFileName,
		  const  std::string & dataFileName);

int reFit(const  std::string & geometryFileName,
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
   // nEntriesProcessed = tpcDeco(geometryFileName, dataFileName);
    nEntriesProcessed = reFit(geometryFileName, dataFileName);
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









void convolve(TH1D *hist1 ,TH1D *hist2_gaus,TH1D *conv_hist  ){





  // Loop over each bin in the new histogram
  for (int i = 1; i <= 512; i++) {
      double sum = 0;
      
      // Calculate the convolution for this bin
      for (int j = 1; j <= 512; j++) {
          if (i-j >= 1 && i-j <= 512) {
              sum += hist1->GetBinContent(j) * hist2_gaus->GetBinContent(i-j);
          }
      }
      
      // Set the bin content of the new histogram to the result of the convolution
      conv_hist->SetBinContent(i, sum);
  }

  //delete hist2_gaus;
  return;

}
 

//using namespace RooFit;
/////////////////////////


int tpcDeco(const  std::string & geometryFileName, const  std::string & dataFileName) {

  // ** EVENT SOURCE ** //
  
  //std::shared_ptr<EventSourceGRAW> myEventSource = std::make_shared<EventSourceGRAW>(geometryFileName);
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

/*
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

*/

  // ** MAIN LOOP ** //
  unsigned int nEntries = myEventSource->numberOfEntries();
  nEntries = 10000; //TEST


  TApplication application("",nullptr,nullptr);
  	gStyle->SetOptStat(0.);

/*
	gStyle->SetOptStat(0.);
	gStyle->SetPadTopMargin(0.);
	gStyle->SetPadRightMargin(0.);
	gStyle->SetPadBottomMargin(0.);
	gStyle->SetPadLeftMargin(0.);
  gStyle->SetPadTickY(0);
  gStyle->SetPadTickX(0);
  gStyle->SetPadBorderMode(0);
  gStyle->SetPadColor(0);
	gStyle->SetPalette(  kGreyScale);*/
//	gStyle->SetPalette(  kGreyScale);

  
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
//	  TH2D *newUStrip = new TH2D("newUStrip","newUStrip",512,0,512, UHist2D->GetNbinsY(),0, UHist2D->GetNbinsY());	
//	  TH2D *newUStripCON = new TH2D("newUStripCON","newUStripCON",512,0,512, UHist2D->GetNbinsY(),0, UHist2D->GetNbinsY());	

    std::shared_ptr<TH2D> VHist2D = event->get2DProjection(projection_type::DIR_TIME_V, filter_type::none, scale_type::raw);
	//  TH2D *newVStrip = new TH2D("newVStrip","newVStrip",512,0,512, VHist2D->GetNbinsY(),0, VHist2D->GetNbinsY());	
	//  TH2D *newVStripCON = new TH2D("newVStripCON","newVStripCON",512,0,512, VHist2D->GetNbinsY(),0, VHist2D->GetNbinsY());	

    std::shared_ptr<TH2D> WHist2D = event->get2DProjection(projection_type::DIR_TIME_W, filter_type::none, scale_type::raw);
	  //TH2D *newWStrip = new TH2D("newWStrip","newWStrip",512,0,512, WHist2D->GetNbinsY(),0, WHist2D->GetNbinsY());	
	  //TH2D *newWStripCON = new TH2D("newWStripCON","newWStripCON",512,0,512, WHist2D->GetNbinsY(),0, WHist2D->GetNbinsY());	


/*
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
    */





    int nbinsx = UHist2D->GetNbinsX();
    
    int nbinsYu = UHist2D->GetNbinsY();
    int nbinsYv = VHist2D->GetNbinsY();
    int nbinsYw = WHist2D->GetNbinsY();


    //const double timeRes_U= gRandom->Gaus(0,2.536); //[timebins 13.1 MeV]
    //const double timeRes_V= gRandom->Gaus(0,2.9425); //[timebins 13.1 MeV]
    //const double timeRes_W= gRandom->Gaus(0,2.793); //[timebins 13.1 MeV]

    // /const double imageRes_U = gRandom->Gaus(0,0.4951); // [strips]
    // /const double imageRes_V = gRandom->Gaus(0,0.4959); // [strips]
    // /const double imageRes_W = gRandom->Gaus(0,0.4528); // [strips]

  // TH2D* copyU1 = (TH2D*)UHist2D->Clone();

    for (int i = 1; i <= nbinsx; ++i) {

      for (int y = 1; y <= nbinsYu; ++y){
       // if(UHist2D->GetBinContent(i,y)<0.01)UHist2D->SetBinContent(i,y,0.001);
        UHist2D->SetBinContent(i,y,UHist2D->GetBinContent(i,y)+gRandom->Gaus(0,1.29821e+02/2)   );
       // UHist2D->SetBinContent(i,y,UHist2D->GetBinContent(i,y)+gRandom->Gaus(2.68196e+02,1.29821e+02)   );
      }

      for (int y = 1; y <= nbinsYv; ++y){
        //if(VHist2D->GetBinContent(i,y)<0.01)VHist2D->SetBinContent(i,y,0.001)
        VHist2D->SetBinContent(i,y,VHist2D->GetBinContent(i,y)+gRandom->Gaus(0,1.31257e+02/2)   );
       // VHist2D->SetBinContent(i,y,VHist2D->GetBinContent(i,y)+gRandom->Gaus(2.68422e+02,1.31257e+02/6)   );

      }


      for (int y = 1; y <= nbinsYw; ++y){
       // if(WHist2D->GetBinContent(i,y)<0.01)WHist2D->SetBinContent(i,y,0.001);
        WHist2D->SetBinContent(i,y,WHist2D->GetBinContent(i,y)+gRandom->Gaus( 0,1.30251e+02/2)   );
       // WHist2D->SetBinContent(i,y,WHist2D->GetBinContent(i,y)+gRandom->Gaus( 2.68416e+02,1.30251e+02/6)   );

      }

    }
   // TH2D* copyU2 = (TH2D*)UHist2D->Clone();





      TH1D *hist2_gausU = new TH1D("hist2_gausU", "hist2_gausU", 512, 0, 512);
      TH1D *hist2_gausV = new TH1D("hist2_gausV", "hist2_gausV", 512, 0, 512);
      TH1D *hist2_gausW = new TH1D("hist2_gausW", "hist2_gausW", 512, 0, 512);


      for (int j = 1; j <= 10000; j++)hist2_gausU->Fill(gRandom->Gaus(0,3.98121));
      for (int j = 1; j <= 10000; j++)hist2_gausV->Fill(gRandom->Gaus(0,4.9392));
      for (int j = 1; j <= 10000; j++)hist2_gausW->Fill(gRandom->Gaus(0,4.9282));
              




      for (int y = 1; y <= nbinsYu; ++y){
        TH1D *conv_hist = new TH1D("conv_hist", "conv_hist",512, 0, 512);
        TH1D *strip = new TH1D("strip", "strip",512, 0, 512);
        
        for (int i = 1; i <= nbinsx; ++i) {

          strip->SetBinContent(i, UHist2D->GetBinContent(i,y));

        }
        convolve(strip,hist2_gausU,conv_hist);
        for (int i = 1; i <= nbinsx; ++i) { UHist2D->SetBinContent(i,y,conv_hist->GetBinContent(i));}

        delete conv_hist;
        delete strip;
      }


      for (int y = 1; y <= nbinsYv; ++y){
        TH1D *conv_hist = new TH1D("conv_hist", "conv_hist",512, 0, 512);
        TH1D *strip = new TH1D("strip", "strip",512, 0, 512);
        
        for (int i = 1; i <= nbinsx; ++i) {

          strip->SetBinContent(i, VHist2D->GetBinContent(i,y));

        }
        convolve(strip,hist2_gausV,conv_hist);
        for (int i = 1; i <= nbinsx; ++i) { VHist2D->SetBinContent(i,y,conv_hist->GetBinContent(i));}

        delete conv_hist;
        delete strip;
      }




      for (int y = 1; y <= nbinsYw; ++y){
        TH1D *conv_hist = new TH1D("conv_hist", "conv_hist",512, 0, 512);
        TH1D *strip = new TH1D("strip", "strip",512, 0, 512);
        
        for (int i = 1; i <= nbinsx; ++i) {

          strip->SetBinContent(i, WHist2D->GetBinContent(i,y));

        }
        convolve(strip,hist2_gausW,conv_hist);
        for (int i = 1; i <= nbinsx; ++i) { WHist2D->SetBinContent(i,y,conv_hist->GetBinContent(i));}

        delete conv_hist;
        delete strip;
      }



   TH1D* sum = SumHistograms(UHist2D,VHist2D,WHist2D);
        
     
/////////////////////////////////////
// SAVE STUFF
/*
TH1D * U_S = new TH1D("U_S","U_S",512,0,512);
TH1D * V_S = new TH1D("V_S","V_S",512,0,512);
TH1D * W_S = new TH1D("W_S","W_S",512,0,512);
//for(int i=1;i<=512;i++){
for(int i=1;i<=512;i++){ // go over every strip


	U_S->SetBinContent(i,UHist2D->GetBinContent(i,66));
	V_S->SetBinContent(i,VHist2D->GetBinContent(i,191));
	W_S->SetBinContent(i,WHist2D->GetBinContent(i,190));



}

TCanvas *KK = new TCanvas("KK","KK",1200,600);
KK->Divide(3,1);
KK->cd(1);
U_S->Draw();
KK->cd(2);
V_S->Draw();
KK->cd(3);
W_S->Draw();
KK->SaveAs("TH2D_MC_timeResTest.root");


*/



/////////////////////////////////////
/*
    UHist2D->SetTitle("");
    VHist2D->SetTitle("");
    WHist2D->SetTitle("");

    UHist2D->GetYaxis()->SetTickLength(0);
    VHist2D->GetYaxis()->SetTickLength(0);
    WHist2D->GetYaxis()->SetTickLength(0);
    UHist2D->GetXaxis()->SetTickLength(0);
    VHist2D->GetXaxis()->SetTickLength(0);
    WHist2D->GetXaxis()->SetTickLength(0);
    UHist2D->GetYaxis()->SetTickLength(0);
    VHist2D->GetYaxis()->SetTickLength(0);
    WHist2D->GetYaxis()->SetTickLength(0);

    UHist2D->SetStats(0);
	  UHist2D->GetXaxis()->SetTitle("");
    UHist2D->GetXaxis()->SetLabelOffset(9999);
    UHist2D->GetXaxis()->SetLabelSize(0);
    UHist2D->GetYaxis()->SetTitle("");
    UHist2D->GetYaxis()->SetLabelOffset(9999);
    UHist2D->GetYaxis()->SetLabelSize(0);
    UHist2D->SetTitle("");
	 
    VHist2D->SetStats(0); 
    VHist2D->GetXaxis()->SetTitle("");
    VHist2D->GetXaxis()->SetLabelOffset(9999);
    VHist2D->GetXaxis()->SetLabelSize(0);
    VHist2D->GetYaxis()->SetTitle("");
    VHist2D->GetYaxis()->SetLabelOffset(9999);
    VHist2D->GetYaxis()->SetLabelSize(0);
    VHist2D->SetTitle("");
		
    WHist2D->SetStats(0);
    WHist2D->GetXaxis()->SetTitle("");
    WHist2D->GetXaxis()->SetLabelOffset(9999);
    WHist2D->GetXaxis()->SetLabelSize(0);
    WHist2D->GetYaxis()->SetTitle("");
    WHist2D->GetYaxis()->SetLabelOffset(9999);
    WHist2D->GetYaxis()->SetLabelSize(0);
    WHist2D->SetTitle("");
    
    UHist2D->GetYaxis()->SetNdivisions(0);
    UHist2D->GetXaxis()->SetNdivisions(0);
    VHist2D->GetYaxis()->SetNdivisions(0);
    VHist2D->GetXaxis()->SetNdivisions(0);
    WHist2D->GetYaxis()->SetNdivisions(0);
    WHist2D->GetXaxis()->SetNdivisions(0);



    TCanvas *canvas1 = new TCanvas("canavs1","canvas1",2*257,2*269);
    canvas1->SetFrameLineColor(0);


    std::string U = "3px/U_vs_time_evt"+std::to_string(iEntry)+".png";
	  std::string V = "3px/V_vs_time_evt"+std::to_string(iEntry)+".png";
	  std::string W = "3px/W_vs_time_evt"+std::to_string(iEntry)+".png";
    */


   


  //////
sum->SetTitle("Total raw time projection");
sum->GetXaxis()->SetTitle("Time bin [arb.u.]");
sum->GetYaxis()->SetTitle("Charge/bin [arb.u.]");

UHist2D->GetZaxis()->SetLabelSize(0.03);
VHist2D->GetZaxis()->SetLabelSize(0.03);
WHist2D->GetZaxis()->SetLabelSize(0.03);

UHist2D->GetZaxis()->SetTitle("Charge/bin [arb.u.]");
VHist2D->GetZaxis()->SetTitle("Charge/bin [arb.u.]");
WHist2D->GetZaxis()->SetTitle("Charge/bin [arb.u.]");

UHist2D->SetTitleOffset(1.4, "Z");
VHist2D->SetTitleOffset(1.4, "Z");
WHist2D->SetTitleOffset(1.4, "Z");

UHist2D->SetTitleOffset(1.7, "y");
VHist2D->SetTitleOffset(1.7, "y");
WHist2D->SetTitleOffset(1.7, "y");
sum->SetTitleOffset(1.7, "y");

UHist2D->SetTitleOffset(0.9, "x");
VHist2D->SetTitleOffset(0.9, "x");
WHist2D->SetTitleOffset(0.9, "x");
sum->SetTitleOffset(0.9, "x");
sum->GetYaxis()->SetTitle("Charge/bin [arb.u.]");

TCanvas *canvas1 = new TCanvas("canvas1","canvas1",1000,1000);
canvas1->SetWindowSize(900, 900);
canvas1->SetBorderSize(0);
canvas1->Divide(2,2);

canvas1->cd(1);
UHist2D->Draw("colz 1");
gPad->SetLeftMargin(0.12);
gPad->SetRightMargin(0.2);
gPad->SetBottomMargin(0.12);
gPad->SetTopMargin(0.1);

canvas1->cd(2);
VHist2D->Draw("colz 1");
gPad->SetLeftMargin(0.12);
gPad->SetRightMargin(0.2);
gPad->SetBottomMargin(0.12);
gPad->SetTopMargin(0.1);

canvas1->cd(3);
WHist2D->Draw("colz 1");
gPad->SetLeftMargin(0.12);
gPad->SetRightMargin(0.2);
gPad->SetBottomMargin(0.12);
gPad->SetTopMargin(0.1);

canvas1->cd(4);
sum->Draw("HIST");
gPad->SetLeftMargin(0.12);

canvas1->Modified();
canvas1->Update(); 
gSystem->ProcessEvents();


/*
//gStyle->SetPadTickY(0);
// gStyle->SetPadTickX(0);

		canvas1->cd(1);//gStyle->SetPadTickY(0);
 //gStyle->SetPadTickX(0);
		UHist2D->Draw("colz");
  //  gPad->SetLeftMargin(0.12);
//gPad->SetRightMargin(0.2);
//gPad->SetBottomMargin(0.12);
//gPad->SetTopMargin(0.1);
    canvas1->Modified();
		canvas1->Update(); gSystem->ProcessEvents();
 //   canvas1->Print((U).c_str());

		canvas1->cd(2);//gStyle->SetPadTickY(0);
 //gStyle->SetPadTickX(0);
		VHist2D->Draw("colz");
  //  gPad->SetLeftMargin(0.12);
//gPad->SetRightMargin(0.2);
//gPad->SetBottomMargin(0.12);
//gPad->SetTopMargin(0.1);
    canvas1->Modified();
		canvas1->Update(); gSystem->ProcessEvents();
  //  canvas1->Print((V).c_str());


    		canvas1->cd(3);//gStyle->SetPadTickY(0);
 //gStyle->SetPadTickX(0);
		WHist2D->Draw("colz");
  //  gPad->SetLeftMargin(0.12);
//gPad->SetRightMargin(0.2);
//gPad->SetBottomMargin(0.12);
//gPad->SetTopMargin(0.1);
    canvas1->Modified();
		//canvas1->Update(); gSystem->ProcessEvents();
 //   canvas1->Print((W).c_str());

//gStyle->SetPadTickY(0);
// gStyle->SetPadTickX(0);

*/

//		canvas1->cd(4);
//    sum->Draw("HIST");
//    gPad->SetLeftMargin(0.12);

	//	WHist2D->Draw("colz ");
 ////   gPad->SetLeftMargin(0.12);
//gPad->SetRightMargin(0.2);
////gPad->SetBottomMargin(0.12);
//gPad->SetTopMargin(0.1);
   // WHist2D->GetXaxis()->SetRangeUser(5,282);
   // WHist2D->GetYaxis()->SetRangeUser(64,105);

   // canvasOTHER->Modified();
	//	canvasOTHER->Update(); gSystem->ProcessEvents();
		
  //  canvasOTHER->cd(4);
  //  sum->Draw("HIST");
  //  gPad->SetLeftMargin(0.12);


   // canvasOTHER->Modified();
	//	canvasOTHER->Update(); gSystem->ProcessEvents();

//  canvasOTHER->SaveAs("MCinspect.root");

  /////


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
  //application.Run();
  return nEntries;
}
/////////////////////////////
////////////////////////////





int reFit(const  std::string & geometryFileName, const  std::string & dataFileName) {

  // ** EVENT SOURCE ** //



  int ID[9]={ 301,923,961,1462,1462,1481,1495,2607,3200};

  
  std::shared_ptr<EventSourceMC> myEventSource = std::make_shared<EventSourceMC>(geometryFileName);
  
  // ** GEOMETRY ** //
  int index = geometryFileName.find("mbar");
  double pressure = stof(geometryFileName.substr(index-3, 3));
  IonRangeCalculator myRangeCalculator(gas_mixture_type::CO2,pressure,293.15);

  // ** MAIN LOOP ** //
  unsigned int nEntries = myEventSource->numberOfEntries();
  nEntries = 10; //TEST
  std::cout<<"pressure -> "<<pressure<<std::endl;
  std::cout<<"nEntries -> "<<nEntries<<std::endl;





  TApplication application("",nullptr,nullptr);

  for(unsigned int iEntry=0;iEntry<9;++iEntry){

    myEventSource->loadFileEntry(ID[iEntry]);    
    auto event = myEventSource->getCurrentEvent();

    std::shared_ptr<TH2D> UHist2D = event->get2DProjection(projection_type::DIR_TIME_U, filter_type::none, scale_type::raw);
    std::shared_ptr<TH2D> VHist2D = event->get2DProjection(projection_type::DIR_TIME_V, filter_type::none, scale_type::raw);
    std::shared_ptr<TH2D> WHist2D = event->get2DProjection(projection_type::DIR_TIME_W, filter_type::none, scale_type::raw);
  

    TCanvas *canvas1 = new TCanvas("canvas1","canvas1",1000,600);
    canvas1->Divide(3,1);

    canvas1->cd(1);
    UHist2D->Draw("colz 1");

    canvas1->cd(2);
    VHist2D->Draw("colz 1");

    canvas1->cd(3);
    WHist2D->Draw("colz 1");


    canvas1->Modified();
    canvas1->Update(); 
    gSystem->ProcessEvents();


    char step;
		std::cout<<"Step forward? y/n"<<std::endl;
		std::cin>>step;
		if(step=='n')return 0 ;



  
  
  
  }


  return 0;
}