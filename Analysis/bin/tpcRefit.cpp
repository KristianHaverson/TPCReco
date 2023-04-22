
#include "TFile.h"
#include "HistoManager.h"
#include "EventSourceROOT.h"
#include <boost/property_tree/json_parser.hpp>
#include <cstdlib>
#include <iostream>
#include <vector>
#include "TTree.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "TString.h"
#include "TStopwatch.h"
#include <boost/program_options.hpp>
#include "IonRangeCalculator.h"
#include "dEdxFitter.h"
#include "TrackBuilder.h"
#include "HistoManager.h"
#include "EventSourceROOT.h"
#ifdef WITH_GET
#include "EventSourceGRAW.h"
#include "EventSourceMultiGRAW.h"
#endif
#include "SigClusterTPC.h"
#include "EventTPC.h"
#include "RecoOutput.h"
#include "RunIdParser.h"
#include "InputFileHelper.h"
#include "MakeUniqueName.h"
#include "colorText.h"
#include <stdio.h>
#include <TROOT.h>
#include <TF1.h>
#include <TH2.h>
#include <TGraph.h>
#include <TLatex.h>
#include <TGraph.h>
#include <TSpectrum.h>
#include <TMath.h>
#include <tuple>
#include <TFile.h>
#include <TTree.h>
#include <TROOT.h>
#include <TCanvas.h>
#include <TObject.h>
#include <TString.h>
#include <iostream>
#include <sstream>
#include <TApplication.h>
#include <TSystem.h>
#include <TStyle.h>
#include <fstream>
#include <iterator>
#include <vector>
#include <HIGS_trees_analysis.h>
#include "LinkDef.h"


#define BOLDRED     "\033[1m\033[31m"       /* Bold Red */
#define BOLDGREEN   "\033[1m\033[32m"       /* Bold Green */
#define BOLDYELLOW  "\033[1m\033[33m"       /* Bold Yellow */
#define BOLDBLUE    "\033[1m\033[34m"       /* Bold Blue */
#define BOLDMAGENTA "\033[1m\033[35m"       /* Bold Magenta */
#define BOLDCYAN    "\033[1m\033[36m"       /* Bold Cyan */
#define BOLDWHITE   "\033[1m\033[37m"       /* Bold White */


boost::program_options::variables_map parseCmdLineArgs(int argc, char **argv){

	boost::program_options::options_description cmdLineOptDesc("Allowed options");
	cmdLineOptDesc.add_options()
		("help", "produce help message")
		("dataFile",  boost::program_options::value<std::string>(), "string - path to data file (OFFLINE) or directory (ONLINE). Overrides the value from the config file. In multi-GRAW mode specify several files separated by commas.")
		("dataFile",  boost::program_options::value<std::string>(), "string - path to data file (OFFLINE) or directory (ONLINE). Overrides the value from the config file. In multi-GRAW mode specify several files separated by commas.")
		("singleAsadGrawFile", boost::program_options::value<bool>(), "bool - Flag enabling multi-GRAW mode. One file stream per each AsAd board.")
		("removePedestal",  boost::program_options::value<bool>(), "bool - Flag to control pedestal removal. Overrides the value from config file.")
		("recoClusterEnable",  boost::program_options::value<bool>(), "bool - Flag to enable RECO cluster.")
		("recoClusterThreshold",  boost::program_options::value<double>(), "double - ADC threshold above pedestal for RECO cluster.")
		("recoClusterDeltaStrips",  boost::program_options::value<int>(), "int - Envelope in strip units around seed hits for RECO cluster.")
		("recoClusterDeltaTimeCells",  boost::program_options::value<int>(), "int - Envelope in time cell units around seed hits for RECO cluster.");

	boost::program_options::variables_map varMap;        
	boost::program_options::store(boost::program_options::parse_command_line(argc, argv, cmdLineOptDesc), varMap);
	boost::program_options::notify(varMap); 

	if (varMap.count("help")) {
		std::cout<<cmdLineOptDesc<<std::endl;
		exit(1);
	}
	return varMap;
}



// oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo //
// oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo //
std::vector<std::vector<double>> read_csv(const std::string& filename){
    std::ifstream file(filename);
    std::vector<std::vector<double>> data;
    std::string line;
    while (std::getline(file, line)){
        std::vector<double> row;
        std::string value;
        for (char ch : line){
            if (ch == ','){
                row.push_back(std::stod(value));
                value.clear();
            }else{
                value += ch;
            }
        }
        if (!value.empty()){
            row.push_back(std::stod(value));
        }
        data.push_back(row);
    }
    return data;
}


// oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo //
// oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo //

std::tuple<double, double, double, double> get2DLine(  TVector3 aStart,  TVector3 aEnd, int aDir, std::shared_ptr<GeometryTPC> myGeometryPtr){


  bool err_flag = false;
  int startX = myGeometryPtr->Pos2timecell(aStart.Z(), err_flag);
  int endX = myGeometryPtr->Pos2timecell(aEnd.Z(), err_flag);

  double directionScale = 1.0/myGeometryPtr->GetStripPitch();
  TVector2 referencePoint = myGeometryPtr->GetReferencePoint();
  TVector2 start2D(aStart.X(), aStart.Y());
  TVector2 end2D(aEnd.X(), aEnd.Y());
  
  start2D -= referencePoint;
  end2D -= referencePoint;
  
  int offset = 1;
  if(aDir==0){
    start2D += 2*referencePoint;
    end2D += 2*referencePoint;
    offset = 1;
  }
  if(aDir==2){
    start2D += 4*referencePoint;
    end2D += 4*referencePoint;
    offset = 6;
  }
  
  int startY = myGeometryPtr->Cartesian2posUVW(start2D, aDir, err_flag)*directionScale + offset;
  int endY = myGeometryPtr->Cartesian2posUVW(end2D, aDir, err_flag)*directionScale + offset;
	      



  return std::make_tuple(startX, startY, endX, endY);
}


// oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo //
// oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo //

int reFit(std::unique_ptr<EventSourceGRAW>& myEventSource, std::string geomFileName, int index,
		TVector3 a1_end ,TVector3 a2_end , TVector3 a3_end , TVector3 vertexPos  ) 
	{

	std::shared_ptr<GeometryTPC> myGeometryPtr = std::make_shared<GeometryTPC>(geomFileName.c_str());


    TLine uLine1 ,uLine2,uLine3;
    TLine vLine1 ,vLine2,vLine3;
    TLine wLine1 ,wLine2,wLine3;


	double U_X_a1_end,  U_X_a2_end,  U_X_a3_end;
	double U_Y_a1_end,  U_Y_a2_end,  U_Y_a3_end;

	double V_X_a1_end,  V_X_a2_end,  V_X_a3_end;
	double V_Y_a1_end,  V_Y_a2_end,  V_Y_a3_end;

	double W_X_a1_end,  W_X_a2_end,  W_X_a3_end;
	double W_Y_a1_end,  W_Y_a2_end,  W_Y_a3_end;
	
	double U_X_vert;  
	double U_Y_vert;  

	double V_X_vert;  	
	double V_Y_vert;  
	
	double W_X_vert;  
	double W_Y_vert;  

	myEventSource->loadFileEntry(index);    
    auto event = myEventSource->getCurrentEvent();


    std::shared_ptr<TH2D> UHist2D = event->get2DProjection(projection_type::DIR_TIME_U, filter_type::island, scale_type::raw);
    std::shared_ptr<TH2D> VHist2D = event->get2DProjection(projection_type::DIR_TIME_V, filter_type::island, scale_type::raw);
    std::shared_ptr<TH2D> WHist2D = event->get2DProjection(projection_type::DIR_TIME_W, filter_type::island, scale_type::raw);

    int binsY_U = UHist2D->GetNbinsY();
    int binsY_V = VHist2D->GetNbinsY();
    int binsY_W = WHist2D->GetNbinsY();
    std::cout<<binsY_U<<std::endl;
    std::cout<<binsY_V<<std::endl;
    std::cout<<binsY_W<<std::endl;

    std::shared_ptr<TH2D> UHist2D_2 = event->get2DProjection(projection_type::DIR_TIME_U, filter_type::island, scale_type::mm);
    std::shared_ptr<TH2D> VHist2D_2 = event->get2DProjection(projection_type::DIR_TIME_V, filter_type::island, scale_type::mm);
    std::shared_ptr<TH2D> WHist2D_2 = event->get2DProjection(projection_type::DIR_TIME_W, filter_type::island, scale_type::mm);

    //oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo//
    //oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo//

   // int nbinsU = UHist2D->GetNbinsY();
 //   int nbinsV = VHist2D->GetNbinsY();
   // i//nt nbinsW = 1*WHist2D->GetNbinsY();

    TCanvas *canvas1 = new TCanvas("canvas1","canvas1",1200,800);
    canvas1->Divide(3,2);
	UHist2D->SetStats(0);
	VHist2D->SetStats(0);
	WHist2D->SetStats(0);

    canvas1->cd(1);
    UHist2D->Draw("colz2 0");

	

	
	std::tie(U_X_vert, U_Y_vert, U_X_a1_end, U_Y_a1_end) = get2DLine(vertexPos,a1_end, 0, myGeometryPtr);
	std::tie(U_X_vert, U_Y_vert, U_X_a2_end, U_Y_a2_end) = get2DLine(vertexPos,a2_end, 0, myGeometryPtr);
	std::tie(U_X_vert, U_Y_vert, U_X_a3_end, U_Y_a3_end) = get2DLine(vertexPos,a3_end, 0, myGeometryPtr);
	uLine1.DrawLine( U_X_vert,   U_Y_vert,  U_X_a1_end,  U_Y_a1_end);
	uLine2.DrawLine( U_X_vert,   U_Y_vert,  U_X_a2_end,  U_Y_a2_end);
	uLine3.DrawLine( U_X_vert,   U_Y_vert,  U_X_a3_end,  U_Y_a3_end);

    // std::cout<<U_X_vert<<" , "<<U_Y_vert<<"   "<<U_X_a1_end<<" ,  "<<U_Y_a1_end<<std::endl;
    // std::cout<<U_X_vert<<" , "<<U_Y_vert<<"   "<<U_X_a2_end<<" ,  "<<U_Y_a2_end<<std::endl;
    // std::cout<<U_X_vert<<" , "<<U_Y_vert<<"   "<<U_X_a3_end<<" ,  "<<U_Y_a3_end<<std::endl;
 
    canvas1->cd(2);
    VHist2D->Draw("colz2 0");
	std::tie(V_X_vert, V_Y_vert, V_X_a1_end, V_Y_a1_end) = get2DLine(vertexPos,a1_end, 1, myGeometryPtr);
	std::tie(V_X_vert, V_Y_vert, V_X_a2_end, V_Y_a2_end) = get2DLine(vertexPos,a2_end, 1, myGeometryPtr);
	std::tie(V_X_vert, V_Y_vert, V_X_a3_end, V_Y_a3_end) = get2DLine(vertexPos,a3_end, 1, myGeometryPtr);
	vLine1.DrawLine(V_X_vert, V_Y_vert, V_X_a1_end, V_Y_a1_end);
	vLine2.DrawLine(V_X_vert, V_Y_vert, V_X_a2_end, V_Y_a2_end);
	vLine3.DrawLine(V_X_vert, V_Y_vert, V_X_a3_end, V_Y_a3_end);


    canvas1->cd(3);
    WHist2D->Draw("colz2 0");
	std::tie(W_X_vert, W_Y_vert, W_X_a1_end, W_Y_a1_end) = get2DLine(vertexPos,a1_end, 2, myGeometryPtr);
	std::tie(W_X_vert, W_Y_vert, W_X_a2_end, W_Y_a2_end) = get2DLine(vertexPos,a2_end, 2, myGeometryPtr);
	std::tie(W_X_vert, W_Y_vert, W_X_a3_end, W_Y_a3_end) = get2DLine(vertexPos,a3_end, 2, myGeometryPtr);
	wLine1.DrawLine(W_X_vert, (W_Y_vert)  ,  W_X_a1_end, (W_Y_a1_end) );
	wLine2.DrawLine(W_X_vert, (W_Y_vert)  ,  W_X_a2_end, (W_Y_a2_end) );
	wLine3.DrawLine(W_X_vert, (W_Y_vert)  ,  W_X_a3_end, (W_Y_a3_end) );
	//oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo//
	//oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo//

    canvas1->cd(4);
    UHist2D->Draw("colz2 0");
	uLine1.DrawLine( U_X_vert,   U_Y_vert,  U_X_a1_end,  U_Y_a1_end);
	uLine2.DrawLine( U_X_vert,   U_Y_vert,  U_X_a2_end,  U_Y_a2_end);
	uLine3.DrawLine( U_X_vert,   U_Y_vert,  U_X_a3_end,  U_Y_a3_end);

    canvas1->cd(5);
    VHist2D->Draw("colz2 0");
	vLine1.DrawLine(V_X_vert, V_Y_vert, V_X_a1_end, V_Y_a1_end);
	vLine2.DrawLine(V_X_vert, V_Y_vert, V_X_a2_end, V_Y_a2_end);
	vLine3.DrawLine(V_X_vert, V_Y_vert, V_X_a3_end, V_Y_a3_end);


    canvas1->cd(6);
    WHist2D->Draw("colz2 0");
	wLine1.DrawLine(W_X_vert, (W_Y_vert)  ,  W_X_a1_end, (W_Y_a1_end) );
	wLine2.DrawLine(W_X_vert, (W_Y_vert)  ,  W_X_a2_end, (W_Y_a2_end) );
	wLine3.DrawLine(W_X_vert, (W_Y_vert)  ,  W_X_a3_end, (W_Y_a3_end) );



    canvas1->Modified();
    canvas1->Update(); 
    gSystem->ProcessEvents();

	//oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo//

    char step;
	std::cout<<"Step forward? y/n"<<std::endl;
	std::cin>>step;
	if(step=='n')return 0 ;



  
  


  return 0;
}

// oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo //
// oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo //
// oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo //
// oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo //

/*
int main(int argc, char **argv){


	//==========================================================================================//
	// ** FETCH FROM JSON ** //
	boost::program_options::variables_map varMap = parseCmdLineArgs(argc, argv);
	boost::property_tree::ptree tree;
	if(argc<1){
		std::cout<<" Usage: tpcGUI config.json"<<std::endl;
		return 0;
	}else {
		std::cout<<"Using configFileName: "<<argv[1]<<std::endl;
		boost::property_tree::read_json(argv[1], tree);
	}
	if (varMap.count("dataFile")) {
		tree.put("dataFile",varMap["dataFile"].as<std::string>());
	}
	if (varMap.count("removePedestal")) {
		tree.put("removePedestal",varMap["removePedestal"].as<bool>());
	}
	if (varMap.count("geometryFile")) {
		tree.put("geometryFile",varMap["geometryFile"].as<std::string>());
	}
	 
	//==========================================================================================//
	// ** GET STUFF ** //
	TApplication application("",nullptr,nullptr);


	std::string dataFileName = tree.get("dataFile","");
	std::string geometryFileName = tree.get("geometryFile","");

	std::unique_ptr<EventSourceGRAW> myEventSource = std::make_unique<EventSourceMultiGRAW>(geometryFileName);
	
	myEventSource->loadDataFile(dataFileName);
	myEventSource->setRemovePedestal(true);
	myEventSource->configurePedestal(tree.find("pedestal")->second);
	myEventSource->getEventFilter().setConditions(tree);
	myEventSource->loadFileEntry(0);

	//==========================================================================================//

    std::string inFile = "/home/kris/Desktop/ELI_fitter/prong3_out.csv";
    auto data = read_csv(inFile); 
    int size = data.size();
    std::cout << "entries in file -> " << size << std::endl;
    
    int graw, index, RunID;
    int eventID;
	
	std::ostringstream oss;
	
	bool DEBUG = false;
	std::string fileKey = "0036";
    for (const auto& col : data){
		
        graw = col[0];
        index = col[1];
        eventID = col[2];
        RunID  = col[3];

		TVector3    a1_end(col[4],col[5],col[6]);
		TVector3    a2_end(col[7],col[8],col[9]);
		TVector3    a3_end(col[10],col[11],col[12]);
		TVector3 vertexPos( col[13], col[14], col[15]);

		
		std::string GRAW_Str;
    	oss << std::setw(4) << std::setfill('0') << graw;
    	GRAW_Str = oss.str();
		oss.str("");


		if(DEBUG==true){
			std::cout<<graw<<std::endl;
			std::cout<<index<<std::endl;
			std::cout<<eventID<<std::endl;
			std::cout<<std::fixed << std::setprecision(0)<< RunID<<std::endl;
			std::cout<<std::endl;
		}


		if(GRAW_Str!=fileKey)std::cout<<"Searching"<<std::endl;
		else{
			reFit(myEventSource,geometryFileName,index , a1_end,a2_end,a3_end,vertexPos);
		}


	
  	}






	return 0;


}

*/

int main(int argc, char **argv){
	/*
  	const double CONV_Kev2MeV = 1./1000.;
	const double Q_VALUE_12C  = - 7367.0 * CONV_Kev2MeV; 

	const double beamE = 13.1;// +0.12; // [MeV]
	double Ecm = beamE + Q_VALUE_12C;// - (3);


	const double AMU        = 931.49410242;
 	const double MASS_8Be   = AMU  * 8.00530510000; 
	const double MASS_4He   = AMU  * 4.00260325413; 

	const double KE_alpha0  = Ecm * (MASS_8Be  /(MASS_4He+MASS_8Be));
	const double KE_BE =      Ecm * (MASS_4He  /(MASS_4He+MASS_8Be));


	
	IonRangeCalculator myRangeCalculator(gas_mixture_type::CO2, 250 ,293.15);
	double ionRangeMM  = myRangeCalculator.getIonRangeMM( pid_type::ALPHA ,KE_alpha0);

	std::cout<<"==================="<<std::endl;
	std::cout<<"beam  E     -> "<< beamE<<std::endl;
	std::cout<<"Q value     -> "<<Q_VALUE_12C <<std::endl;
	std::cout<<"CMS   E     -> "<<Ecm <<std::endl;
	std::cout<<"KE_alpha0   -> "<<KE_alpha0 <<std::endl;
	std::cout<<"KE_BE       -> "<<KE_BE <<std::endl;
	std::cout<<"==================="<<std::endl;
	std::cout<<"ion range a0    -> "<< ionRangeMM<<std::endl;
	std::cout<<"==================="<<std::endl;
*/

	double aLen1  = 10.7535  ;
	double aLen2  = 43.936 ;
	double aLen3  = 30.6014 ;
	IonRangeCalculator myRangeCalculator(gas_mixture_type::CO2, 250 ,293.15);

	double E = myRangeCalculator.getIonEnergyMeV( pid_type::ALPHA ,aLen1 );
	E += myRangeCalculator.getIonEnergyMeV( pid_type::ALPHA ,aLen2 );
	E += myRangeCalculator.getIonEnergyMeV( pid_type::ALPHA ,aLen3 );
	std::cout << "Total E -> "<< E <<std::endl;
	std::cout << "Total D -> "<<aLen1+aLen2+aLen3 <<std::endl;
	
}