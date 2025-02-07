#include "dEdxFitter.h"
#include "colorText.h"

#include "TFitResultPtr.h"
#include "Math/MinimizerOptions.h"
#include <TMath.h>
#include "TRandom3.h"

TGraph* dEdxFitter::braggGraph_alpha = new TGraph("dEdx_corr_alpha_10MeV_CO2_250mbar.dat", "%lg %lg");
TGraph* dEdxFitter::braggGraph_12C = new TGraph("dEdx_corr_12C_5MeV_CO2_250mbar.dat", "%lg %lg");
double dEdxFitter::nominalPressure = 250.0;

double dEdxFitter::currentPressure = 190.0;
////////////////////////////////////////////////
////////////////////////////////////////////////
dEdxFitter::dEdxFitter(double aPressure){

  setPressure(aPressure);

  braggGraph_alpha->SetBit(TGraph::kIsSortedX);
  braggGraph_12C->SetBit(TGraph::kIsSortedX);

  alpha_ionisation = new TF1("alpha_ionisation", bragg_alpha, 0, 600.0);
  carbon_ionisation = new TF1("carbon_ionisation",bragg_12C, 0,  200.0);

  carbon_alpha_model = new TF1("carbon_alpha_model",bragg_12C_alpha, -20, 350.0, 7);

  carbon_alpha_model->SetParName(0, "sigma");
  carbon_alpha_model->SetParName(1, "vertexOffset");
  carbon_alpha_model->SetParName(2, "alphaOffset");
  carbon_alpha_model->SetParName(3, "carbonOffset");
  carbon_alpha_model->SetParName(4, "alphaScale");
  carbon_alpha_model->SetParName(5, "carbonScale");
  carbon_alpha_model->SetParName(6, "commonScale");

  carbon_alpha_model->SetParLimits(0, 0.1, 3.0);
  carbon_alpha_model->SetParLimits(6, 10, 1000);
  
  carbon_alpha_model->FixParameter(4, 1.0);
  carbon_alpha_model->FixParameter(5, 1.0);  
 
  alpha_model = new TF1(*carbon_alpha_model);
  alpha_model->SetName("alpha_model");
  alpha_model->FixParameter(1, 0.0);
  alpha_model->FixParameter(3, 0.0);
  alpha_model->FixParameter(4, 0.0);
  alpha_model->FixParameter(5, 0.0);
 
  reset();
}
////////////////////////////////////////////////
////////////////////////////////////////////////
void dEdxFitter::setPressure(double aPressure) {
  
  currentPressure = aPressure;

  minVtxOffset = -5;
  minAlphaOffset = 0;
  minCarbonOffset = 0;

  //190 mbar, 10 MeV alpha, 2.5 MeV C
  if(std::abs(nominalPressure-190)<1E-3){
    maxAlphaOffset = (385.99)*(nominalPressure/currentPressure); 
    maxCarbonOffset = (28.64)*(nominalPressure/currentPressure);
  }
  //250 mbar, 10 MeV alpha, 5 MeV C
  else if(std::abs(nominalPressure-250)<1E-3){
    maxAlphaOffset = (297.23)*(nominalPressure/currentPressure);
    maxCarbonOffset = (23.43)*(nominalPressure/currentPressure);
  }
  else{
    std::cout<<KRED<<"dEdxFitter: nominal pressure: "<<RST<<nominalPressure
	     <<" does not corresond to any dEdx data files."<<std::endl;
    exit(0);
  }
}
////////////////////////////////////////////////
////////////////////////////////////////////////
void dEdxFitter::reset(){

  carbon_alpha_model->SetRange(-20, maxAlphaOffset+maxCarbonOffset);
  carbon_alpha_model->SetParLimits(1, 0.0, maxVtxOffset);
  
  carbon_alpha_model->SetParLimits(2, minAlphaOffset, maxAlphaOffset);
  carbon_alpha_model->SetParLimits(3, minCarbonOffset, maxCarbonOffset);
  carbon_alpha_model->SetParameters(1.0,
				    (minVtxOffset+maxVtxOffset)/2.0,
				    (minAlphaOffset+maxAlphaOffset)/2.0,
				    (minCarbonOffset+maxCarbonOffset)/2.0,
				    1.0, 1.0, 0.006);

  alpha_model->SetRange(-20, maxAlphaOffset);
  alpha_model->SetParLimits(1, minVtxOffset, maxVtxOffset);
  alpha_model->SetParLimits(2, minAlphaOffset, maxAlphaOffset);
  alpha_model->SetParameters(1.0,
			     (minVtxOffset+maxVtxOffset)/2.0,
			     (minAlphaOffset+maxAlphaOffset)/2.0,
			     0.0,
			     1.0, 0.0, 0.006);
  
  theFitResult = TFitResult();
  theFittedModel = alpha_model;
  theFittedHisto = emptyHisto;
  bestFitEventType = pid_type::UNKNOWN;
}
////////////////////////////////////////////////
////////////////////////////////////////////////
double dEdxFitter::bragg_alpha(double *x, double *params) {
  
  double sigma = params[0];
  double vertex_pos = params[1];
  double shift = params[2];  
  double n_sigma = 3.0;

  double pressure_scale_factor = currentPressure/nominalPressure;
  double x_shifted = (x[0] - vertex_pos) + shift;
  double x_scaled = pressure_scale_factor*x_shifted;

  double total_dEdx_length = 300;
  double a = pressure_scale_factor*shift;
  double b = total_dEdx_length;
  double t = (x[0] - vertex_pos)*pressure_scale_factor + a;

  double bragg_at_n_sigma = braggGraph_alpha->Eval(a+n_sigma*sigma*pressure_scale_factor);
  double p0 = bragg_at_n_sigma;
  double p1 = 0.0;
  double p2 = 0.0;

  double smeared_edge = 0;
  smeared_edge += 0.5*TMath::Erf((t-a)/sqrt(2)/sigma)*(p2*(t*t + sigma*sigma) + p1*t + p0);
  smeared_edge += 1.0/sqrt(2*M_PI)*sigma*exp(-(a-t)*(a-t)/(2*sigma*sigma))*(p2*(a+t) + p1);
  smeared_edge -= 0.5*TMath::Erf((t-b)/sqrt(2)/sigma)*(p2*(t*t + sigma*sigma) + p1*t + p0);
  smeared_edge -= 1.0/sqrt(2*M_PI)*sigma*exp(-(b-t)*(b-t)/(2*sigma*sigma))*(p2*(b+t) + p1);

  double bragg = braggGraph_alpha->Eval(x_scaled)*(x_scaled<total_dEdx_length)*(x_scaled>0)*(x_shifted-shift>0);
  return smeared_edge*(x_shifted-shift<n_sigma*sigma) + bragg*(x_shifted-shift>n_sigma*sigma);
}
////////////////////////////////////////////////
////////////////////////////////////////////////
double dEdxFitter::bragg_12C(double *x, double *params) {

  double sigma = params[0];
  double vertex_pos = params[1];
  double shift = params[3];  
  double n_sigma = 20.0;

  double pressure_scale_factor = currentPressure/nominalPressure;
  double x_shifted = (vertex_pos-x[0]) + shift;
  double x_scaled = pressure_scale_factor*x_shifted;

  double total_dEdx_length = 23.5;
  double a = pressure_scale_factor*shift;
  double b = total_dEdx_length;
  double t = (vertex_pos - x[0])*pressure_scale_factor + a;

  double p0 = 0.0, p1=0.0, p2 = 0.0;
  if(t<15){
    p0 = 254.284;
    p1 = 16.6512;
    p2 = -1.40697;
  }
  else if(t>15){
    p0 = 533.018;
    p1 = -23.1911;
    p2 = 0.00861351;
  }

  double smeared_edge = 0;
  smeared_edge += 0.5*TMath::Erf((t-a)/sqrt(2)/sigma)*(p2*(t*t + sigma*sigma) + p1*t + p0);
  smeared_edge += 1.0/sqrt(2*M_PI)*sigma*exp(-(a-t)*(a-t)/(2*sigma*sigma))*(p2*(a+t) + p1);
  smeared_edge -= 0.5*TMath::Erf((t-b)/sqrt(2)/sigma)*(p2*(t*t + sigma*sigma) + p1*t + p0);
  smeared_edge -= 1.0/sqrt(2*M_PI)*sigma*exp(-(b-t)*(b-t)/(2*sigma*sigma))*(p2*(b+t) + p1);

  double bragg = braggGraph_12C->Eval(x_scaled)*(x_scaled<total_dEdx_length)*(x_scaled>0)*(x_shifted-shift>0);

  return smeared_edge*(x_shifted-shift<n_sigma*sigma) + bragg*(x_shifted-shift>n_sigma*sigma);
}
////////////////////////////////////////////////
////////////////////////////////////////////////
double dEdxFitter::bragg_12C_alpha(double *x, double *params) {

  double value = 0.0;
  double alpha_scale = params[4];
  double carbon_scale = params[5];
  double common_scale = params[6];
  
  value = alpha_scale*bragg_alpha(x, params);
  value += carbon_scale*bragg_12C(x, params);
                 
  return common_scale*value;
}
////////////////////////////////////////////////
////////////////////////////////////////////////
TH1F dEdxFitter::reflectHisto(const TH1F &aHisto) const{

  TH1F hReflected(aHisto);
  hReflected.Reset();
 
  double value=0.0;
  for(int iBin=0;iBin<=aHisto.GetNbinsX()+1;++iBin){
    value = aHisto.GetBinContent(iBin);
    hReflected.SetBinContent(aHisto.GetNbinsX()-iBin, value);
    value = aHisto.GetBinError(iBin);
    hReflected.SetBinError(aHisto.GetNbinsX()-iBin, value);
    
  }
  return hReflected;
}
////////////////////////////////////////////////
////////////////////////////////////////////////
TFitResult dEdxFitter::fitHypothesis(TF1 *fModel, TH1F & aHisto){

  TFitResult theResult;
  if(!aHisto.GetEntries()) return theResult;

  double tkLength = aHisto.GetXaxis()->GetXmax();
  maxVtxOffset = maxCarbonOffset;
  minAlphaOffset = std::max(0.0, maxAlphaOffset - tkLength);
  minCarbonOffset = std::max(0.0, maxCarbonOffset - tkLength);

  int fitCounter = 0;
  TFitResultPtr theResultPtr;
  do{
    reset();
    theResultPtr = aHisto.Fit(fModel,"BRWWS");
    if(!theResultPtr.Get()) break;
    ++fitCounter;
  }while(!theResultPtr->IsValid() && fitCounter<3);
  
  if(theResultPtr.Get()) theResult = *theResultPtr.Get();
  else{
    std::cout<<KRED<<"No fit result"<<RST<<std::endl;
  }
  return theResult;
}
////////////////////////////////////////////////
////////////////////////////////////////////////
TFitResult dEdxFitter::fitHisto(const TH1F & aHisto){

  ////////////C12+alpha hypothesis
  TH1F fittedHisto_for_C12_alpha = aHisto;
  bool reflection_for_C12_alpha = false;
  int maxBin = aHisto.GetMaximumBin();
  if(maxBin>aHisto.GetNbinsX()/2.0){
    fittedHisto_for_C12_alpha = reflectHisto(aHisto);
    reflection_for_C12_alpha = true;
  }
  TFitResult carbon_alphaResult = fitHypothesis(carbon_alpha_model, fittedHisto_for_C12_alpha);

  ////////////alpha hypothesis
  bool reflection_for_alpha = false;
  TH1F fittedHisto_for_alpha = aHisto;
  maxBin = aHisto.GetMaximumBin();
  if(maxBin<aHisto.GetNbinsX()/2.0){
    fittedHisto_for_alpha = reflectHisto(aHisto);
    reflection_for_alpha = true;
  }
  TFitResult alphaResult = fitHypothesis(alpha_model, fittedHisto_for_alpha);

  if(alphaResult.MinFcnValue()<carbon_alphaResult.MinFcnValue()){
    theFitResult = alphaResult;
    theFittedModel = alpha_model;
    theFittedHisto = fittedHisto_for_alpha;
    isReflected = reflection_for_alpha;
    bestFitEventType = pid_type::ALPHA;
  }
  else{
    theFitResult = carbon_alphaResult;
    theFittedModel = carbon_alpha_model;
    theFittedHisto = fittedHisto_for_C12_alpha;
    isReflected = reflection_for_C12_alpha;
    bestFitEventType = pid_type::C12_ALPHA;
  }
  theFittedModel->SetParameters(theFitResult.Parameters().data());
  return theFitResult;
}
////////////////////////////////////////////////
////////////////////////////////////////////////
double dEdxFitter::getVertexOffset() const{

  return theFittedModel->GetParameter("vertexOffset");  
}
////////////////////////////////////////////////
////////////////////////////////////////////////
double dEdxFitter::getAlphaRange() const{

  if(bestFitEventType==pid_type::UNKNOWN) return 0.0;
    
  double alphaOffset = theFittedModel->GetParameter("alphaOffset");
  double alphaRange = maxAlphaOffset - alphaOffset;
  return alphaRange;
}
////////////////////////////////////////////////
////////////////////////////////////////////////
double dEdxFitter::getCarbonRange() const{

  if(bestFitEventType!=pid_type::C12_ALPHA) return 0.0;

  double carbonOffset = theFittedModel->GetParameter("carbonOffset");
  double carbonRange = maxCarbonOffset - carbonOffset;
  return carbonRange;
}
////////////////////////////////////////////////
////////////////////////////////////////////////
