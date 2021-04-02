#include <cstdlib>
#include <iostream>

#include <TGResourcePool.h>
#include <TGFrame.h>
#include <TGLabel.h>
#include <TGTextEntry.h>
#include <TGNumberEntry.h>
#include <TGCanvas.h>
#include <TGTableLayout.h>
#include <TGFontDialog.h>
#include <TFrame.h>

#include "colorText.h"
#include "MarkersManager.h"
#include "EntryDialog.h"
#include "MainFrame.h"
#include "HistoManager.h"
#include "ScrollFrame.h"
#include "CommonDefinitions.h"

/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
MarkersManager::MarkersManager(const TGWindow * p, MainFrame * aFrame)
 : TGCompositeFrame(p, 10, 10, kVerticalFrame), fParentFrame(aFrame){

   SetCleanup(kDeepCleanup);

   fTopFrame = new TGVerticalFrame(this, 300, 300);
   TGLayoutHints aLayoutHints(kLHintsTop | kLHintsLeft | kLHintsExpandY |
			      kLHintsShrinkX|kLHintsShrinkY |
			      kLHintsFillX|kLHintsFillY, 2, 2, 2, 2);
   AddFrame(fTopFrame, &aLayoutHints);
   
   TGGroupFrame *aHeaderFrame = new TGGroupFrame(fTopFrame, "Segment selection");
   fTopFrame->AddFrame(aHeaderFrame, new TGLayoutHints(kLHintsExpandX, 2, 2, 1, 1));

   TGTextButton* aButton = new TGTextButton(aHeaderFrame,"Add", M_ADD_SEGMENT);
   ULong_t aColor;
   gClient->GetColorByName("yellow", aColor);
   aButton->ChangeBackground(aColor);
   
   aHeaderFrame->AddFrame(aButton, new TGLayoutHints(kLHintsLeft, 2, 2, 1, 1));
   aButton->Connect("Clicked()","MarkersManager",this,"DoButton()");
   /*
   TGLabel *aLabel = new TGLabel(aHeaderFrame,"U");
   aHeaderFrame->AddFrame(aLabel, new TGLayoutHints(kLHintsLeft, 2, 2, 1, 1));

   aLabel = new TGLabel(aHeaderFrame,"V");
   aHeaderFrame->AddFrame(aLabel, new TGLayoutHints(kLHintsLeft, 2, 2, 1, 1));

   aLabel = new TGLabel(aHeaderFrame,"W");
   aHeaderFrame->AddFrame(aLabel, new TGLayoutHints(kLHintsLeft, 2, 2, 1, 1));
   */
   /*
   fMarkerGCanvas = new TGCanvas(fTopFrame, 300, 300);
   TGCompositeFrame *aMarkerContainer = new TGCompositeFrame(fMarkerGCanvas->GetViewPort(), kVerticalFrame);
   fMarkerGCanvas->SetContainer(aMarkerContainer);
   fTopFrame->AddFrame(fMarkerGCanvas, new TGLayoutHints(kLHintsExpandX, 2, 2, 1, 1));
   fTopFrame->Layout();
   */
   /*
   for(int iMarker=0;iMarker<4;++iMarker){
     addMarkerFrame(iMarker);
   }
   */

   initialize();
}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
MarkersManager::~MarkersManager(){

  delete fMarkerGCanvas;
}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
void MarkersManager::initialize(){

  fMarkersContainer.resize(3);
  fLinesContainer.resize(3);
  
  firstMarker = 0;
  secondMarker = 0;

}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
void MarkersManager::reset(){

   std::for_each(fMarkersContainer.begin(), fMarkersContainer.end(),
		 [](TMarker *item){if(item){delete item; item = 0;}});

  if(firstMarker) delete firstMarker;
  if(secondMarker) delete secondMarker;
  firstMarker = 0;
  secondMarker = 0;

  clearLines();
}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
void MarkersManager::clearLines(){

std::for_each(fLinesContainer.begin(), fLinesContainer.end(),
	      [](TLine *item){if(item){delete item; item = 0;}});
  
}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
/*
void MarkersManager::addHeaderFrame(){
TGHorizontalFrame *aHorizontalFrame = new TGHorizontalFrame(fMarkerGCanvas->GetContainer(), 200, 30);
  
}
*/
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
void MarkersManager::addMarkerFrame(int iMarker){

  TGHorizontalFrame *aHorizontalFrame = new TGHorizontalFrame(fMarkerGCanvas->GetContainer(), 200, 30);
  TGLayoutHints *aLayoutHints = new TGLayoutHints(kLHintsLeft, 2, 2, 2, 2);
  TGCompositeFrame *aCompositeFrame = (TGCompositeFrame*)fMarkerGCanvas->GetContainer();
  aCompositeFrame->AddFrame(aHorizontalFrame, aLayoutHints);

  float value = 1.0;
  TGNumberEntry *aNumberEntry = new TGNumberEntry(aHorizontalFrame, value, 5, 0,
						  TGNumberFormat::EStyle::kNESRealTwo);
  aNumberEntry->Connect("ValueSet(Long_t)","MainFrame",fParentFrame,"ProcessMessage(Long_t)");
  aNumberEntry->Associate(this);
  aHorizontalFrame->AddFrame(aNumberEntry, aLayoutHints);

  aNumberEntry = new TGNumberEntry(aHorizontalFrame, value, 5, 0,
				   TGNumberFormat::EStyle::kNESRealTwo);
  aNumberEntry->Connect("ValueSet(Long_t)","MainFrame",fParentFrame,"ProcessMessage(Long_t)");
  aNumberEntry->Associate(this);
  aHorizontalFrame->AddFrame(aNumberEntry, aLayoutHints);

  aNumberEntry = new TGNumberEntry(aHorizontalFrame, value, 5, 0,
				   TGNumberFormat::EStyle::kNESRealTwo);
  aNumberEntry->Connect("ValueSet(Long_t)","MainFrame",fParentFrame,"ProcessMessage(Long_t)");
  aNumberEntry->Associate(this);
  aHorizontalFrame->AddFrame(aNumberEntry, aLayoutHints);

  aNumberEntry = new TGNumberEntry(aHorizontalFrame, value, 5, 0,
				   TGNumberFormat::EStyle::kNESRealTwo);
  aNumberEntry->Connect("ValueSet(Long_t)","MainFrame",fParentFrame,"ProcessMessage(Long_t)");
  aNumberEntry->Associate(this);
  aHorizontalFrame->AddFrame(aNumberEntry, aLayoutHints);
 
  fMarkerGCanvas->Layout();
}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
void MarkersManager::processClickCoordinates(int iDir, float x, float y){

  if(iDir<0 || iDir>=(int)fMarkersContainer.size() || fMarkersContainer.at(iDir)) return;  
  if(firstMarker){ x = firstMarker->GetX(); }

  int iMarkerColor = 2;
  int iMarkerStyle = 8;
  int iMarkerSize = 1;
  TMarker *aMarker = new TMarker(x, y, iMarkerStyle);
  aMarker->SetMarkerColor(iMarkerColor);
  aMarker->SetMarkerSize(iMarkerSize);
  aMarker->Draw();
  fMarkersContainer.at(iDir) = aMarker;
  if(!firstMarker){
    firstMarker = aMarker;
    drawFixedTimeLines(iDir, x);
  }
  else{
    clearLines();
    y = 30;//FIXME - calculate from other markers
    int missingMarkerDir = findMissingMarkerDir();
    aMarker = new TMarker(x, y, iMarkerStyle);
    aMarker->SetMarkerColor(iMarkerColor);
    aMarker->SetMarkerSize(iMarkerSize);
    std::string padName = "Histograms_"+std::to_string(missingMarkerDir+1);
    TPad *aPad = (TPad*)gROOT->FindObject(padName.c_str());
    aPad->cd();
    aMarker->Draw();
    fMarkersContainer[missingMarkerDir] = aMarker;    
  }
}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
void MarkersManager::drawFixedTimeLines(int iDir, double time){

  clearLines();
  int aColor = 1;
  TLine aLine(time, 0, time, 0);
  aLine.SetLineColor(aColor);
  aLine.SetLineWidth(2);
  
  for(int iDirTmp=0;iDirTmp<3;++iDirTmp){
    std::string padName = "Histograms_"+std::to_string(iDirTmp+1);
    TPad *aPad = (TPad*)gROOT->FindObject(padName.c_str());
    if(!aPad) continue;
    aPad->cd();
    TFrame *hFrame = (TFrame*)aPad->GetListOfPrimitives()->At(0);
    if(!hFrame) continue;
    double minY = hFrame->GetY1();
    double maxY = hFrame->GetY2();
    fLinesContainer[iDirTmp] = aLine.DrawLine(time, minY, time, maxY);
  }  
}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
int MarkersManager::findMissingMarkerDir(){

 for(int iDirTmp=0;iDirTmp<3;++iDirTmp){
    TMarker *item = fMarkersContainer.at(iDirTmp);
    if(!item) return iDirTmp;
    }
    return -1;
}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
double MarkersManager::getMissingYCoordinate(unsigned int missingMarkerDir){

  ///Not implemented yet.
  return 0.0;
  
}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
void MarkersManager::HandleMarkerPosition(Int_t event, Int_t x, Int_t y, TObject *sel){
  
  TObject *select = gPad->GetSelected();
  std::string objName = "";
  if(select) objName = std::string(select->GetName());
  if(event == kButton1 && objName.find("vs_time")!=std::string::npos){
    int iDir = (objName.find("U")!=std::string::npos) +
      2*(objName.find("V")!=std::string::npos) +
      3*(objName.find("W")!=std::string::npos) - 1;

    TVirtualPad *aCurrentPad = gPad->GetSelectedPad();
    aCurrentPad->cd();
    float localX = aCurrentPad->AbsPixeltoX(x);
    float localY = aCurrentPad->AbsPixeltoY(y);
    processClickCoordinates(iDir, localX, localY);		     
  }
  return;
}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
Bool_t MarkersManager::HandleButton(Int_t id){
   switch (id) {
   case M_ADD_SEGMENT:
    {
      std::cout<<KRED<<"Button!"<<RST<<std::endl;
    }
    break;
   }
   return kTRUE;
}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
void MarkersManager::DoButton(){
 TGButton* button = (TGButton*)gTQSender;
   UInt_t button_id = button->WidgetId();
   HandleButton(button_id);
 }
////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
