#include <cstdlib>
#include <iostream>

#include <TGFrame.h>
#include <TGLabel.h>
#include <TGFontDialog.h>

#include <EntryDialog.h>
#include <MainFrame.h>
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
EntryDialog::EntryDialog(const TGWindow * p, MainFrame * aFrame)
 : TGCompositeFrame(p, 10, 10, kVerticalFrame), theMainFrame(aFrame){

   SetCleanup(kDeepCleanup);

   datasetInfoFrame = new TGHorizontalFrame(this, 300, 100);
   fileInfoFrame = new TGHorizontalFrame(this, 300, 100);
   modeInfoFrame = new TGHorizontalFrame(this, 300, 100);
   
   TGLayoutHints *aLayoutHints = new TGLayoutHints(kLHintsTop | kLHintsLeft  |
						 kLHintsShrinkX|kLHintsShrinkY |
						 kLHintsFillX|kLHintsFillY, 2, 2, 2, 2);

   AddFrame(modeInfoFrame, new TGLayoutHints(kLHintsFillX, 2, 2, 2, 2));
   AddFrame(fileInfoFrame, new TGLayoutHints(kLHintsFillX, 2, 2, 2, 2));
   AddFrame(datasetInfoFrame, aLayoutHints);
   
   TGGroupFrame *totalEventsFrame = new TGGroupFrame(datasetInfoFrame, "Events in the file:");
   TGGroupFrame *currentEventFrame = new TGGroupFrame(datasetInfoFrame, "Event id.:");
   TGGroupFrame *currentEntryFrame = new TGGroupFrame(datasetInfoFrame, "File entry:");
   TGGroupFrame *fileNameFrame = new TGGroupFrame(fileInfoFrame, "Processing file:");
   TGGroupFrame *modeFrame = new TGGroupFrame(modeInfoFrame, "Mode:");
   
   totalEventsLabel = new TGLabel(totalEventsFrame, "No input.");
   currentEventLabel = new TGLabel(currentEventFrame, "No input.");
   currentEntryLabel = new TGLabel(currentEntryFrame, "No input.");
   std::string tmp = "No input.";
   tmp.resize(fileNameLineLength,' ');
   fileNameLabel = new TGLabel(fileNameFrame,tmp.c_str());

   std::string mode = "NONE";
   modeLabel = new TGLabel(modeFrame, mode.c_str());
   ULong_t iColor;
   gClient->GetColorByName("red", iColor);
   modeLabel->SetTextColor(iColor); 

   aLayoutHints = new TGLayoutHints(kLHintsTop | kLHintsLeft, 5, 5, 5, 5);
   totalEventsFrame->AddFrame(totalEventsLabel, aLayoutHints);
   currentEventFrame->AddFrame(currentEventLabel, aLayoutHints);
   currentEntryFrame->AddFrame(currentEntryLabel, aLayoutHints);
   fileNameFrame->AddFrame(fileNameLabel, aLayoutHints);
   modeFrame->AddFrame(modeLabel, aLayoutHints);

   aLayoutHints = new TGLayoutHints(kLHintsTop | kLHintsLeft, 2, 2, 1, 1);
   datasetInfoFrame->AddFrame(totalEventsFrame, aLayoutHints);
   datasetInfoFrame->AddFrame(currentEventFrame, aLayoutHints);
   datasetInfoFrame->AddFrame(currentEntryFrame, aLayoutHints);
   fileInfoFrame->AddFrame(fileNameFrame, aLayoutHints);
   modeInfoFrame->AddFrame(modeFrame, aLayoutHints);
}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
EntryDialog::~EntryDialog(){

  delete datasetInfoFrame;
  delete totalEventsLabel;
  delete currentEventLabel;
  delete currentEntryLabel;

  delete fileInfoFrame;
  delete fileNameLabel;
}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
void EntryDialog::initialize(){ }
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
void EntryDialog::updateEventNumbers(unsigned int nTotalEvents,
                                     unsigned int iCurrentEvent,
				     unsigned int iCurrentEntry){

  totalEventsLabel->SetText(Form("%u", nTotalEvents));
  currentEventLabel->SetText(Form("%u",iCurrentEvent));
  currentEntryLabel->SetText(Form("%u",iCurrentEntry));  
  datasetInfoFrame->Layout();
}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
void EntryDialog::updateFileName(const std::string & fileName){

  std::string fileNameWithBreaks = fileName;
  size_t previousBreakPoint = 0;
  for(size_t iPos=0;iPos<fileNameWithBreaks.size();){
    iPos = fileNameWithBreaks.find("/",iPos+1);
    bool longPartFromStart = fileNameLineLength-iPos+previousBreakPoint<10;
    if(longPartFromStart){
      fileNameWithBreaks.insert(iPos+1,"\n");
      previousBreakPoint = iPos;
    }
  }

  size_t iPos = fileNameWithBreaks.find_last_of("/");
  bool longLineToEnd = fileNameWithBreaks.length() - iPos>fileNameLineLength;
  if(longLineToEnd){
    fileNameWithBreaks.insert(iPos+1,"\n");
  }
  
  fileNameLabel->SetText(fileNameWithBreaks.c_str());  
  fileInfoFrame->Layout();
}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
void EntryDialog::updateModeLabel(const std::string & aMode){

  modeLabel->SetText(aMode.c_str());
  modeInfoFrame->Layout();
}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
Bool_t EntryDialog::ProcessMessage(Long_t msg, Long_t parm1, Long_t /*parm2*/){
   switch (GET_MSG(msg)) {
   case kC_COMMAND:
      {
         switch (GET_SUBMSG(msg)) {
         case kCM_BUTTON:
            {
               switch (parm1) {
                  // exit button
               case 1:
                  {
                     std::cout<<"EntryDialog::ProcessMessage(): msg: "<<msg<<std::endl;
                     break;
                  }
                  // set button
               case 2:
                  {
                     break;
                  }
               }
               break;
            }
         }
         break;
      }
   }
   return kTRUE;
}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
