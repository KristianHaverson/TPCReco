#include "RunConditions.h"
#include "colorText.h"

/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
void RunConditions::setDriftVelocity(double v) { vdrift = v;}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
void RunConditions::setSamplingRate(double r) { sampling_rate = r;}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
void RunConditions::setTriggerDelay(double d) { trigger_delay = d;}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
double RunConditions::getDriftVelocity() const { return vdrift;}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
double RunConditions::getSamplingRate() const { return sampling_rate;}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
double RunConditions::getTriggerDelay() const { return trigger_delay;}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
std::ostream & operator << (std::ostream &out, const RunConditions &aConditions){

  out<<KBLU<<"Run conditions: "<<RST<<std::endl;
  out<<KBLU<<" drift velocity [cm /us]: "<<RST<<aConditions.getDriftVelocity()<<std::endl;
  out<<KBLU<<" sampling rate: [MHz]   : "<<RST<<aConditions.getSamplingRate()<<std::endl;
  out<<KBLU<<" trigger delay: [MHz]   : "<<RST<<aConditions.getTriggerDelay()<<std::endl;

  return out;
}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
