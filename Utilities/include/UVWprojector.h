#ifndef __UVWPROJECTOR_H__
#define __UVWPROJECTOR_H__

// Class for 3D event projection to UVW strips.
// VERSION: 05 May 2018

#include <cstdlib>
#include <vector>
#include <map>

#include "root/include/TROOT.h"
#include "root/include/TH3D.h"
#include "root/include/TH2D.h"
#include "root/include/TH1.h"
#include "root/include/TH2Poly.h"

#include "GeometryTPC.h"

//class UVWprojector : public TH2Poly {
class UVWprojector {

 public:
  
  // empty constructor is required by TObject
  //  inline UVWprojector() { UVWprojector(nullptr,0,0,0); } 
  virtual ~UVWprojector();
  
  // Setter methods 
  
  UVWprojector(GeometryTPC *geo, int n=100, int nx=25, int ny=25);
  void SetAreaNpoints(int n); // change number of points used for random probing of TH2Poly bins
  void SetEvent3D(TH3D &h3); // 3D ionization map: (x [mm], y [mm], z [mm], Q [arb.u.])
  void SetEvent2D(TH2D &h2); // 2D ionization map: (x [mm], y [mm], Q [arb.u.])
  inline void SetDebug(bool flag) { _debug = flag; }

  // Getter methods
  
  TH1D    *GetStripProfile_TH1D(projection dir); // Get TH1D of time-integrated strip projection (SELECTED DIRECTION)
  TH2Poly *GetStripProfile_TH2Poly();     // Get TH2Poly of time-integrated strip projection (ALL STRIPS)
  TH2D    *GetStripVsTime_TH2D(projection dir);  // Get TH2D of strip vs time projection (SELECTED DIRECTION)  
  inline int GetAreaNpoints() { return area_npoints; }
  inline double GetEventIntegral() { return (input_hist==nullptr ? 0.0 : input_hist->Integral()); }
  
  // Nested helper class definitions
  struct BinFracMap {
    std::map<int /* TH1/TH2 bin's index */, double /* bin's fraction */ > FracMap; 
  };

 protected:

  // Setter methods

  bool InitAreaMapping();
  bool InitTimeMapping();
  virtual void AddBinContent(int32_t bin, double val);
  virtual void SetBinContent(int32_t bin, double val);
		  
  // Getter methods

  bool CheckBinsXY(TH3D *h1, TH3D *h2);
  bool CheckBinsXY(TH2D *h1, TH2D *h2);
  bool CheckBinsZ(TH3D *h1, TH3D *h2);
  
  int area_npoints;
  bool isOK_AreaMapping;
  bool isOK_TimeMapping;
  GeometryTPC *geo_ptr; // pointer to the existing TPC geometry
  TH1 *input_hist; // input histogram to be projected (can be TH3D/TH3F or TH2D/TH2F)
  bool is_input_2D; // is the event input histogram of TH2D or TH3D type?
  
  std::map<MultiKey2 /* TH2 bin index [1..NX*NY] */, BinFracMap, multikey2_less> fAreaFractionMap;
  std::map<int /* TH1 bin index [1..NX] */, BinFracMap> fTimeFractionMap; 
  
 private: 
  bool _debug;
  
  //  ClassDef(UVWprojector,1)
  
};

#endif
