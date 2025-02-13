#ifndef __GEOMETRYTPC_H__
#define __GEOMETRYTPC_H__

// TPC geometry class.
// VERSION: 11 June 2022

#include <cstdlib>
#include <cstddef> 
#include <vector>
#include <map>
#include <memory>
#include <string>

#include "TROOT.h"
#include "TVector2.h"
#include "TVector3.h"
#include "TH2Poly.h"
#include "TGraph.h"
#include "TH3D.h"
#include "MultiKey.h"
#include "CommonDefinitions.h"
#include "RunConditions.h"
#include "GeometryStats.h"
#include "UtilsMath.h"
#include "StripTPC.h"

#define FPN_CH   3     // FPN channel type index
#define ERROR    -1    // error result indicator

class StripTPC;

// UVW strip geometry defined as a class

//class GeometryTPC : public TObject {
class GeometryTPC {

  friend class UVWprojector;

 private:
  mutable RunConditions runConditions;
  GeometryStats geometryStats;
  bool initOK;                          // was geometry initialized properly?
  int COBO_N;                           // total # of COBO boards in the system
  int AGET_Nchips;                      // # of AGET chips per ASAD board
  int AGET_Nchan;                       // # of channels in AGET chip (without FPN channels)
  int AGET_Nchan_fpn;                   // # of FPN channels in AGET chip
  int AGET_Nchan_raw;                   // # of total channels in AGET chip (including FPN channels)
  int AGET_Ntimecells;                  // # of time cells (buckets) in AGET chip
  std::map<int, int>         stripN;    // pair=(directory_idx, number_of_strips=xxx,xxx,xxx,4*ASAD_N*COBO_N)                          
  std::map<int, std::string> dir2name;  // pair=(directory_idx, group_name="U","V","W","FPN")
  std::map<std::string, int> name2dir;  // pair=(group_name="U","V","W","FPN", directory_idx)
  
  std::map<MultiKey4, std::shared_ptr<StripTPC> > mapByAget;     // key=(COBO_idx[0-1], ASAD_idx[0-3], AGET_idx[0-3], channel_idx[0-63])
  std::map<MultiKey4, std::shared_ptr<StripTPC> > mapByAget_raw; // key=(COBO_idx[0-1], ASAD_idx[0-3], AGET_idx [0-3], raw_channel_idx [0-67] )
  std::map<MultiKey3, std::shared_ptr<StripTPC> > mapByStrip;    // key=(STRIP_DIRECTION [0-2], STRIP_SECTION [0-2], STRIP_NUMBER [1-1024])
  std::map<int, int> ASAD_N;       // pair=(COBO_idx, number of ASAD boards)
  std::vector<int> FPN_chanId;     // FPN channels in AGET chips
  double pad_size;                 // in [mm]
  double pad_pitch;                // in [mm]
  double strip_pitch;              // in [mm]
  TVector2 reference_point;        // XY offset in [mm] of the REFERENCE POINT used to define the geometry
  TVector2 empty;                  //dummy vector returned for invalid input data
  std::map<int, TVector2> strip_unit_vec; // XY unit vector in [mm] for a given family of strips pointing towards ascending pad numbers of the strip
  std::map<int, TVector2> pitch_unit_vec; // XY unit vector in [mm] for a given family of strips pointing towards ascending strip numbers 
  std::map<int, std::shared_ptr<StripTPC> > fStripMap; // maps TH2Poly bin to a given StripTPC object
  double drift_zmin;               // lower drift cage acceptance limit along Z-axis [mm] (closest to readout PCB)
  double drift_zmax;               // upper drift cage acceptance limit along Z-axis [mm] (farthest from readout PCB)
  TH2Poly* tp;                     // for internal storage of arbitrary strip shapes
  int grid_nx;                     // partition size of TH2Poly in X-dir
  int grid_ny;                     // partition size of TH2Poly in Y-dir
  bool isOK_TH2Poly;               // is TH2Poly already initialized?
  bool _debug;                     // debug/verbose info flag
  TH2PolyBin* tp_convex;           // for internal storage of the convex hull for UVW active area

  // Setter methods 
  
  bool Load(const char *fname);                 // loads geometry from TXT config file
  bool LoadAnalog(std::istream &f);                            //subrutine. Loads analog channels from geometry TXT config file
  bool InitTH2Poly();                           // define bins for the underlying TH2Poly histogram

  void SetTH2PolyStrip(int ibin, std::shared_ptr<StripTPC> s);  // maps TH2Poly bin to a given StripTPC object

  bool InitActiveAreaConvexHull(TGraph *g);     // calculates convex hull from cloud of 2D points [mm]
  
 public:
  void Debug();
  GeometryTPC() { ; }  // empty constructor required by TObject
  //  virtual ~GeometryTPC();
  
  // Setter methods 
  GeometryTPC(const char* fname, bool debug=false);
  void SetTH2PolyPartition(int nx, int ny); // change cartesian binning of the underlying TH2Poly
  inline int GetTH2PolyPartitionX() const{ return grid_nx; }
  inline int GetTH2PolyPartitionY() const{ return grid_ny; }
  inline void SetDebug(bool flag) { _debug = flag; }

  void setDriftVelocity(double v);
  void setSamplingRate(double r);
  void setTriggerDelay(double d);
  
  // Getter methods
  const RunConditions & getRunConditions() const { return runConditions;}

  inline TH2Poly *GetTH2Poly() const{ return tp; }   // returns pointer to the underlying TH2Poly
  std::shared_ptr<StripTPC> GetTH2PolyStrip(int ibin)const;          // returns pointer to StripTPC object corresponding to TH2Poly bin 
  
  inline bool IsOK() const{ return initOK; }
  
  //returns total number of strips in sections
  int GetDirNstrips(int dir)const;
  int GetDirNstrips(std::string name)const;
  int GetDirNstrips(const char *name)const;
  int GetDirNstrips(std::shared_ptr<StripTPC> s)const;

  inline SectionIndexList GetDirSectionIndexList(int dir) const {return geometryStats.GetDirSectionIndexList(dir);}
  inline int GetDirNStrips(int dir,int section) const {return geometryStats.GetDirNStrips(dir,section);}
  inline int GetDirNSections(int dir) const {return GetDirSectionIndexList(dir).size();}
  inline int GetDirStripNSections(int dir, int num) const {return geometryStats.GetDirStripNSections(dir,num);}
  StripSectionBoundaryList GetStripSectionBoundaryList(int dir, int num) const {return geometryStats.GetStripSectionBoundaryList(dir,num);}
  inline int GetDirMinStrip(int dir,int section) const {return geometryStats.GetDirMinStrip(dir,section);}
  inline int GetDirMaxStrip(int dir,int section) const {return geometryStats.GetDirMaxStrip(dir,section);}
  inline int GetDirNStripsMerged(int dir) const {return geometryStats.GetDirNStripsMerged(dir);}
  inline int GetDirMinStripMerged(int dir) const {return geometryStats.GetDirMinStripMerged(dir);}
  inline int GetDirMaxStripMerged(int dir) const {return geometryStats.GetDirMaxStripMerged(dir);}

  inline int GetAgetNchips() const { return AGET_Nchips; }
  inline int GetAgetNchannels() const { return AGET_Nchan; }
  inline int GetAgetNchannels_raw() const { return AGET_Nchan_raw; }
  inline int GetAgetNchannels_fpn() const { return AGET_Nchan_fpn; }
  inline int GetAgetNtimecells() const { return AGET_Ntimecells; }
  inline int GetAsadNboards(int COBO_idx) const { return (ASAD_N.find(COBO_idx)==ASAD_N.end() ? 0 : ASAD_N.at(COBO_idx)); }
  inline int GetAsadNboards() const { int n=0; for(int icobo=0; icobo<COBO_N; icobo++) { n+=ASAD_N.at(icobo); } return n; }
  inline int GetCoboNboards() const { return COBO_N; }

  int GetDirIndex(const char *name)const; 
  int GetDirIndex(std::string name)const;
  int GetDirIndex(int global_channel_idx)const;
  int GetDirIndex(int COBO_idx, int ASAD_idx, int AGET_idx, int channel_idx)const;
  int GetDirIndex_raw(int global_raw_channel_idx)const;
  int GetDirIndex_raw(int COBO_idx, int ASAD_idx, int AGET_idx, int raw_channel_idx)const;

  const char* GetDirName(int dir)const;
  const char* GetStripName(std::shared_ptr<StripTPC> s)const;

  std::shared_ptr<StripTPC> GetStripByAget(int COBO_idx, int ASAD_idx, int AGET_idx, int channel_idx) const;         // valid range [0-1][0-3][0-3][0-63]
  std::shared_ptr<StripTPC> GetStripByGlobal(int global_channel_idx) const;                                          // valid range [0-1023]
  std::shared_ptr<StripTPC> GetStripByAget_raw(int COBO_idx, int ASAD_idx, int AGET_idx, int raw_channel_idx) const; // valid range [0-1][0-3][0-3][0-67]
  std::shared_ptr<StripTPC> GetStripByGlobal_raw(int global_raw_channel_idx) const;                                  // valid range [0-(1023+4*ASAD_N*COBO_N)]
  std::shared_ptr<StripTPC> GetStripByDir(int dir, int section, int num) const;                                      // valid range [0-2][0-2][1-1024]

  // various helper functions for calculating local/global normal/raw channel index
  int Aget_normal2raw(int channel_idx)const;                      // valid range [0-63]
  int Aget_raw2normal(int raw_channel_idx)const;                  // valid range [0-67]
  int Aget_fpn2raw(int FPN_idx)const;                             // valid range [0-3]

  int Asad_normal2raw(int aget_idx, int channel_idx)const;        // valid range [0-3][0-63]
  int Asad_normal2raw(int asad_channel_idx)const;                 // valid range [0-255]
  int Asad_normal2normal(int aget_idx, int channel_idx)const;     // valid range [0-3][0-63]
  int Asad_raw2normal(int aget_idx, int raw_channel_idx)const;    // valid range [0-3][0-67]
  int Asad_raw2normal(int asad_raw_channel_idx)const;             // valid range [0-271]
  int Asad_raw2raw(int aget_idx, int raw_channel_idx)const;       // valid range [0-3][0-67]

  int Global_normal2raw(int COBO_idx, int ASAD_idx, int aget_idx, int channel_idx)const;      // valid range [0-1][0-3][0-3][0-63]
  int Global_normal2raw(int glb_channel_idx)const;                                            // valid range [0-1024]
  int Global_normal2normal(int COBO_idx, int ASAD_idx, int aget_idx, int channel_idx)const;   // valid range [0-1][0-3][0-3][0-63]

  int Global_raw2normal(int COBO_idx, int ASAD_idx, int aget_idx, int raw_channel_idx)const;  // valid range [0-1][0-3][0-3][0-67]
  int Global_raw2normal(int glb_raw_channel_idx)const;                                        // valid range [0-(1023+4*ASAD_N*COBO_N)]
  int Global_raw2raw(int COBO_idx, int ASAD_idx, int aget_idx, int raw_channel_idx)const;     // valid range [0-1][0-3][0-3][0-67]
  int Global_fpn2raw(int COBO_idx, int ASAD_idx, int aget_idx, int FPN_idx)const;             // valid range [0-1][0-3][0-3][0-3]

  int Global_strip2normal(StripTPC *s) const;
  int Global_strip2normal(int dir, int section, int num)const;                                // valid range [0-2][0-2][1-1024]
  int Global_strip2raw(StripTPC *s)const;
  int Global_strip2raw(int dir, int section, int num)const;                                   // valid range [0-2][0-2][1-1024]

  bool GetCrossPoint(std::shared_ptr<StripTPC> strip1,
		     std::shared_ptr<StripTPC> strip2,
		     TVector2 &point)const;
  
  bool GetCrossPoint(StripTPC *strip1,
		     StripTPC *strip2,
		     TVector2 &point)const;
  
  bool GetUVWCrossPointInMM(int dir1, double UVW_pos1, int dir2, double UVW_pos2, TVector2 &point)const;
  bool MatchCrossPoint(StripTPC *strip1, StripTPC *strip2, StripTPC *strip3, double radius, TVector2 &point)const;
  bool MatchCrossPoint(std::shared_ptr<StripTPC> strip1,
		       std::shared_ptr<StripTPC> strip2,
		       std::shared_ptr<StripTPC> strip3,
		       double radius, TVector2 &point)const;

  inline double GetPadSize() const{ return pad_size; } // [mm]
  inline double GetPadPitch() const{ return pad_pitch; } // [mm]
  inline double GetStripPitch() const{ return strip_pitch; } // [mm]

  inline double GetDriftVelocity() const{ return runConditions.getDriftVelocity(); } // [cm/us]
  inline double GetSamplingRate() const{ return runConditions.getSamplingRate(); } // [MHz]
  inline double GetTriggerDelay() const{ return runConditions.getTriggerDelay(); } // [us]
  inline double GetTimeBinWidth() const{ return  1.0/GetSamplingRate()*GetDriftVelocity()*10.0; } // [mm]

  inline TVector2 GetReferencePoint() const{ return reference_point; } // XY ([mm],[mm])
  TVector2 GetStripUnitVector(int dir) const; // XY ([mm],[mm])
  TVector2 GetStripPitchVector(int dir) const; // XY ([mm],[mm])
  TVector3 GetStripPitchVector3D(int dir) const; // XYZ ([mm],[mm],0)  
  inline double GetDriftCageZmin() const{ return drift_zmin; } // [mm]
  inline double GetDriftCageZmax() const{ return drift_zmax; } // [mm]

  double Strip2posUVW(int dir, int section, int number, bool &err_flag)const; // [mm] (signed) distance of projection of (X=0, Y=0) point from projection of the central line of the (existing) strip on the strip pitch axis for a given direction
  double Strip2posUVW(int dir, int number, bool &err_flag)const; // [mm] (signed) distance of projection of (X=0, Y=0) point from projection of the central line of the (existing) strip on the strip pitch axis for a given direction
  double Strip2posUVW(std::shared_ptr<StripTPC> strip, bool &err_flag)const; // [mm] (signed) distance of projection of (X=0, Y=0) point from projection of the central line of the (existing) strip on the strip pitch axis for a strip given direction
  bool IsStripDirReversed(int dir) const;
  
  double Cartesian2posUVW(double x, double y, int dir, bool &err_flag)const; // [mm] (signed) distance of projection of (X=0, Y=0) point from projection of a given (X,Y) point on the strip pitch axis for a given direction
  double Cartesian2posUVW(TVector2 pos, int dir, bool &err_flag)const; // [mm] (signed) distance of projection of (X=0, Y=0) point from projection of a given (X,Y) point on the strip pitch axis for a given direction
  
  double Timecell2pos(double position_in_cells, bool &err_flag)const; // [mm] output: position along Z-axis
  double Pos2timecell(double z, bool &err_flag)const; // output: time-cell number, valid range [0-511]

  TGraph GetActiveAreaConvexHull(double vetoBand=0); // get convex hull [mm] of the entire UVW active area
                                                     // with (optionally) excluded outer VETO band [mm]
  bool IsInsideActiveVolume(TVector3 point); // checks if 3D point [mm] has X,Y inside
                                             // UVW active area and Z within [zmin, zmax] range
  bool IsInsideActiveArea(TVector2 point); // checks if 2D point [mm] is inside UVW active area
  bool IsInsideElectronicsRange(double z); // checks if Z coordinate [mm] is inside Z-slice covered by the GET electronics
  bool IsInsideElectronicsRange(TVector3 point); // checks 3D point [mm] is inside Z-slice covered by the GET electronics

  std::tuple<double, double> rangeX() const; //min/max X [mm] cartesian coordinates covered by UVW active area
  std::tuple<double, double> rangeY() const; //min/max Y [mm] cartesian coordinates covered by UVW active area
  std::tuple<double, double> rangeZ() const; //min/max Z [mm] cartesian coordinates covered by drift cage
  std::tuple<double, double, double, double> rangeXY() const; //min/max X and Y [mm] cartesian coordinates covered by UVW active area
  std::tuple<double, double, double, double, double, double> rangeXYZ() const; //min/max X and Y [mm] cartesian coordinates covered by UVW active area
                                                         // and min/max Z [mm] cartesian coordinates covered by drift cage
  std::tuple<double, double> rangeStripSectionInMM(int dir, int section) const; // [mm] min/max (signed) distance between projection of outermost strip's central axis and projection of the origin (X=0,Y=0) point on the U/V/W pitch axis for a given direction (per section)
  std::tuple<double, double> rangeStripDirInMM(int dir) const; // [mm] min/max (signed) distance between projection of outermost strip's central axis and projection of the origin (X=0,Y=0) point on the U/V/W pitch axis for a given direction (all sections)

  TH3D *Get3DFrame(int rebin_space, int rebin_time) const; //frame for plotting 3D reconstruction
  
  TH2D *GetXY_TestUV(TH2D *h=NULL); // auxillary functions for x-check 
  TH2D *GetXY_TestVW(TH2D *h=NULL); // auxillary functions for x-check 
  TH2D *GetXY_TestWU(TH2D *h=NULL); // auxillary functions for x-check
  
  bool operator==(const GeometryTPC&) const;
  inline bool operator!=(const GeometryTPC& B) const { return !(*this==B); }
  static const int outside_section{GeometryStats::outside_section}; // index for dummy sections outside of UVW active area used to truncate toy MC generated signals
  //  ClassDef(GeometryTPC,1)
};
#define __GEOMETRYTPC_H__
#endif
