#ifndef __LightPathCalculator_hh__
#define __LightPathCalculator_hh__

#include <dirent.h>
#include <stdio.h>

#include <algorithm>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

#include "G4PhysicalConstants.hh"
#include "RAT/DB.hh"
#include "TGraph.h"
#include "TMath.h"
#include "TObject.h"
#include "TVector3.h"

namespace RAT {
// Light Path Types
// 'IAO'   - Inner Vessel -> Acrylic -> Outer Vessel -> PMT
// 'AO'    - Acrylic -> Outer Vessel -> PMT
// 'OAIAO' - Outer Vessel -> Acrylic -> Inner Vessel -> Acrylic -> Outer Vessel -> PMT
// 'AIAO'  - Acrylic -> Inner Vessel -> Acrylic -> Outer Vessel -> PMT
// 'OAO'   - Outer Vessel -> Acrylic -> Outer Vessel -> PMT
// 'O'     - Outer Vessel -> PMT
// 'Null'  - Not Intitialized

enum eLightPathType { IAO, AO, OAIAO, AIAO, OAO, O, Null };

class LightPathCalculator {
 public:
  /////////////////////////////////
  ////////     METHODS     ////////
  /////////////////////////////////

  // Default Constructor
  LightPathCalculator(std::string materialIV, std::string materialAcrylic, std::string materialOV);

  // Constructor enforcing exact refractive indices rather than pulling from ratDB table
  LightPathCalculator(Double_t refIV, Double_t refAcrylic, Double_t refOV);

  // Part of constructor that is the same for both constructors
  void SetValues();

  void ResetValues();

  // Main analysis loop. This calculates the distances in each material and the angle at each
  // intersection within the precision of 'precision' [mm] Default precision is 0.0 mm, which
  // automatically invokes the straight line approximation (no refraction). This is numerically
  // calculated via the RTSafe method, which is in a separate function below.
  void CalcByPosition(TVector3 startPos, TVector3 endPos, double wavelength, double precision = 0.0);

  // Calculates the distances for light paths which start in the IV
  // @param[in] eventPos The starting location of the light path (inside the AV)
  // @param[out] pmtPos The end location of the light path (outside the AV)
  Bool_t CalculateDistancesIV(const TVector3& eventPos, const TVector3& pmtPos);

  // Calculates the distances for light paths which start in the OV.
  // @param[in] eventPos The starting location of the light path (outside of the AV)
  // @param[out] pmtPos The end location of the light path (outside the AV)
  Bool_t CalculateDistancesOV(const TVector3& eventPos, const TVector3& pmtPos);

  // Determine which region the point is in: IV, acrylic, OV, or past-PMT
  std::string GetRegion(const TVector3& point);

  /////////////////////////////////
  ////////     GETTERS     ////////
  /////////////////////////////////

  // Light Path Type
  std::string GetLightPathType() { return fLightPathTypeMap[fLightPathType]; }

  // Eos Geometry
  Double_t GetIVCylRadius() const { return fIVCylRadius; }
  Double_t GetIVCylHeight() const { return fIVCylHeight; }
  Double_t GetIVCapRadius() const { return fIVCapRadius; }
  Double_t GetAcrylicThickness() const { return fAcrylicThickness; }

  // PMT Geometry
  Double_t GetBarrelPMTRadius() const { return fBarrelPMTRadius; }
  Double_t GetBarrelPMTHeight() const { return fBarrelPMTHeight; }
  Double_t GetTopPMTRadius() const { return fTopPMTRadius; }
  Double_t GetBotPMTRadius() const { return fBotPMTRadius; }

  // Refractive Indices
  Double_t GetIVRefIndex() const { return fIVRefIndex; }
  Double_t GetAcrylicRefIndex() const { return fAcrylicRefIndex; }
  Double_t GetOVRefIndex() const { return fOVRefIndex; }

  // Light Path Vectors
  TVector3 GetIncidentVecOnPMT() const { return fIncidentVecOnPMT; }
  TVector3 GetInitialLightVec() const { return fInitialLightVec; }

  // Positions
  TVector3 GetStartPos() const { return fStartPos; }
  TVector3 GetEndPos() const { return fEndPos; }
  TVector3 GetLightPathEndPos() const { return fLightPathEndPos; }
  Double_t GetPMTTargetTheta() const { return fPMTTargetTheta; }

  // Path Status Flags
  Bool_t GetIsTIR() const { return fIsTIR; }
  Bool_t GetResvHit() const { return fResvHit; }
  Bool_t GetStraightLine() const { return fStraightLine; }

  // Refraction Points
  TVector3 GetPointOnAcrylic1st() const { return fPointOnAcrylic1st; }
  TVector3 GetPointOnAcrylic2nd() const { return fPointOnAcrylic2nd; }
  TVector3 GetPointOnAcrylic3rd() const { return fPointOnAcrylic3rd; }
  TVector3 GetPointOnAcrylic4th() const { return fPointOnAcrylic4th; }

  // Distances
  Double_t GetDistInIV() const { return fDistInIV; }
  Double_t GetDistInAcrylic() const { return fDistInAcrylic; }
  Double_t GetDistInOV() const { return fDistInOV; }

  // Energy
  Double_t GetEnergy() const { return fEnergy; }

  // Solid Angle & Incident Angle
  Double_t GetSolidAngle() const { return fSolidAngle; }
  Double_t GetCosThetaAvg() const { return fCosThetaAvg; }

  // Fresnel Coefficients
  Double_t GetFresnelTCoeff() const { return fFresnelTCoeff; }
  Double_t GetFresnelRCoeff() const { return fFresnelRCoeff; }

  /////////////////////////////////
  ////////     SETTERS     ////////
  /////////////////////////////////

  // IV Geometry
  void SetIVCylRadius(const Double_t radius) { fIVCylRadius = radius; }
  void SetIVCylHeight(const Double_t height) { fIVCylHeight = height; }
  void SetIVCapRadius(const Double_t radius) { fIVCapRadius = radius; }
  void SetAcrylicThickness(const Double_t thickness) { fAcrylicThickness = thickness; }

  // PMT Geometry
  void SetBarrelPMTRadius(const Double_t radius) { fBarrelPMTRadius = radius; }
  void SetBarrelPMTHeight(const Double_t height) { fBarrelPMTHeight = height; }
  void SetTopPMTRadius(const Double_t radius) { fTopPMTRadius = radius; }
  void SetBotPMTRadius(const Double_t radius) { fBotPMTRadius = radius; }

  // Refractive Indices
  void SetIVRefIndex(const Double_t index) { fIVRefIndex = index; }
  void SetAcrylicRefIndex(const Double_t index) { fAcrylicRefIndex = index; }
  void SetOVRefIndex(const Double_t index) { fOVRefIndex = index; }

  // Light Path Vectors
  void SetIncidentVecOnPMT(const TVector3& vec) { fIncidentVecOnPMT = vec; }
  void SetInitialLightVec(const TVector3& vec) { fInitialLightVec = vec; }

  // Positions
  void SetStartPos(const TVector3& pos) { fStartPos = pos; }
  void SetEndPos(const TVector3& pos) { fEndPos = pos; }
  void SetLightPathEndPos(const TVector3& pos) { fLightPathEndPos = pos; }
  void SetPMTTargetTheta(const Double_t theta) { fPMTTargetTheta = theta; }

  // Path Status Flags
  void SetIsTIR(const Bool_t isTIR) { fIsTIR = isTIR; }
  void SetResvHit(const Bool_t resvHit) { fResvHit = resvHit; }
  void SetStraightLine(const Bool_t straightLine) { fStraightLine = straightLine; }

  // Refraction Points
  void SetPointOnAcrylic1st(const TVector3& point) { fPointOnAcrylic1st = point; }
  void SetPointOnAcrylic2nd(const TVector3& point) { fPointOnAcrylic2nd = point; }
  void SetPointOnAcrylic3rd(const TVector3& point) { fPointOnAcrylic3rd = point; }
  void SetPointOnAcrylic4th(const TVector3& point) { fPointOnAcrylic4th = point; }

  // Distances
  void SetDistInIV(const Double_t dist) { fDistInIV = dist; }
  void SetDistInAcrylic(const Double_t dist) { fDistInAcrylic = dist; }
  void SetDistInOV(const Double_t dist) { fDistInOV = dist; }

  // Energy
  void SetEnergy(const Double_t energy) { fEnergy = energy; }

  // Solid Angle & Incident Angle
  void SetSolidAngle(const Double_t solidAngle) { fSolidAngle = solidAngle; }
  void SetCosThetaAvg(const Double_t cosThetaAvg) { fCosThetaAvg = cosThetaAvg; }

  // Fresnel Coefficients
  void SetFresnelTCoeff(const Double_t coeff) { fFresnelTCoeff = coeff; }
  void SetFresnelRCoeff(const Double_t coeff) { fFresnelRCoeff = coeff; }

 private:
  eLightPathType fLightPathType;  // Light path type, based on what regions of the detector the path enters
  std::map<eLightPathType, std::string> fLightPathTypeMap;  // Map containing a descriptor for the light path type

  // Private Variables describing the Eos geometry
  Double_t fIVCylRadius;       // Radius of the IV cylindrical region
  Double_t fIVCylHeight;       // Height of the IV cylindrical region
  Double_t fIVCapRadius;       // Radius of the spherical caps on the top/bottom
  Double_t fIVCapHeight;       //
  Double_t fAcrylicThickness;  // Thickness of the AV

  // Private Variables describing the PMT geometry
  Double_t fBarrelPMTRadius;  // Radius of the barrel PMT region
  Double_t fBarrelPMTHeight;  // Height of the barrel PMT region
  Double_t fTopPMTRadius;     // Radius of the top PMT array, measured from top of cylinder
  Double_t fBotPMTRadius;     // Radius of the bottom PMT array, measured from bottom of cylinder

  Double_t fIVRefIndex;          // The value of the IV refractive index used for this path
  TGraph fIVRefIndexGraph;       // IV refractive index vector from optics table
  Double_t fAcrylicRefIndex;     // The value of the acrylic vessel refractive index used for this path
  TGraph fAcrylicRefIndexGraph;  // Acrylic refractive index vector from optics table
  Double_t fOVRefIndex;          // The value of the OV refractive index used for this path
  TGraph fOVRefIndexGraph;       // OV refractive index vector from optics table

  TVector3 fIncidentVecOnPMT;  // Final light path direction (unit normalised)
  TVector3 fInitialLightVec;   // Initial light path direction (unit normalised)

  TVector3 fStartPos;         // Start position of the light path (event position)
  TVector3 fEndPos;           // Required end position of the light path (PMT position)
  TVector3 fLightPathEndPos;  // Calculated end position of the light path
  Double_t fPMTTargetTheta;   // The target PMT theta angle for the light path

  Bool_t fIsTIR;         // TRUE: Total Internal Reflection encountered FALSE: It wasn't
  Bool_t fResvHit;       // TRUE: Difficult path to resolve and calculate FALSE: It wasn't
  Bool_t fStraightLine;  // TRUE: Light Path is a straight line approximation FALSE: It isn't

  // Refraction points
  // Depending on Light Path, could intersect the acrylic up to four times (counting inner and outer points separately)

  TVector3 fPointOnAcrylic1st;  // Point on acrylic where light path first hits the acrylic
  TVector3 fPointOnAcrylic2nd;  // Point on acrylic where light path hits the acrylic a second time
  TVector3 fPointOnAcrylic3rd;  // Point on acrylic where light path hits the acrylic a third time
  TVector3 fPointOnAcrylic4th;  // Point on acrylic where light path hits the acrylic a fourth time

  Double_t fDistInIV;       // Distance in the inner vessel region
  Double_t fDistInAcrylic;  // Distance in the acrylic
  Double_t fDistInOV;       // Distance in the OV

  Double_t fEnergy;  // The value of the wavelength in MeV

  Double_t fSolidAngle;   // The solid angle subtended by the PMT for this light path
  Double_t fCosThetaAvg;  // Average incident angle on the PMT for this path, only calculated after a call to
                          // CalculateSolidAngle

  Double_t fFresnelTCoeff;  // The combined Fresnel TRANSMISSION coefficient for this path
  Double_t fFresnelRCoeff;  // The combined Fresnel REFLECTIVITY coefficient for this path

  Double_t fLoopCeiling;    // Iteration Ceiling for algortithm loop
  Double_t fFinalLoopSize;  // Final loop value which meets locality conditions
  Double_t fPathPrecision;  // The accepted path proximity/tolerance to the PMT location [mm]
};
}  // namespace RAT

#endif