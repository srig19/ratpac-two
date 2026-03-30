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
/// Light Path Types
/// 'IAO'   - Inner Vessel -> Acrylic -> Outer Vessel -> PMT
/// 'AO'    - Acrylic -> Outer Vessel -> PMT
/// 'OAIAO' - Outer Vessel -> Acrylic -> Inner Vessel -> Acrylic -> Outer Vessel -> PMT
/// 'AIAO'  - Acrylic -> Inner Vessel -> Acrylic -> Outer Vessel -> PMT
/// 'OAO'   - Outer Vessel -> Acrylic -> Outer Vessel -> PMT
/// 'O'     - Outer Vessel -> PMT
/// 'Null'  - Not Intitialized

enum eLightPathType { IAO, AO, OAIAO, AIAO, OAO, O, Null };

class LightPathCalculator {
 public:
  /////////////////////////////////
  ////////     METHODS     ////////
  /////////////////////////////////

  /// Default Constructor
  LightPathCalculator(std::string materialIV, std::string materialAcrylic, std::string materialOV);

  /// Constructor enforcing exact refractive indices rather than pulling from ratDB table
  LightPathCalculator(Double_t refIV, Double_t refAcrylic, Double_t refOV);

  /// Part of constructor that is the same for both constructors
  void SetValues();

  void ResetValues();

  /// Main analysis loop. This calculates the distances in each material and the angle at each
  /// intersection within the precision of 'precision' [mm] Default precision is 0.0 mm, which
  /// automatically invokes the straight line approximation (no refraction). This is numerically
  /// calculated via the RTSafe method, which is in a separate function below.
  void CalcByPosition(TVector3 startPos, TVector3 endPos, double wavelength, double precision = 0.0);

  /// Calculates the distances for light paths which start in the IV
  /// @param[in] eventPos The starting location of the light path (inside the AV)
  /// @param[out] pmtPos The end location of the light path (outside the AV)
  Bool_t CalculateDistancesIV(const TVector3& eventPos, const TVector3& pmtPos);

  /// Calculates the distances for light paths which start in the OV.
  /// @param[in] eventPos The starting location of the light path (outside of the AV)
  /// @param[out] pmtPos The end location of the light path (outside the AV)
  Bool_t CalculateDistancesOV(const TVector3& eventPos, const TVector3& pmtPos);

  /// Determine which region the point is in: IV, acrylic, OV, or past-PMT
  std::string DetermineRegion(const TVector3& point);

  /// Determine the vector normal to the acrylic at the provided point.
  /// The point should be within some specified tolerance to the acrylic boundaries.
  /// Tolerance [mm] for how close a point must be to the acrylic edges for this to work.
  /// NEEDS TO BE ITERATIVELY CHANGED
  TVector3 GetNormalVector(const TVector3& point, double tolerance = 3.0);

  /// Calculate refracted photon vector (unit normalised) using Snell's Law
  ///
  /// @param[in] incidentVec The unit vector incident on the surface surface
  /// @param[in] incientSurfVec The unit surface vector of the incident surface
  /// @param[in] incRIndex The incident refractive index
  /// @param[in] refRIndex The refractive index of the refractive media
  ///
  /// @return The refracted unit vector. Same as incidentVec if (fStraightLine) is true.
  TVector3 PathRefraction(const TVector3& incidentVec, const TVector3& incidentSurfVec, const Double_t incRIndex,
                          const Double_t refRIndex);

  /// Calculate the refracted path. This is called when a RTSafe fails or
  /// if a straight line is requested by the user from the beginning.
  /// If (fStraightLine) is true, this performs no refractions.
  /// @param[in] initDir The initial unit vector pointing from the starting position for the light path direction
  void PathCalculation(const TVector3& initDir);

  /// Calculate the point where the light path, specified by the parameters, intersects either the inner or outer edge
  /// of the acrylic. This is robust, can be used for either IV->Inner Edge, Inner Edge->Outer Edge, or Outer Edge->PMT
  ///
  /// @param[in] initPos The initial point of the light path
  /// @param[in] initDir The direction of the light from initPos (Unit Vector)
  /// @param[in] outerEdge [True] indicates if we want the intersection point with the outer edge of the acrylic. Inner
  /// edge if False
  ///
  /// @return The point where the light path intersects either edge of the acrylic.
  TVector3 LightPathCalculator::IntersectAcrylic(const TVector3& initPos, const TVector3& initDir,
                                                 const bool& outerEdge);

  //////////////////////////////////////////////////
  ////////          UTILITY ROUTINES        ////////
  //////////////////////////////////////////////////

  /// Utility Routines for refraction between IV/ Acrylic / OV
  /// 1-3: Inside and Outside light path types (IAO and O)
  /// 4-5: Outside light path types only (OAIAO)
  /// ThetaResidual: The difference between fPMTTargetTheta and the sum
  /// of ( Theta1st + Theta2nd + Theta3rd )
  /// or ( Theta1st + Theta2nd + Theta3rd + Theta4th + Theta5th )
  /// depending on light path type

  /// For a typical path there are various angles defined. A path is calculated
  /// in the 2D-plane containing BOTH the start(source) position and the end(pmt)
  /// position. This ensures that the path calculated is the minimum refracted
  /// path.
  /// To begin, a specific cordinate system is set up for each path, the (x,y,z) directions of this
  /// coordinate system are defined as follows:

  /// x : The direction defined by the radial vector pointing from the origin (centre of IV)
  ///     to the source position.
  /// z : The direction perpendicular to both the radial vector from the origin (centre of IV)
  ///     to the source position, and the radial vector from the origin (centre of IV) to the
  ///     end position. Mathematically speaking, the z-unit direction defines the 2D-plane
  ///     in which the path is calculated
  /// y : The direction defined by the cross product of 'z CROSS x'. The y direction
  ///     is therefore perpendicular to both x and z directions, but lies in the same
  ///     2D-plane as the x direction.

  /// These are utility routines which calculate the angle of refraction between
  /// the material boundaries in the detector. For a typical path whose start position is
  /// inside the IV, there will be three sets of angles:

  /// 1. The angle between the source position and the IV and Acrylic interface point
  ///    as viewed from the origin (centre of IV).
  /// 2. The angle between the first interface point (see above [1.]) and the Acrylic/OV interface
  ///    point as viewed from the origin (centre of IV).
  /// 3. The angle between the second interface point (see above [2.]) and the end position
  ///    of the path (PMT position).

  /// In the following, 'theta' is the test value for the initial direction of the path
  /// which is minimised against in RTSafe(). It is the initial and defining angle which
  /// ultimate determines the paths course throughout the rest of the detector and consequently
  /// the values of Theta1, Theta2, Theta3 etc. which follow as a result.

  /// The angle between the source position and the first AV intersection point
  /// as viewed from the origin (centre of AV)
  ///
  /// @param[in] theta (passed by reference).
  /// @return The angle between the source position and the first intersection point
  /// for this value of theta.
  Double_t Theta1st(const Double_t& theta);

  /// The derivative on this first angle (Theta1) with respect to 'theta'
  ///
  /// @param[in] theta (passed by reference).
  /// @return The derivative on this first angle (Theta1) with respect to 'theta'
  Double_t DTheta1st(const Double_t& theta);

  /// The angle between the first AV intersection point and the second AV interseciton point
  /// as viewed from the origin (centre of AV)
  ///
  /// @param[in] theta (passed by reference).
  /// @return The angle between the first AV intersection point and the second AV interseciton
  /// point as viewed from the origin (centre of AV)
  Double_t Theta2nd(const Double_t& theta);

  /// The derivative on this second angle (Theta2) with respect to 'theta'
  ///
  /// @param[in] theta (passed by reference).
  /// @return The derivative on this second angle (Theta2) with respect to 'theta'
  Double_t DTheta2nd(const Double_t& theta);

  /// The angle between the second AV intersection point and either the PMT position (source
  /// inside the AV) or the third intersection point (source outside the AV).
  ///
  /// @param[in] theta (passed by reference).
  /// @return The angle between the second AV intersection point and either the
  /// PMT position (source inside the AV) or the third intersection point
  /// (source outside the AV).
  Double_t Theta3rd(const Double_t& theta);

  /// The derivative on this third angle (Theta3) with respect to 'theta'
  ///
  /// @param[in] theta (passed by reference).
  /// @return The derivative on this third angle (Theta3) with respect to 'theta'
  Double_t DTheta3rd(const Double_t& theta);

  /// The angle between the third AV intersection point and the fourth
  /// intersection point (source outside the AV).
  ///
  /// @param[in] theta (passed by reference).
  /// @return The angle between the third AV intersection point and the fourth
  /// intersection point (source outside the AV).
  Double_t Theta4th(const Double_t& theta);

  /// The derivative on this fourth angle (Theta4) with respect to 'theta'
  ///
  /// @param[in] theta (passed by reference).
  /// @return The derivative on this fourth angle (Theta4) with respect to 'theta'
  Double_t DTheta4th(const Double_t& theta);

  /// The angle between the fourth AV intersection point and the fifth
  /// intersection point (source outside the AV).
  ///
  /// @param[in] theta (passed by reference).
  /// @return The angle between the third AV intersection point and the fifth
  /// intersection point (source outside the AV).
  Double_t Theta5th(const Double_t& theta);

  /// The derivative on this fifth angle (Theta5) with respect to 'theta'
  ///
  /// @param[in] theta (passed by reference).
  /// @return The derivative on this fifth angle (Theta5) with respect to 'theta'
  Double_t DTheta5th(const Double_t& theta);

  /// Calculate the residual between the target angle between the source and PMT position
  /// and the calculated value.
  ///
  /// @param[in] theta (passed by reference).
  /// @return The residual between the target angle between the source and PMT position
  /// and the calculated value.
  Double_t ThetaResidual(const Double_t& theta);

  /// Calculate the derivative on this residual with respect to 'theta'
  ///
  /// @param[in] theta (passed by reference).
  /// @return The derivative on this residual with respect to 'theta'
  Double_t DThetaResidual(const Double_t& theta);

  /// Using a combination of Newton-Raphson and bisection methods,
  /// this function returns the root of the function 'LightPathCalculator::Func'
  /// defined on the domain [x1, x2] to within an acceptable accuracy of +/- xAcc
  /// It returns the minimised (optimal) value of 'theta'
  ///
  /// @param[in] x1 The minimum possible value of 'theta' required.
  /// @param[in] x2 The maximum possible value of 'theta' required.
  /// @param[in] xAcc The acceptable precision allowed on this value of 'theta'
  /// @return The minimised, optimal value of 'theta' for this path.
  Double_t RTSafe(Double_t x1, Double_t x2, Double_t xAcc = 1.0e-3);

  /// Utility function used by 'RTSafe()' to perform the minimiation for
  /// the optimal value of 'theta'.
  ///
  /// @param[in] theta The test value of 'theta' for this path.
  /// @param[in,out] funcVal Computation of 'ThetaResidual()' based on 'theta'
  /// @param[in,out] dFuncVal The derivative on the above value with respect to 'theta'
  void FuncD(Double_t theta, Double_t& funcVal, Double_t& dFuncVal);

  /////////////////////////////////
  ////////     GETTERS     ////////
  /////////////////////////////////

  /// Light Path Type
  std::string GetLightPathType() { return fLightPathTypeMap[fLightPathType]; }

  /// Eos Geometry
  Double_t GetIVCylRadius() const { return fIVCylRadius; }
  Double_t GetIVCylHeight() const { return fIVCylHeight; }
  Double_t GetIVCapRadius() const { return fIVCapRadius; }
  Double_t GetAcrylicThickness() const { return fAcrylicThickness; }

  /// PMT Geometry
  Double_t GetBarrelPMTRadius() const { return fBarrelPMTRadius; }
  Double_t GetBarrelPMTHeight() const { return fBarrelPMTHeight; }
  Double_t GetTopPMTRadius() const { return fTopPMTRadius; }
  Double_t GetBotPMTRadius() const { return fBotPMTRadius; }

  /// Refractive Indices
  Double_t GetIVRefIndex() const { return fIVRefIndex; }
  Double_t GetAcrylicRefIndex() const { return fAcrylicRefIndex; }
  Double_t GetOVRefIndex() const { return fOVRefIndex; }

  /// Light Path Vectors
  TVector3 GetIncidentVecOnPMT() const { return fIncidentVecOnPMT; }
  TVector3 GetInitialLightVec() const { return fInitialLightVec; }

  /// Positions
  TVector3 GetStartPos() const { return fStartPos; }
  TVector3 GetEndPos() const { return fEndPos; }
  TVector3 GetLightPathEndPos() const { return fLightPathEndPos; }
  Double_t GetPMTTargetTheta() const { return fPMTTargetTheta; }

  /// Path Status Flags
  Bool_t GetIsTIR() const { return fIsTIR; }
  Bool_t GetResvHit() const { return fResvHit; }
  Bool_t GetStraightLine() const { return fStraightLine; }

  /// Refraction Points
  TVector3 GetPointOnAcrylic1st() const { return fPointOnAcrylic1st; }
  TVector3 GetPointOnAcrylic2nd() const { return fPointOnAcrylic2nd; }
  TVector3 GetPointOnAcrylic3rd() const { return fPointOnAcrylic3rd; }
  TVector3 GetPointOnAcrylic4th() const { return fPointOnAcrylic4th; }

  /// Distances
  Double_t GetDistInIV() const { return fDistInIV; }
  Double_t GetDistInAcrylic() const { return fDistInAcrylic; }
  Double_t GetDistInOV() const { return fDistInOV; }

  /// Energy
  Double_t GetEnergy() const { return fEnergy; }

  /// Solid Angle & Incident Angle
  Double_t GetSolidAngle() const { return fSolidAngle; }
  Double_t GetCosThetaAvg() const { return fCosThetaAvg; }

  /// Fresnel Coefficients
  Double_t GetFresnelTCoeff() const { return fFresnelTCoeff; }
  Double_t GetFresnelRCoeff() const { return fFresnelRCoeff; }

  /////////////////////////////////
  ////////     SETTERS     ////////
  /////////////////////////////////

  /// IV Geometry
  void SetIVCylRadius(const Double_t radius) { fIVCylRadius = radius; }
  void SetIVCylHeight(const Double_t height) { fIVCylHeight = height; }
  void SetIVCapRadius(const Double_t radius) { fIVCapRadius = radius; }
  void SetAcrylicThickness(const Double_t thickness) { fAcrylicThickness = thickness; }

  /// PMT Geometry
  void SetBarrelPMTRadius(const Double_t radius) { fBarrelPMTRadius = radius; }
  void SetBarrelPMTHeight(const Double_t height) { fBarrelPMTHeight = height; }
  void SetTopPMTRadius(const Double_t radius) { fTopPMTRadius = radius; }
  void SetBotPMTRadius(const Double_t radius) { fBotPMTRadius = radius; }

  /// Refractive Indices
  void SetIVRefIndex(const Double_t index) { fIVRefIndex = index; }
  void SetAcrylicRefIndex(const Double_t index) { fAcrylicRefIndex = index; }
  void SetOVRefIndex(const Double_t index) { fOVRefIndex = index; }

  /// Light Path Vectors
  void SetIncidentVecOnPMT(const TVector3& vec) { fIncidentVecOnPMT = vec; }
  void SetInitialLightVec(const TVector3& vec) { fInitialLightVec = vec; }

  /// Positions
  void SetStartPos(const TVector3& pos) { fStartPos = pos; }
  void SetEndPos(const TVector3& pos) { fEndPos = pos; }
  void SetLightPathEndPos(const TVector3& pos) { fLightPathEndPos = pos; }
  void SetPMTTargetTheta(const Double_t theta) { fPMTTargetTheta = theta; }

  /// Path Status Flags
  void SetIsTIR(const Bool_t isTIR) { fIsTIR = isTIR; }
  void SetResvHit(const Bool_t resvHit) { fResvHit = resvHit; }
  void SetStraightLine(const Bool_t straightLine) { fStraightLine = straightLine; }

  /// Refraction Points
  void SetPointOnAcrylic1st(const TVector3& point) { fPointOnAcrylic1st = point; }
  void SetPointOnAcrylic2nd(const TVector3& point) { fPointOnAcrylic2nd = point; }
  void SetPointOnAcrylic3rd(const TVector3& point) { fPointOnAcrylic3rd = point; }
  void SetPointOnAcrylic4th(const TVector3& point) { fPointOnAcrylic4th = point; }

  /// Distances
  void SetDistInIV(const Double_t dist) { fDistInIV = dist; }
  void SetDistInAcrylic(const Double_t dist) { fDistInAcrylic = dist; }
  void SetDistInOV(const Double_t dist) { fDistInOV = dist; }

  /// Energy
  void SetEnergy(const Double_t energy) { fEnergy = energy; }

  /// Solid Angle & Incident Angle
  void SetSolidAngle(const Double_t solidAngle) { fSolidAngle = solidAngle; }
  void SetCosThetaAvg(const Double_t cosThetaAvg) { fCosThetaAvg = cosThetaAvg; }

  /// Fresnel Coefficients
  void SetFresnelTCoeff(const Double_t coeff) { fFresnelTCoeff = coeff; }
  void SetFresnelRCoeff(const Double_t coeff) { fFresnelRCoeff = coeff; }

 private:
  eLightPathType fLightPathType;  // Light path type, based on what regions of the detector the path enters
  std::map<eLightPathType, std::string> fLightPathTypeMap;  // Map containing a descriptor for the light path type

  /// Private Variables describing the Eos geometry
  Double_t fIVCylRadius;       // Radius of the IV cylindrical region
  Double_t fIVCylHeight;       // Height of the IV cylindrical region
  Double_t fIVCapRadius;       // Radius of the spherical caps on the top/bottom
  Double_t fIVCapOffset;       // Offset between the center of the spherical caps and the origin of the IV
  Double_t fAcrylicThickness;  // Thickness of the acrylic

  /// Private Variables describing the PMT geometry
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
  Double_t fUseOpticsTables;

  TVector3 fIncidentVecOnPMT;  // Final light path direction (unit normalised)
  TVector3 fInitialLightVec;   // Initial light path direction (unit normalised)

  TVector3 fStartPos;           // Start position of the light path (event position)
  TVector3 fEndPos;             // Required end position of the light path (PMT position)
  TVector3 fLightPathEndPos;    // Calculated end position of the light path
  Double_t fPMTTargetTheta;     // The target PMT theta angle for the light path
  std::string fStartingRegion;  // The region of the starting position.

  Bool_t fIsTIR;         // TRUE: Total Internal Reflection encountered FALSE: It wasn't
  Bool_t fResvHit;       // TRUE: Difficult path to resolve and calculate FALSE: It wasn't
  Bool_t fStraightLine;  // TRUE: Light Path is a straight line approximation FALSE: It isn't

  /// Refraction points
  /// Depending on Light Path, could intersect the acrylic up to four times (counting inner and outer points separately)

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