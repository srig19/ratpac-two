#ifndef __RAT_LightPathCalculator_hh__
#define __RAT_LightPathCalculator_hh__

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
  LightPathCalculator(const std::string& materialIV, const std::string& materialAcrylic, const std::string& materialOV);

  /// Constructor enforcing exact refractive indices rather than pulling from ratDB table
  LightPathCalculator(Double_t refIV, Double_t refAcrylic, Double_t refOV);

  /// Part of constructor that is the same for both constructors. Eos geometry related variables
  void SetValues();

  // Resets the variables that are calculated and filled during "CalcByPosition"
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
  TVector3 PathRefraction(const TVector3& incidentVec, const TVector3& incidentPoint, const Double_t incRIndex,
                          const Double_t refRIndex);

  /// Calculate the refracted path. If (fStraightLine) is true, this performs no refractions.
  void PathCalculation();

  /// Calculate the point where the light path, specified by the parameters, intersects either the inner or outer edge
  /// of the acrylic. This is robust, can be used for either IV->Inner Edge, Inner Edge->Outer Edge, or Outer Edge->PMT
  ///
  /// @param[in] initPos The initial point of the light path
  /// @param[in] initDir The direction of the light from initPos (Unit Vector)
  /// @param[in] outerEdge [True] indicates if we want the intersection point with the outer edge of the acrylic. Inner
  /// edge if False
  ///
  /// @return The point where the light path intersects either edge of the acrylic.
  TVector3 IntersectAcrylic(const TVector3& initPos, const TVector3& initDir, const bool& outerEdge);

  // Wrapper function to avoid TGraph using extrapolation beyond the bounding points.
  // It's possible that using the default extrapolation returns unphysical (negative) values
  /// @param[in] graph The graph containing xy datasets loaded from the db
  /// @param[in] energy The wavelength in nm
  ///
  /// @return The refractive index in the scintillator for this wavelength (energy)
  Double_t InterpolateTGraph(const TGraph& graph, const Double_t wl);

  // Function to print out the private member variables
  void PrintPrivateVariables();

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
  Double_t fUseOpticsTables;     // Flag that dictates if LPC is using the ratdb optics information

  std::string fIVMaterial;       // IV Material, initialized if the material constructor is used
  std::string fAcrylicMaterial;  // Acrylic Material, initialized if the material constructor is used
  std::string fOVMaterial;       // OV Material, initialized if the material constructor is used

  TVector3 fIncidentVecOnPMT;  // Final light path direction (unit normalised)
  TVector3 fInitialLightVec;   // Initial light path direction (unit normalised)

  TVector3 fStartPos;           // Start position of the light path (event position)
  TVector3 fEndPos;             // Required end position of the light path (PMT position)
  TVector3 fLightPathEndPos;    // Calculated end position of the light path
  std::string fStartingRegion;  // The region of the starting position.

  Double_t fPathPrecision;  // The accepted path proximity/tolerance to the PMT location [mm]
  Bool_t fIsTIR;            // TRUE: Total Internal Reflection encountered FALSE: It wasn't
  Bool_t fWithinError;      // TRUE: Difficult path to resolve and calculate FALSE: It wasn't
  Bool_t fStraightLine;     // TRUE: Light Path is a straight line approximation FALSE: It isn't

  /// Refraction points
  /// Depending on Light Path, could intersect the acrylic up to four times (counting inner and outer points separately)

  TVector3 fPointOnAcrylic1st;  // Point on acrylic where light path first hits the acrylic
  TVector3 fPointOnAcrylic2nd;  // Point on acrylic where light path hits the acrylic a second time
  TVector3 fPointOnAcrylic3rd;  // Point on acrylic where light path hits the acrylic a third time
  TVector3 fPointOnAcrylic4th;  // Point on acrylic where light path hits the acrylic a fourth time

  Double_t fDistInIV;       // Distance in the inner vessel region
  Double_t fDistInAcrylic;  // Distance in the acrylic
  Double_t fDistInOV;       // Distance in the OV

  Double_t fWavelength;  // The value of the wavelength in nm

  Double_t fSolidAngle;   // The solid angle subtended by the PMT for this light path
  Double_t fCosThetaAvg;  // Average incident angle on the PMT for this path, only calculated after a call to
                          // CalculateSolidAngle

  Double_t fFresnelTCoeff;  // The combined Fresnel TRANSMISSION coefficient for this path
  Double_t fFresnelRCoeff;  // The combined Fresnel REFLECTIVITY coefficient for this path
};
}  // namespace RAT

#endif
