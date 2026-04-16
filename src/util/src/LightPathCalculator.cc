//#include "/home/srikar-temp/ratpac-two/src/util/include/RAT/LightPathCalculator.hh"
#include "RAT/LightPathCalculator.hh"

namespace RAT {
// Constructor where user assigns the material for each region. Refractive indices are pulled from optics tables
LightPathCalculator::LightPathCalculator(const std::string& materialIV, const std::string& materialAcrylic,
                                         const std::string& materialOV) {
  std::cout << "Entering constructor" << std::endl;
  DB* db = DB::Get();
  fUseOpticsTables = true;
  fIVMaterial = materialIV;
  fAcrylicMaterial = materialAcrylic;
  fOVMaterial = materialOV;

  DBLinkPtr opticsIV;
  std::string optionIV;
  std::vector<double> valOneIV;
  std::vector<double> valTwoIV;

  DBLinkPtr opticsAcrylic;
  std::string optionAcrylic;
  std::vector<double> valOneAcrylic;
  std::vector<double> valTwoAcrylic;

  DBLinkPtr opticsOV;
  std::string optionOV;
  std::vector<double> valOneOV;
  std::vector<double> valTwoOV;
  std::cout << "After intializations" << std::endl;

  try {
    opticsIV = db->GetLink("OPTICS", fIVMaterial);
    optionIV = opticsIV->GetS("RINDEX_option");
    valOneIV = opticsIV->GetDArray("RINDEX_value1");
    valTwoIV = opticsIV->GetDArray("RINDEX_value2");
  } catch (DBNotFoundError& e) {
    std::cout << "LightPathCalculator: Failed to retrieve optics information for IV: " << materialIV << std::endl;
    return;
  }
  try {
    opticsAcrylic = db->GetLink("OPTICS", fAcrylicMaterial);
    optionAcrylic = opticsAcrylic->GetS("RINDEX_option");
    valOneAcrylic = opticsAcrylic->GetDArray("RINDEX_value1");
    valTwoAcrylic = opticsAcrylic->GetDArray("RINDEX_value2");
  } catch (DBNotFoundError& e) {
    std::cout << "LightPathCalculator: Failed to retrieve optics information for Acrylic: " << materialAcrylic
              << std::endl;
    return;
  }
  try {
    opticsOV = db->GetLink("OPTICS", fOVMaterial);
    optionOV = opticsOV->GetS("RINDEX_option");
    valOneOV = opticsOV->GetDArray("RINDEX_value1");
    valTwoOV = opticsOV->GetDArray("RINDEX_value2");
  } catch (DBNotFoundError& e) {
    std::cout << "LightPathCalculator: Failed to retrieve optics information for OV: " << materialOV << std::endl;
    return;
  }

  std::cout << "After linking" << std::endl;

  const double hbarc = 197.32705e-6;  // conversion from energy in MeV to nm
  for (int i = 0; i < 3; i++) {
    std::string option;
    std::vector<double> valOne;
    std::vector<double> valTwo;
    if (i == 0) {
      option = optionIV;
      valOne = valOneIV;
      valTwo = valTwoIV;
    } else if (i == 1) {
      option = optionAcrylic;
      valOne = valOneAcrylic;
      valTwo = valTwoAcrylic;
    } else if (i == 2) {
      option = optionOV;
      valOne = valOneOV;
      valTwo = valTwoOV;
    }

    if (valOne.size() != valTwo.size()) {
      std::cout << "LightPathCalculator: Optics table vectors are not the same size." << std::endl;
      return;
    }

    int pointNum = 0;
    for (size_t j = 0; j < valOne.size(); j++) {
      if (option == "wavelength") {
        if (i == 0) {
          fIVRefIndexGraph.SetPoint(pointNum, valOne[j], valTwo[j]);
        } else if (i == 1) {
          fAcrylicRefIndexGraph.SetPoint(pointNum, valOne[j], valTwo[j]);
        } else if (i == 2) {
          fOVRefIndexGraph.SetPoint(pointNum, valOne[j], valTwo[j]);
        }
      } else if (option == "energy") {
        if (i == 0) {
          fIVRefIndexGraph.SetPoint(pointNum, (2 * TMath::Pi() * hbarc) / valOne[j], valTwo[j]);
        } else if (i == 1) {
          fAcrylicRefIndexGraph.SetPoint(pointNum, (2 * TMath::Pi() * hbarc) / valOne[j], valTwo[j]);
        } else if (i == 2) {
          fOVRefIndexGraph.SetPoint(pointNum, (2 * TMath::Pi() * hbarc) / valOne[j], valTwo[j]);
        }
      }
      pointNum++;
    }
  }

  SetValues();
}

// Constructor where user chooses a single refractive index for each material
LightPathCalculator::LightPathCalculator(const Double_t& refIV, const Double_t& refAcrylic, const Double_t& refOV) {
  std::cout << "Entering constructor" << std::endl;
  fUseOpticsTables = false;

  fIVRefIndex = refIV;
  fAcrylicRefIndex = refAcrylic;
  fOVRefIndex = refOV;

  SetValues();
}

// Common aspects of constructor between both constructor options above
void LightPathCalculator::SetValues() {
  std::cout << "Setting values" << std::endl;
  /*DB* db = DB::Get();
  std::cout << "Got DB" << std::endl;
  db->LoadDefaults();
  std::cout << "Loaded defaults" << std::endl;
  db->Load("~/EosSimulations/ratdb/Eos/Eos.geo");
  std::cout << "Loaded geometry file" << std::endl;
  //db->Load (db->GetLink("DETECTOR")->GetS("geo_file"));
  //  Setting geometry specific values
  DBLinkPtr geoPtr = DB::Get()->GetLink("GEO", "eos_inner");
  fIVCylRadius = geoPtr->GetD("r_max");
  fIVCylHeight = geoPtr->GetD("size_z");  // This is HALF the height
  fIVCapRadius = geoPtr->GetD("top_radius");*/

  // For now, the ratdb links will be simplified because this isn't working. COME BACK TO IT

  fIVCylRadius = 888.8;
  fIVCylHeight = 676.91;
  fIVCapRadius = 1699.51;

  fIVCapOffset = TMath::Abs(fIVCapRadius - fIVCylHeight);
  fAcrylicThickness = 25.4;  // Acrylic is 1 inch thick
  fBarrelPMTRadius = 1049.6;
  fBarrelPMTHeight = 658.8;  // This is HALF the height
  fTopPMTRadius = 1070.0;
  fBotPMTRadius = 1070.0;

  // Making the Light Path Type Map
  fLightPathTypeMap[IAO] = "IV->Acrylic->OV";
  fLightPathTypeMap[AIAO] = "Acrylic->IV->Acrylic->OV";
  fLightPathTypeMap[AO] = "Acrylic->OV";
  fLightPathTypeMap[O] = "OV";
  fLightPathTypeMap[OAO] = "OV->Acrylic->OV";
  fLightPathTypeMap[OAIAO] = "OV->Acrylic->IV->Acrylic->OV";
  fLightPathTypeMap[Null] = "Light Path Type Un-initialized";
  std::cout << "Finished setting values" << std::endl;
}

void LightPathCalculator::ResetValues() {
  // Reset the private member variables that are modified by CalcByPosition by applying non-physical dummy values

  fLightPathType = Null;

  fIncidentVecOnPMT.SetXYZ(0.0, 0.0, 0.0);
  fInitialLightVec.SetXYZ(0.0, 0.0, 0.0);

  fStartPos.SetXYZ(-1e-6, -1e-6, -1e-6);
  fEndPos.SetXYZ(-1e-6, -1e-6, -1e-6);
  fLightPathEndPos.SetXYZ(-1e-6, -1e-6, -1e-6);
  fStartingRegion = "";

  fPathPrecision = -10.0;
  fIsTIR = false;
  fWithinError = false;
  fStraightLine = false;

  fPointOnAcrylic1st.SetXYZ(-1e-6, -1e-6, -1e-6);
  fPointOnAcrylic2nd.SetXYZ(-1e-6, -1e-6, -1e-6);
  fPointOnAcrylic3rd.SetXYZ(-1e-6, -1e-6, -1e-6);
  fPointOnAcrylic4th.SetXYZ(-1e-6, -1e-6, -1e-6);

  fDistInIV = -10.0;
  fDistInAcrylic = -10.0;
  fDistInOV = -10.0;

  fWavelength = 0.0;

  fSolidAngle = -100.0;
  fCosThetaAvg = -2.0;

  fFresnelTCoeff = -10.0;
  fFresnelRCoeff = -10.0;
}

void LightPathCalculator::CalcByPosition(TVector3 eventPos, TVector3 pmtPos, double wavelength, double precision) {
  ResetValues();
  std::cout << "Starting calcbyposition" << std::endl;
  // Check if refractive indices tables were properly loaded
  if ((fIVRefIndexGraph.GetN() * fAcrylicRefIndexGraph.GetN() * fOVRefIndexGraph.GetN() == 0) &&
      (fIVRefIndex * fAcrylicRefIndex * fOVRefIndex == 0)) {
    std::cout << "LPC::CalcByPosition: One or more refractive indices were not properly assigned. Check if "
                 "constructor was called before CalcByPosition"
              << std::endl;
    return;
  }

  // Check that none of the start or end positions are 'nan or 'inf'
  Double_t startPosMag = eventPos.Mag();
  if (std::isnan(startPosMag) || std::isinf(startPosMag)) {
    std::cout << "LightPathCalculator::CalcByPosition: The start position is nan/inf, ensure you use a valid start "
                 "position!"
              << std::endl;
    return;
  }
  Double_t endPosMag = pmtPos.Mag();
  if (std::isnan(endPosMag) || std::isinf(endPosMag)) {
    std::cout
        << "LightPathCalculator::CalcByPosition: The end position is nan/inf, ensure you use a valid end position!"
        << std::endl;
    return;
  }

  // Set the start and end position requirements of the path, and tolerance
  fStartPos = eventPos;
  fEndPos = pmtPos;
  fPathPrecision = precision;
  fInitialLightVec = (fEndPos - fStartPos).Unit();
  fWavelength = wavelength;

  // Determine the starting region
  fStartingRegion = DetermineRegion(eventPos);
  std::cout << "Starting region: " << fStartingRegion << std::endl;
  // Check to make sure starting region was determined correctly
  if (fStartingRegion == "NULL") {
    std::cout << "LightPathCalculator::DetermineRegion: Starting region not found correctly with position ("
              << eventPos.X() << ", " << eventPos.Y() << ", " << eventPos.Z() << ")" << std::endl;
    return;
  }

  if (fUseOpticsTables) {
    // Calculate the refractive indices for the given wavlength
    fIVRefIndex = InterpolateTGraph(fIVRefIndexGraph, wavelength);
    fAcrylicRefIndex = InterpolateTGraph(fAcrylicRefIndexGraph, wavelength);
    fOVRefIndex = InterpolateTGraph(fOVRefIndexGraph, wavelength);
  }  // Otherwise, stick with the user-loaded values

  // Check the refractive indices are set properly
  if (std::isnan(fIVRefIndex * fAcrylicRefIndex * fOVRefIndex) ||
      std::isinf(fIVRefIndex * fAcrylicRefIndex * fOVRefIndex) || (fIVRefIndex * fAcrylicRefIndex * fOVRefIndex == 0)) {
    std::cout << "LightPathCalculator::CalcByPosition: Refractive indices are not set properly" << std::endl;
    return;
  }

  // Check to see if the straight line approximation is required
  if (fPathPrecision == 0.0) {
    // Set the straight line path boolean (true)
    fStraightLine = true;
    // Path calculation without any refraction
    PathCalculation();
    return;
  } else {
    // Explicitly set straight line path flag to false and continue with calculation with refraction
    fStraightLine = false;
  }

  // The path result, whether it was calculated successfully, within error and no TIR, (true) or it failed (false)
  Bool_t pathResult = false;

  // Check to see if the starting position is in the IV region...
  if (fStartingRegion == "IV") {
    std::cout << "Starting IV Calc " << std::endl;
    pathResult = CalculateDistancesIV(fStartPos, fEndPos);
  }
  // ...or inside the acrylic (rare cases - return a straight line)...
  else if (fStartingRegion == "Acrylic") {
    std::cout << "LightPathCalculator::CalcByPosition: Start Position is in the Acrylic. LPC does not have this "
                 "functionality "
                 "yet."
              << std::endl;
    pathResult = false;
    return;
  }
  // ...or in the OV region outside.
  else if (fStartingRegion == "OV") {
    std::cout
        << "LightPathCalculator::CalcByPosition: Start Position is in the OV. LPC does not have this functionality "
           "yet."
        << std::endl;
    pathResult = false;
    return;
  } else if (fStartingRegion == "Past PMT") {
    pathResult = false;
    std::cout << "LightPathCalculator::CalcByPosition: Event Position is beyond the PMT volume" << std::endl;
    return;
  } else {
    pathResult = false;
    std::cout << "LightPathCalculator::CalcByPosition: Impossible start region occurred: " << fStartingRegion
              << std::endl;
    return;
  }
  std::cout << "Path result: " << pathResult << std::endl;
  if (!pathResult) {
    // Check to see if bad path result is due to TIR
    if (fIsTIR) {
      // If it is due to TIR, calculate the path without assuming refraction
      fStraightLine = true;
      if (fStartingRegion == "IV") {
        pathResult = CalculateDistancesIV(fStartPos, fEndPos);
      } else if (fStartingRegion == "Acrylic") {
        std::cout << "LightPathCalculator::CalcByPosition: Start Position is in the Acrylic. LPC does not have this "
                     "functionality "
                     "yet."
                  << std::endl;
        return;
      } else if (fStartingRegion == "OV") {
        std::cout << "LightPathCalculator::CalcByPosition: Start Position is in the OV. LPC does not have this "
                     "functionality "
                     "yet."
                  << std::endl;
        return;
        // pathResult = CalculateDistancesOV(fStartPos, fEndPos);
      } else if (fStartingRegion == "Past PMT") {
        std::cout << "LightPathCalculator::CalcByPosition: Event Position is beyond the PMT volume" << std::endl;
        return;
      } else {
        std::cout << "LightPathCalculator::CalcByPosition: Impossible start region occurred: " << fStartingRegion
                  << std::endl;
        return;
      }
    }
    // If not, it's because the error was greater than the user-specified precision. Bummer :(
  }
  std::cout << "Finished iv calc" << std::endl;
  return;
}

Bool_t LightPathCalculator::CalculateDistancesIV(const TVector3& eventPos, const TVector3& pmtPos) {
  fLightPathType = IAO;  // IV->Acrylic->OV->PMT

  PathCalculation();

  fDistInIV = (fPointOnAcrylic1st - fStartPos).Mag();
  fDistInAcrylic = (fPointOnAcrylic2nd - fPointOnAcrylic1st).Mag();
  fDistInOV = (fLightPathEndPos - fPointOnAcrylic2nd).Mag();
  if (fIsTIR || !fWithinError)
    return false;
  else
    return true;
}

Bool_t LightPathCalculator::CalculateDistancesOV(const TVector3& eventPos, const TVector3& pmtPos) {
  std::cout << "LightPathCalculator::CalculateDistancesOV: Not yet implemented" << std::endl;
  return false;
}  // Placeholder function

void LightPathCalculator::PathCalculation() {
  if (fStartingRegion == "IV") {
    // Set Light Path Type IV->Acrylic->OV->PMT
    fLightPathType = IAO;

    // Find the point where the light path intersects the inner edge of the acrylic from the IV.
    fPointOnAcrylic1st = IntersectAcrylic(fStartPos, fInitialLightVec, false);
    std::cout << "first point: " << fPointOnAcrylic1st.X() << ", " << fPointOnAcrylic1st.Y() << ", "
              << fPointOnAcrylic1st.Z() << std::endl;

    // Determine the refracted vector after refraction with the inner edge of the acrylic.
    TVector3 vecOne = PathRefraction(fInitialLightVec.Unit(), fPointOnAcrylic1st, fIVRefIndex, fAcrylicRefIndex);
    std::cout << "vecOne: " << vecOne.X() << ", " << vecOne.Y() << ", " << vecOne.Z() << std::endl;

    // Find the point where the light path intersects the outer edge of the acrylic from the inner edge.
    fPointOnAcrylic2nd = IntersectAcrylic(fPointOnAcrylic1st, vecOne, true);
    std::cout << "second point: " << fPointOnAcrylic2nd.X() << ", " << fPointOnAcrylic2nd.Y() << ", "
              << fPointOnAcrylic2nd.Z() << std::endl;

    // Determine the refracted vector after refraction with the outer edge of the acrylic.
    TVector3 vecTwo = PathRefraction(vecOne.Unit(), fPointOnAcrylic2nd, fAcrylicRefIndex, fOVRefIndex);
    std::cout << "vecTwo: " << vecTwo.X() << ", " << vecTwo.Y() << ", " << vecTwo.Z() << std::endl;

    // Work backwards from the PMT position to determine where vecTwo should originate from the outer edge of the
    // acrylic to hit the designated end pos.
    TVector3 calcPointOnAcrylic2nd = IntersectAcrylic(fEndPos, (-1 * vecTwo), true);
    std::cout << "calc point2nd: " << calcPointOnAcrylic2nd.X() << ", " << calcPointOnAcrylic2nd.Y() << ", "
              << calcPointOnAcrylic2nd.Z() << std::endl;

    // Difference between the fPointOnAcrylic2nd and calcPointOnAcrylic2nd is the ERROR
    fLightPathEndPos = fEndPos + (fPointOnAcrylic2nd - calcPointOnAcrylic2nd);
    std::cout << "calc end pos: " << fLightPathEndPos.X() << ", " << fLightPathEndPos.Y() << ", "
              << fLightPathEndPos.Z() << std::endl;

    // Check if the error is within the user-requested precision. If not, set the path result to be false.
    if ((fPointOnAcrylic2nd - calcPointOnAcrylic2nd).Mag() < fPathPrecision)
      fWithinError = true;
    else
      fWithinError = false;
    std::cout << "within error: " << fWithinError << std::endl;

  } else if (fStartingRegion == "Acrylic") {
    std::cout << "LightPathCalculator::PathCalculation: This functionality is not yet available." << std::endl;
    return;
  } else if (fStartingRegion == "OV") {
    std::cout << "LightPathCalculator::PathCalculation: This functionality is not yet available." << std::endl;
    return;
  } else if (fStartingRegion == "Past PMT") {
    std::cout << "LightPathCalculator::PathCalculation: Starting position is beyond the PMTs." << std::endl;
    return;
  } else {
    std::cout << "LightPathCalculator::PathCalculation: Impossible conditional block occured with start region "
              << fStartingRegion << std::endl;
    return;
  }
}

TVector3 LightPathCalculator::IntersectAcrylic(const TVector3& initPos, const TVector3& initDir,
                                               const bool& outerEdge) {
  Double_t startingRho =
      TMath::Sqrt((initPos.X() * initPos.X()) +
                  (initPos.Y() * initPos.Y()));  // TMath::Sqrt(std::pow(initPos.X(), 2) + std::pow(initPos.Y(), 2));
  TVector3 intersectionPoint;
  Double_t offset;
  Double_t xCross;
  Double_t yCross;
  Double_t zCross;
  Double_t acrylicCylRad = fIVCylRadius;
  Double_t acrylicCapRad = fIVCapRadius;

  if (outerEdge) {
    // If true, then we are trying to find the intersection with the outer edge of the AV. So adjust the radii by the
    // acrylic thickness
    acrylicCylRad += fAcrylicThickness;
    acrylicCapRad += fAcrylicThickness;
  }

  // First, check the cases where the direction is entirely vertical. Those paths won't cross the acrylic boundaries, so
  // the code after this conditional won't work
  if (initDir.X() == 0 && initDir.Y() == 0) {
    if (initDir.Z() < 0) {
      // Crossing occurs at the lower spherical cap
      zCross = -TMath::Sqrt((fIVCapRadius * fIVCapRadius) - (initPos.X() * initPos.X()) - (initPos.Y() * initPos.Y()));
    } else if (initDir.Z() > 0) {
      // Crossing occurs at the upper spherical cap
      zCross = TMath::Sqrt((fIVCapRadius * fIVCapRadius) - (initPos.X() * initPos.X()) - (initPos.Y() * initPos.Y()));
    }
    // Because the path is entirely vertical, the x, y values at the intersection point don't change.
    intersectionPoint.SetXYZ(initPos.X(), initPos.Y(), zCross);
    std::cout << "vertical point: " << initPos.X() << ", " << initPos.Y() << ", " << zCross << std::endl;
    return intersectionPoint;
  }

  // The path as a function of time is described as fStartPos + fInitialLightVec * t
  // First, we check if the cylindrical radius, rho(t), crosses the acrylic when the height is below the acrylic
  // cylinder height. rho(t) = sqrt((rho_0x^2 + v_0x^2 * t)^2 + (rho_0y^2 + v_0y^2 * t)^2) = fIVCylRadius. This becomes
  // a quadratic equation in t, with the following coefficients

  Double_t aCoeff = (initDir.X() * initDir.X()) +
                    (initDir.Y() * initDir.Y());  // std::pow(initDir.X(), 2) + std::pow(initDir.Y(), 2);
  Double_t bCoeff = 2 * (initPos.X() * initDir.X() + initPos.Y() * initDir.Y());
  Double_t cCoeff = (startingRho * startingRho) -
                    (acrylicCylRad * acrylicCylRad);  // std::pow(startingRho, 2) - std::pow(acrylicCylRad, 2);

  Double_t tCross;
  if (DetermineRegion(initPos) == "OV") {
    tCross = (-1 * bCoeff - TMath::Sqrt((bCoeff * bCoeff) - (4 * aCoeff * cCoeff))) / (2 * aCoeff);
  } else {
    tCross = (-1 * bCoeff + TMath::Sqrt((bCoeff * bCoeff) - (4 * aCoeff * cCoeff))) / (2 * aCoeff);
  }

  /*Double_t tCross2 =
      (-1 * bCoeff - TMath::Sqrt((bCoeff * bCoeff) - (4 * aCoeff * cCoeff))) /
      (2 * aCoeff);
  Double_t tCross;
  if (tCross2 > 0) tCross = tCross1;
  else if (tCross1 > 0) tCross = tCross2;
  else {
    std::cout << "LightPathCalculator::IntersectAcrylic: Reached a impossible conditional block" << std::endl;
    return TVector3(0.0, 0.0, 0.0);
  }*/
  // The procedure above has two possible roots, which is especially relevant when the initial point for this function
  // is outside of the cylinder. We pick the smallest positive root, corresponding to the first intersection

  zCross = initPos.Z() + (initDir.Z() * tCross);

  std::cout << aCoeff << ", " << bCoeff << ", " << cCoeff << ", " << tCross << ", " << zCross << std::endl;
  // Checking z crossing height
  if (TMath::Abs(zCross) <= fIVCylHeight) {
    std::cout << "made it here " << std::endl;
    // The crossing height is within the cylinder part of the acrylic. This makes life a lot easier.
    xCross = initPos.X() + (initDir.X() * tCross);
    yCross = initPos.Y() + (initDir.Y() * tCross);
    std::cout << xCross << ", " << yCross << ", " << zCross << std::endl;
    intersectionPoint.SetXYZ(xCross, yCross, zCross);
    return intersectionPoint;
  }
  // Otherwise, the crossing will occur at one of the spherical caps. Redo this calculation with the spherical radius
  // instead, but need to shift the zero of the coordinate system to be the center of the cylinder
  else if (zCross > fIVCylHeight) {
    std::cout << "made it here 2 " << std::endl;
    // Crossing top cap
    offset = fIVCapOffset;
  } else if (zCross < fIVCylHeight) {
    std::cout << "made it here 3" << std::endl;
    // Crossing bot cap
    offset = -1 * fIVCapOffset;
  } else {
    // This should NOT be reached, but have a check anyway. If a path is horizontal and starts within the cylinder, it
    // should end within the cylinder, meaning the first conditional is satisfied.
    std::cout << "LightPathCalculator::IntersectAcrylic: Reached impossible conditional block. Start Position is ("
              << initPos.X() << ", " << initPos.Y() << ", " << initPos.Z() << "), and initial direction is ("
              << initDir.X() << ", " << initDir.Y() << ", " << initDir.Z() << std::endl;
    return TVector3(0, 0, 0);
  }

  std::cout << "starting modification" << std::endl;
  TVector3 modStartPos;
  modStartPos.SetXYZ(initPos.X(), initPos.Y(),
                     initPos.Z() + offset);  // Moving the z position to match the new relevant origins

  std::cout << "modified position: " << modStartPos.X() << ", " << modStartPos.Y() << ", " << modStartPos.Z()
            << std::endl;
  // Now, we repeat the process, but looking for the moment when the path crosses the spherical cap
  Double_t startingR = modStartPos.Mag();
  aCoeff = (initDir.X() * initDir.X()) + (initDir.Y() * initDir.Y()) + (initDir.Z() * initDir.Z());
  bCoeff = 2 * (modStartPos.X() * initDir.X() + modStartPos.Y() * initDir.Y() + modStartPos.Z() * initDir.Z());
  cCoeff = (startingR * startingR) - (acrylicCapRad * acrylicCapRad);

  if (DetermineRegion(initPos) == "OV") {
    tCross = (-1 * bCoeff - TMath::Sqrt((bCoeff * bCoeff) - (4 * aCoeff * cCoeff))) / (2 * aCoeff);
  } else {
    tCross = (-1 * bCoeff + TMath::Sqrt((bCoeff * bCoeff) - (4 * aCoeff * cCoeff))) / (2 * aCoeff);
  }
  xCross = modStartPos.X() + (initDir.X() * tCross);
  yCross = modStartPos.Y() + (initDir.Y() * tCross);
  zCross = modStartPos.Z() + (initDir.Z() * tCross);

  std::cout << aCoeff << ", " << bCoeff << ", " << cCoeff << std::endl
            << tCross << ", " << xCross << "," << yCross << ", " << zCross << std::endl;

  intersectionPoint.SetXYZ(xCross, yCross, zCross - offset);  // Undoing the offset from shifting the center
  return intersectionPoint;
}

TVector3 LightPathCalculator::PathRefraction(const TVector3& incidentVec, const TVector3& incidentPoint,
                                             const Double_t incRIndex, const Double_t refRIndex) {
  if (fStraightLine) {
    std::cout << "using straight line pr" << std::endl;
    // No need to calculate refraction
    return incidentVec;
  }

  // Find the normal vector at the provided intersection point
  TVector3 normVec = GetNormalVector(incidentPoint);
  std::cout << "norm vec: " << normVec.X() << ", " << normVec.Y() << ", " << normVec.Z() << std::endl;

  const Double_t ratioRI = incRIndex / refRIndex;
  std::cout << ratioRI << std::endl;
  const Double_t cosTheta1 = (normVec.Unit()).Dot(incidentVec.Unit());  // Incident angle [Snell's Law]
  std::cout << "costheta1: " << cosTheta1 << std::endl;
  if (cosTheta1 > 1.0) {
    std::cout << "LightPathCalculator::PathRefraction: cosTheta1 is greater than 1 somehow" << std::endl;
    PrintPrivateVariables();
  }
  const Double_t cosTheta2 =
      TMath::Sqrt(1 - (ratioRI * ratioRI) * (1 - (cosTheta1 * cosTheta1)));  // Refracted Angle (Snell's Law)
  std::cout << "costheta2: " << cosTheta2 << std::endl;

  // Initialize the refracted photon vector
  TVector3 refractedVec(0.0, 0.0, 0.0);

  // Check for Total Internal Reflection (TIR) (equivalent to a 'nan' radicand i.e. a negative squareroot)
  if (std::isnan(cosTheta2)) {
    fIsTIR = true;
    // Set the refracted vec to the straight equivalent
    refractedVec = incidentVec;
  }

  else {
    refractedVec = (ratioRI * incidentVec) + ((ratioRI * cosTheta1) - cosTheta2) * -1 * normVec;
  }

  std::cout << "refractedVec: " << refractedVec.X() << ", " << refractedVec.Y() << ", " << refractedVec.Z()
            << std::endl;
  // Ensure the refracted vector is unit normalised
  return refractedVec.Unit();
}

TVector3 LightPathCalculator::GetNormalVector(const TVector3& point, double tolerance) {
  // GetNormalVector should only be called using points that are on the acrylic or very close to the acrylic.

  // Check that the point is not beyond the PMTs before proceeding
  if (DetermineRegion(point) == "Past PMT") {
    std::cout << "LightPathCalculator::GetNormalVector: Provided point is beyond the PMT region." << std::endl;
    PrintPrivateVariables();
    return TVector3(0.0, 0.0, 0.0);
  }

  double cylRad = TMath::Sqrt((point.X() * point.X()) + (point.Y() * point.Y()));
  double sphereRad = point.Mag();
  double height = point.Z();
  TVector3 normalVector;

  // The difference between intersection with the inner and outer edges only matter for this tolerance check.
  // Because the Eos vesssel is cylindrical and spherical the normal vectors will be equal on the inner and outer
  // edges.
  if ((TMath::Abs(cylRad - fIVCylRadius) <= tolerance) ||
      (TMath::Abs(cylRad - (fIVCylRadius + fAcrylicThickness)) <= tolerance) ||
      (TMath::Abs(sphereRad - fIVCapRadius) <= tolerance) ||
      (TMath::Abs(sphereRad - (fIVCapRadius + fAcrylicThickness)) <=
       tolerance)) {  // Checks to see if point is on or close to the inner or outer edge of the acrylic
    if (TMath::Abs(height) <= fIVCylHeight) {
      normalVector.SetXYZ(point.X(), point.Y(), 0.0);
      normalVector = normalVector.Unit();
    } else if (height > 0) {
      normalVector.SetXYZ(point.X(), point.Y(), point.Z() + fIVCapOffset);
      normalVector = normalVector.Unit();
    } else if (height < 0) {
      normalVector.SetXYZ(point.X(), point.Y(), point.Z() - fIVCapOffset);
      normalVector = normalVector.Unit();
    } else {
      std::cout << "LightPathCalculator::GetNormalVector: Impossible conditional block reached. Point is (" << point.X()
                << ", " << point.Y() << ", " << point.Z() << ")" << std::endl;
      return TVector3(0.0, 0.0, 0.0);
    }
  } else {
    std::cout << "LightPathCalculator::GetNormalVector: Provided point is too far from either acrylic edge."
              << std::endl;
    return TVector3(0.0, 0.0, 0.0);
  }
  return normalVector;
}

std::string LightPathCalculator::DetermineRegion(const TVector3& point) {
  double cylRad = TMath::Sqrt((point.X() * point.X()) + (point.Y() * point.Y()));
  double sphereRad = point.Mag();
  double height = point.Z();

  if (cylRad < fIVCylRadius) {
    if (TMath::Abs(height) < fIVCylHeight) {
      return "IV";
    } else if (sphereRad < fIVCapRadius) {
      return "IV";  // Technically satisfied by previous, but nice to separate
    } else if (sphereRad < (fIVCapRadius + fAcrylicThickness)) {
      return "Acrylic";
    } else {
      if (height > 0) {
        if (height < fTopPMTRadius) {
          return "OV";
        } else {
          return "Past PMT";
        }
      } else if (height < 0) {
        if (TMath::Abs(height) > fBotPMTRadius) {
          return "OV";
        } else {
          return "Past PMT";
        }
      } else {
        std::cout << "LightPathCalculator::DetermineRegion: Reached a forbidden conditional. Point is (" << point.X()
                  << ", " << point.Y() << ", " << point.Z() << ")" << std::endl;
        return "NULL";
      }
    }
  } else if (cylRad < (fIVCylRadius + fAcrylicThickness)) {
    if (TMath::Abs(height) < fIVCylHeight) {
      return "Acrylic";
    } else {
      if (height > 0) {
        if (height < fTopPMTRadius) {
          return "OV";
        } else {
          return "Past PMT";
        }
      } else if (height < 0) {
        if (TMath::Abs(height) > fBotPMTRadius) {
          return "OV";
        } else {
          return "Past PMT";
        }
      } else {
        std::cout << "LightPathCalculator::DetermineRegion: Reached a forbidden conditional. Point is (" << point.X()
                  << ", " << point.Y() << ", " << point.Z() << ")" << std::endl;
        return "NULL";
      }
    }
  } else if (cylRad <= fBarrelPMTRadius) {
    if (TMath::Abs(height) < fBarrelPMTHeight) {
      return "OV";
    } else {
      if (height > 0) {
        if (height < fTopPMTRadius) {
          return "OV";
        } else {
          return "Past PMT";
        }
      } else if (height < 0) {
        if (TMath::Abs(height) > fBotPMTRadius) {
          return "OV";
        } else {
          return "Past PMT";
        }
      } else {
        std::cout << "LightPathCalculator::DetermineRegion: Reached a forbidden conditional. Point is (" << point.X()
                  << ", " << point.Y() << ", " << point.Z() << ")" << std::endl;
        return "NULL";
      }
    }
  } else {
    return "Past PMT";
  }
}

Double_t LightPathCalculator::InterpolateTGraph(const TGraph& graph, const Double_t wl) {
  if (graph.GetN() == 0) {
    // Cannot use a graph that's empty!
    Log::Die("LightPathCalculator::InterpolateTGraph: Cannot interpolate an empty TGraph.");
  }

  Double_t* wl_values = graph.GetX();
  if (wl < wl_values[0]) return wl_values[0];
  if (wl > wl_values[graph.GetN() - 1]) return wl_values[graph.GetN() - 1];
  return graph.Eval(wl);
}

void LightPathCalculator::PrintPrivateVariables() {
  // This function is used to print out all the various private member variables, useful when an "impossible" or
  // "improper" condition is encountered
  std::cout << "-----------------------------------------------" << std::endl
            << "Light Path Type: " << fLightPathTypeMap[fLightPathType] << std::endl
            << "Starting position: (" << fStartPos.X() << ", " << fStartPos.Y() << ", " << fStartPos.Z() << ")"
            << std::endl
            << "Starting region: " << fStartingRegion << std::endl
            << "Ending position (user input): " << fEndPos.X() << ", " << fEndPos.Y() << ", " << fEndPos.Z() << ")"
            << std::endl
            << "Ending position (calculated): " << fLightPathEndPos.X() << ", " << fLightPathEndPos.Y() << ", "
            << fLightPathEndPos.Z() << ")" << std::endl
            << "-----------------------------------------------" << std::endl
            << "Refractive indices: " << fIVRefIndex << ", " << fAcrylicRefIndex << ", " << fOVRefIndex << std::endl
            << "Materials: " << fIVMaterial << ", " << fAcrylicMaterial << ", " << fOVMaterial << std::endl
            << "-----------------------------------------------" << std::endl
            << "Initial Light Vector: (" << fInitialLightVec.X() << ", " << fInitialLightVec.Y() << ", "
            << fInitialLightVec.Z() << ")" << std::endl
            << "First Point on Acrylic: (" << fPointOnAcrylic1st.X() << ", " << fPointOnAcrylic1st.Y() << ", "
            << fPointOnAcrylic1st.Z() << ")" << std::endl
            << "Second Point on Acrylic: (" << fPointOnAcrylic2nd.X() << ", " << fPointOnAcrylic2nd.Y() << ", "
            << fPointOnAcrylic2nd.Z() << ")" << std::endl
            << "Third Point on Acrylic: (" << fPointOnAcrylic3rd.X() << ", " << fPointOnAcrylic3rd.Y() << ", "
            << fPointOnAcrylic3rd.Z() << ")" << std::endl
            << "Fourth Point on Acrylic: (" << fPointOnAcrylic4th.X() << ", " << fPointOnAcrylic4th.Y() << ", "
            << fPointOnAcrylic4th.Z() << ")" << std::endl
            << "Distances in IV, Acrylic, OV: " << fDistInIV << ", " << fDistInAcrylic << ", " << fDistInOV << std::endl
            << "-----------------------------------------------" << std::endl;
}

}  // namespace RAT
