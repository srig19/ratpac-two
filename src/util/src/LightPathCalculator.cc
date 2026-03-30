#include "RAT/LightPathCalculator.hh"

ClassImp(RAT::LightPathCalculator)

    namespace RAT {
  // Constructor where user assigns the material for each region. Refractive indices are pulled from optics tables
  LightPathCalculator::LightPathCalculator(std::string materialIV, std::string materialAcrylic,
                                           std::string materialOV) {
    fUseOpticsTables = true;

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

    try {
      opticsIV = DB::Get()->GetLink("OPTICS", materialIV);
      optionIV = opticsIV->GetS("RINDEX_option");
      valOneIV = opticsIV->GetDArray("RINDEX_value1");
      valTwoIV = opticsIV->GetDArray("RINDEX_value2");
    } catch (DBNotFoundError& e) {
      Log::Die("LightPathCalculator: Failed to retrieve optics information for IV: " + materialIV);
    }
    try {
      opticsAcrylic = DB::Get()->GetLink("OPTICS", materialAcrylic);
      optionAcrylic = opticsAcrylic->GetS("RINDEX_option");
      valOneAcrylic = opticsAcrylic->GetDArray("RINDEX_value1");
      valTwoAcrylic = opticsAcrylic->GetDArray("RINDEX_value2");
    } catch (DBNotFoundError& e) {
      Log::Die("LightPathCalculator: Failed to retrieve optics information for Acrylic: " + materialAcrylic);
    }
    try {
      opticsOV = DB::Get()->GetLink("OPTICS", materialOV);
      optionOV = opticsOV->GetS("RINDEX_option");
      valOneOV = opticsOV->GetDArray("RINDEX_value1");
      valTwoOV = opticsOV->GetDArray("RINDEX_value2");
    } catch (DBNotFoundError& e) {
      Log::Die("LightPathCalculator: Failed to retrieve optics information for OV: " + materialOV);
    }

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
        Log::Die("LightPathCalculator: Optics bable vectors are not the same size.");
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
          fIVRefIndexGraph.SetPoint(pointNum, (2 * TMath::Pi() * hbarc) / valOne[j], valTwo[j]);
        } else if (i == 1) {
          fAcrylicRefIndexGraph.SetPoint(pointNum, (2 * TMath::Pi() * hbarc) / valOne[j], valTwo[j]);
        } else if (i == 2) {
          fOVRefIndexGraph.SetPoint(pointNum, (2 * TMath::Pi() * hbarc) / valOne[j], valTwo[j]);
        }
        pointNum++;
      }
    }

    SetValues();
  }

  // Constructor where user chooses a single refractive index for each material
  LightPathCalculator::LightPathCalculator(Double_t refIV, Double_t refAcrylic, Double_t refOV) {
    fUseOpticsTables = false;

    fIVRefIndex = refIV;
    fAcrylicRefIndex = refAcrylic;
    fOVRefIndex = refOV;

    SetValues();
  }

  // Common aspects of constructor between both constructor options above
  void LightPathCalculator::SetValues() {
    // Setting geometry specific values
    DBLinkPtr geoPtr = DB::Get()->GetLink("GEO", "eos_vessel");
    fIVCylRadius = geoPtr->GetD("r_max");
    fIVCylHeight = geoPtr->GetD("size_z");  // This is HALF the height
    fIVCapRadius = geoPtr->GetD("top_radius") Double_t fIVCapOffset = TMath::Abs(fIVCapRadius - fIVCylHeight);
    fAcrylicThickness = 25.4;  // Acrylic is 1 inch thick
    fBarrelPMTRadius = 1049.6;
    fBarrelPMTHeight = 658.8;  // This is HALF the height
    fTopPMTRadius = 1070.0;
    fBotPMTRadius = 1070.0;
  }

  void LightPathCalculator::ResetValues() {
    // Reset the private member variables that are modified by CalcByPosition by applying non-physical dummy values

    fPathPrecision = -10.0;

    fIncidentVecOnPMT.SetXYZ(0.0, 0.0, 0.0);
    fInitialLightVec.SetXYZ(0.0, 0.0, 0.0);

    fStartPos.SetXYZ(-10.0e6, -10.0e6, -10.0e6);
    fEndPos.SetXYZ(-10.0e6, -10.0e6, -10.0e6);
    fLightPathEndPos.SetXYZ(-10.0e6, -10.0e6, -10.0e6);

    fIsTIR = false;
    fResvHit = false;
    fStraightLine = false;

    fPointOnAcrylic1st.SetXYZ(-10.0e6, -10.0e6, -10.0e6);
    fPointOnAcrylic2nd.SetXYZ(-10.0e6, -10.0e6, -10.0e6);
    fPointOnAcrylic3rd.SetXYZ(-10.0e6, -10.0e6, -10.0e6);
    fPointOnAcrylic4th.SetXYZ(-10.0e6, -10.0e6, -10.0e6);

    fLightPathType = Null;

    fEnergy = 0.0;

    fSolidAngle = -100.0;
    fCosThetaAvg = -2.0;

    fFresnelTCoeff = -10.0;
    fFresnelRCoeff = -10.0;
  }

  void LightPathCalculator::CalcByPosition(TVector3 eventPos, TVector3 pmtPos, double wavelength, double precision) {
    ResetValues();

    // Check if refractive indices tables were properly loaded
    if ((fIVRefIndexGraph.GetN() * fAcrylicRefIndexGraph.GetN() * fOVRefIndexGraph.GetN() == 0) &&
        (fIVRefIndex * fAcrylicRefIndex * fOVRefIndex == 0)) {
      Log::Die(
          "LPC::CalcByPosition: One or more refractive indices were not properly assigned. Check if constructor was "
          "called before CalcByPosition");
    }

    // Check that none of the start or end positions are 'nan or 'inf'
    Double_t startPosMag = eventPos.Mag();
    if (std::isnan(startPosMag) || std::isinf(startPosMag)) {
      Log::Die(
          "LightPathCalculator::CalcByPosition: The start position is nan/inf, ensure you use a valid start position!");
    }
    Double_t endPosMag = pmtPos.Mag();
    if (std::isnan(endPosMag) || std::isinf(endPosMag)) {
      Log::Die(
          "LightPathCalculator::CalcByPosition: The end position is nan/inf, ensure you use a valid end position!");
    }

    // Set the start and end position requirements of the path, and tolerance
    fStartPos = eventPos;
    fEndPos = pmtPos;
    fPathPrecision = precision;
    // Determine the starting region
    fStartingRegion = GetRegion(eventPos);

    if (fUseOpticsTables) {
      // Calculate the refractive indices for the given wavlength
      fIVRefIndex = InterpolateTGraph(fIVRefIndexGraph, wavelength);
      fAcrylicRefIndex = InterpolateTGraph(fAcrylicRefIndexGraph, wavelength);
      fOVRefIndex = InterpolateTGraph(fOVRefIndexGraph, wavelength);
    }  // Otherwise, stick with the user-loaded values

    // Check the refractive indices are set properly
    if (std::isnan(fIVRefIndex * fAcrylicRefIndex * fOVRefIndex) ||
        std::isinf(fIVRefIndex * fAcrylicRefIndex * fOVRefIndex) ||
        (fIVRefIndex * fAcrylicRefIndex * fOVRefIndex == 0)) {
      Log::Die("LightPathCalculator:: Refractive indices are not set properly");
    }

    // Check to see if the straight line approximation is required
    if (fPathPrecision == 0.0) {
      // Set the straight line path boolean (true)
      fStraightLine = true;
      // Path calculation without any refraction
      PathCalculation((fEndPos - fStartPos).Unit());
      return;
    } else {
      // Explicitly set straight line path flag to false and continue with calculation with refraction
      fStraightLine = false;
    }

    // If the start position magnitude is 0.0, then set it to something very small so that all the vector algebra works
    if (startPosMag == 0.0) {
      fStartPos.SetXYZ(1.0e-4, 0.0, 0.0);
    }
    if (endPosMag == 0.0) {
      fEndPos.SetXYZ(1.0e-4, 0.0, 0.0);
    }

    // The path result, whether it was calculated successfully (true) or it failed (false)
    Bool_t pathResult = false;

    // Check to see if the starting position is in the IV region...
    if (startRegion == "IV") {
      pathResult = CalculateDistancesIV(fStartPos, fEndPos);
    }
    // ...or inside the acrylic (rare cases - return a straight line)...
    else if (startRegion == "Acrylic") {
      pathResult = false;
    }
    // ...or in the OV region outside.
    else if (startRegion == "OV") {
      Log::Die(
          "LightPathCalculator::CalcByPosition: Start Position is in the OV. LPC does not have this functionality "
          "yet.");
      // pathResult = CalculateDistancesOV(fStartPos, fEndPos);
    } else if (startRegion == "Past PMT") {
      pathResult = false;
      Log::Die("LightPathCalculator::CalcByPosition: Event Position is beyond the PMT volume");
    } else {
      Log::Die("LightPathCalculator::CalcByPosition: Impossible start region occurred: " + startRegion);
    }

    // Check to see the status of the pathResult boolean.
    // If it is false, it is likely one of the above methods
    // either encountered Total Internal Reflection [TIR] (fIsTIR = true)
    // or the calculated end position of the path is far from the required
    // localityVal (fResvHit = true)

    /*if (pathResult == false){
        fStraightLine = true;

        // Safegaurd to keep the original total internal reflection and locality end position statuses
        Bool_t tmpTIR = false;
        Bool_t tmpResv = false;

        if ( fIsTIR == true ){ tmpTIR = true; }
        if ( fResvHit == true ){ tmpResv = true; }

    }*/
  }

  Bool_t CalculateDistancesIV(const TVector3& eventPos, const TVector3& pmtPos) {
    LightPathType = IAO;  // IV->Acrylic->OV->PMT
    fStartPos = eventPos;
    fEndPos = pmtPos;

    // xUnit - Points along the radial direction from the
    //         origin to the start position.
    // zUnit - Vector perpendicular to the plane
    //         containing the origin, start and pmt position.
    // yUnit - Vector perpendicular to the starting position vector,
    //         located in the plane of the start and pmt positions.
    TVector3 xUnit = fStartPos.Unit();
    TVector3 zUnit = (xUnit.Cross(fEndPos)).Unit();
    TVector3 yUnit = (zUnit.Cross(xUnit)).Unit();

    // The required target angle between the starting position vector
    // and the PMT radial vector. This is required for LightPathCalculator::RTSafe()
    fPMTTargetTheta = TMath::ATan2((fEndPos.Cross(fStartPos)).Mag(), fEndPos * fStartPos);

    // 'theta': the calculated angle between the starting position vector
    //  and the PMT radial vector. This is compared to the 'fPMTTargetTheta'
    //  value and minimised against in RTSafe
    //  ( and inherently, ::ThetaResidual, ::DThetaResidual and ::FuncD )
    Double_t theta = RTSafe(0.0, pi, 1.0e-3);

    // RTSafe implicitly computes the angles between each of the interface
    // points for the path. It accounts for the refractive effects between
    // each of these boundaries (see Theta1st(), Theta2nd() etc.)
    // So it knows about Snell's law and thus if an instance of total internal
    // reflection is encountered.

    // ...otherwise we can continue to define the refracted path intersection points.

    // Test for total internal reflection whilst calculating the below
    // ThetaX values to see if the fitted 'theta' value is invalid.
    fIsTIR = false;
  }

  Bool_t CalculateDistancesOV(const TVector3& eventPos, const TVector3& pmtPos) {
    return false;
  }  // Placeholder function

  void PathCalculation(const TVector3& initDir) {
    fInitialLightVec = initDir.Unit();

    if (fStartingRegion == "IV") {
      // Set Light Path Type IV->Acrylic->OV->PMT
      fLightPathType = IAO;

      // Find the point where the light path intersects the inner edge of the acrylic from the IV.
      fPointOnAcrylic1st = IntersectAcrylic(fStartPos, fInitialLightVec, false);

      // Determine the refracted vector after refraction with the inner edge of the acrylic.
      TVector vecOne = PathRefraction(fInitialLightVec.Unit(), )

    } else if (fStartingRegion == "Acrylic") {
      Log::Die("LightPathCalculator::PathCalculation: This functionality is not yet available.");
    } else if (fStartingRegion == "OV") {
      Log::Die("LightPathCalculator::PathCalculation: This functionality is not yet available.");
    } else if (fStartingRegion == "Past PMT") {
      Log::Die("LightPathCalculator::PathCalculation: Starting position is beyond the PMTs.");
    } else {
      Log::Die("LightPathCalculator::PathCalculation: Impossible conditional block occured with start region " +
               fStartingRegion);
    }
  }

  Double_t LightPathCalculator::ThetaResidual(const Double_t theta) {
    // Compute the difference between the PMT's fPMTTargetTheta
    // and ThetaX( theta ) for the current estimate of starting angle
    // from the source theta.
    // Here thetaX( theta ) is the final angle in the path type
    // e.g. for IV -> Acrylic -> OV; ThetaX = Theta3
    // e.g. for OV -> Acrylic -> IV -> Acrylic -> OV; ThetaX = Theta5

    Double_t retValue = 0.0;
    if (fLightPathType == IAO) {
      retValue = fPMTTargetTheta - (Theta3rd(theta) + Theta2nd(theta) + Theta1st(theta));
    } else if (fLightPathType == OAIAO) {
      // SHOULD NOT REACH HERE YET
      retValue =
          fPMTTargetTheta - (Theta5th(theta) + Theta4th(theta) + Theta3rd(theta) + Theta2nd(theta) + Theta1st(theta));
    } else {  // Shouldn't reach here, but handle just in case and return the light path type
      Log::Die("LightPathCalculator::ThetaResidual: Error! Unknown light path encountered with light path type: " +
               fLightPathTypeMap[fLightPathType])
    }
    return retValue;
  }

  Double_t LightPathCalculator::DThetaResidual(const Double_t theta) {
    // Compute the derivative in the difference between ... (see above)
    Double_t retDValue = 0.0;

    if (fLightPathType == IAO) {
      retDValue = (-1.0 * (DTheta3rd(theta) + DTheta2nd(theta) + DTheta1st(theta)));
    } else if (fLightPathType == OAIAO) {
      retDValue =
          (-1.0 * (DTheta5th(theta) + DTheta4th(theta) + DTheta3rd(theta) + DTheta2nd(theta) + DTheta1st(theta)));

    } else {  // Shouldn't reach here, but handle just in case and return the light path type
      Log::Die("LightPathCalculator::ThetaDResidual: Error! Unknown light path encountered with light path type: " +
               fLightPathTypeMap[fLightPathType])
    }
    return retDValue;
  }

  TVector3 LightPathCalculator::IntersectAcrylic(const TVector3& initPos, const TVector3& initDir,
                                                 const bool& outerEdge) {
    Double_t startingRho = TMath::Sqrt(pow(initPos.X(), 2) + pow(initPos.Y(), 2));
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

    // The path as a function of time is described as fStartPos + fInitialLightVec * t
    // First, we check if the cylindrical radius, rho(t), crosses the acrylic when the height is below the acrylic
    // cylinder height. rho(t) = sqrt((rho_0x^2 + v_0x^2 * t)^2 + (rho_0y^2 + v_0y^2 * t)^2) = fIVCylRadius This becomes
    // a quadratic equation in t, with the following coefficients

    Double_t aCoeff = pow(initDir.X(), 2) + pow(initDir.Y(), 2);
    Double_t bCoeff = 2 * (initPos.X() * initDir.X() + initPos.Y() * initDir.Y());
    Double_t cCoeff = pow(startingRho, 2) - pow(acrylicCylRad, 2);

    Double_t tCross = (-1 * bCoeff + TMath::Sqrt(pow(bCoeff, 2) - 4 * aCoeff * cCoeff)) / (2 * aCoeff);
    zCross = initPos.Z() + (initDir.Z() * tCross);
    // Checking z crossing height
    if (TMath::Abs(zCross) <= fIVCylHeight) {
      // The crossing hieght is within the cylinder part of the acrylic. This makes life a lot easier.
      xCross = initPos.X() + (initDir.X() * tCross);
      yCross = initPos.Y() + (initDir.Y() * tCross);
      intersectionPoint.SetXYZ(xCross, yCross, zCross);
      return intersectionPoint;
    }
    // Otherwise, the crossing will occur at one of the spherical caps. Redo this calculation with the spherical radius
    // instead, but need to shift the zero of the coordinate system to be the center of the cylinder
    else if (zCross > fIVCylHeight) {
      // Crossing top cap
      offset = fIVCapOffset;
    } else if (zCross < fIVCylHeight) {
      // Crossing bot cap
      offset = -1 * fIVCapOffset;
    } else {
      // This should NOT be reached, but have a check anyway. If a path is horizontal and starts within the cylinder, it
      // should end within the cylinder, meaning the first conditional is satisfied.
      Log::Die("LightPathCalculator::IVToInnerEdge: Reached impossible conditional block. Start Position is (" +
               fStartPos.X() + ", " + fStartPos.Y() + ", " + fStartPos.Z() + "), and initial direction is (" +
               fInitialLightVec.X() + ", " + fInitialLightVec.Y() + ", " + fInitialLightVec.Z());
    }

    TVector3 modStartPos;
    modStartPos.SetXYZ(initPos.X(), initPos.Y(),
                       initPos.Z() + offset);  // Moving the z position to match the new relevant origin

    // Now, we repeat the process, but looking for the moment when the path crosses the spherical cap
    Double_t startingR = modStartPos.Mag();
    aCoeff = pow(initDir.X(), 2) + pow(initDir.Y(), 2) + pow(initDir.Z(), 2);
    bCoeff = 2 * (modStartPos.X() * initDir.X() + modStartPos.Y() * initDir.Y() + modStartPos.Z() * initDir.Z());
    cCoeff = pow(startingR, 2) - pow(acrylicCapRad, 2);

    tCross = (-1 * bCoeff + TMath::Sqrt(pow(bCoeff, 2) - 4 * aCoeff * cCoeff)) / (2 * aCoeff);
    xCross = modStartPos.X() + (initDir.X() * tCross);
    yCross = modStartPos.Y() + (initDir.Y() * tCross);
    zCross = modStartPos.Z() + (initDir.Z() * tCross);

    intersectionPoint.SetXYZ(xCross, yCross, zCross - offset);  // Undoing the offset from shifting the center
    return intersectionPoint;
  }

  void LightPathCalculator::FuncD(Double_t theta, Double_t & funcVal, Double_t & dFuncVal) {
    funcVal = ThetaResidual(theta);
    dFuncVal = DThetaResidual(theta);
  }

  TVector3 LightPathCalculator::PathRefraction(const TVector3& incidentVec, const TVector3& incidentPoint,
                                               const Double_t incRIndex, const Double_t refRIndex) {
    if (fStraightLine)
      // No need to calculate refraction
      return incidentVec;

    // Find the normal vector at the provided intersection point
    TVector3 normVec = GetNormalVector(incidentPoint);

    const Double_t ratioRI = incRIndex / refRIndex;
    const Double_t cosTheta1 = normVec.Dot(-1.0 * incidentVec);  // Incident angle [Snell's Law]
    const Double_t cosTheta2 =
        TMath::Sqrt(1 - TMath::Power(ratioRI, 2) * (1 - TMath::Power(cosTheta1, 2)));  // Refracted Angle [Snell's Law]

    // Initialise the refracted photon vector
    TVector3 refractedVec(0.0, 0.0, 0.0);

    // Check for Total Internal Reflection (TIR) (equivalent to a 'nan' radicand i.e. a negative squareroot)
    if (std::isnan(cosTheta2)) {
      fIsTIR = true;
      // Set the refracted vec to the straight equivalent
      refractedVec = incidentVec;
    }

    // Define the refracted vector
    else if (cosTheta1 >= 0.0) {
      refractedVec = (ratioRI * incidentVec) + ((ratioRI * cosTheta1) - cosTheta2) * incidentSurfVec;
    }

    else {
      refractedVec = (ratioRI * incidentVec) - ((ratioRI * cosTheta1) - cosTheta2) * normVec;
    }

    // Ensure the refracted vector is unit normalised
    return refractedVec.Unit();
  }

  TVector3 GetNormalVector(const TVector3& point, double tolerance) {
    // GetNormalVector should only be called using points that are on the acrylic or very close to the acrylic.

    // Check that the point is not beyond the PMTs before proceeding
    if (GetRegion(point) == "Past PMT") {
      Log::Die("LightPathCalculator::GetNormalVector: Provided point is beyond the PMT region. Point is (" +
               to_string(point.X()) + ", " + to_string(point.Y()) + ", " + to_string(point.Z()) + ")");
    }

    double cylRad = TMath::Sqrt(pow(point.X(), 2) + pow(point.Y(), 2));
    double sphereRad = point.Mag();
    double height = point.Z();
    TVector3 normalVector;

    // The difference between intersection with the inner and outer edges only matter for this tolerance check.
    // Because the Eos vesssel is cylindrical and spherical the normal vectors will be equal on the inner and outer
    // edges.
    if ((TMath::Abs(cylRad - fIVCylRadius) <= tolerance) ||
        (TMath::Abs(cylRad - (fIVCylRadius + fAcrylicThickness)) <=
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
        Log::Die("LightPathCalculator::GetNormalVector: Impossible conditional block reached. Point is (" +
                 to_string(point.X()) + ", " + to_string(point.Y()) + ", " + to_string(point.Z()) + ")");
      }
    } else {
      Log::Die("LightPathCalculator::GetNormalVector: Provided point is too far from either acrylic edge.");
    }
    return normalVector;
  }

  std::string DetermineRegion(const TVector3& point) {
    double cylRad = TMath::Sqrt(point.X() * *2 + point.Y() * *2);
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
          if (height > fBotPMTRadius) {
            return "OV";
          } else {
            return "Past PMT";
          }
        } else {
          Log::Die("LightPathCalculator::GetRegion: Reached a forbidden conditional. Point is (" +
                   to_string(point.X()) + ", " + to_string(point.Y()) + ", " + to_string(point.Z()) + ")");
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
          if (height > fBotPMTRadius) {
            return "OV";
          } else {
            return "Past PMT";
          }
        } else {
          Log::Die("LightPathCalculator::GetRegion: Reached a forbidden conditional. Point is (" +
                   to_string(point.X()) + ", " + to_string(point.Y()) + ", " + to_string(point.Z()) + ")");
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
          if (height > fBotPMTRadius) {
            return "OV";
          } else {
            return "Past PMT";
          }
        } else {
          Log::Die("LightPathCalculator::GetRegion: Reached a forbidden conditional. Point is (" +
                   to_string(point.X()) + ", " + to_string(point.Y()) + ", " + to_string(point.Z()) + ")");
        }
      }
    } else {
      return "Past PMT";
    }
  }

}  // namespace RAT