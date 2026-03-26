#include "RAT/LightPathCalculator.hh"

ClassImp(RAT::LightPathCalculator)

    namespace RAT {
  // Constructor where user assigns the material for each region. Refractive indices are pulled from optics tables
  LightPathCalculator::LightPathCalculator(std::string materialIV, std::string materialAcrylic,
                                           std::string materialOV) {
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
    fIVRefIndex = refIV;
    fAcrylicRefIndex = refAcrylic;
    fOVRefIndex = refOV;

    SetValues();
  }

  // Common aspects of constructor between both constructor options above
  void LightPathCalculator::SetValues() {
    // Setting geometry specific values
    DBLinkPtr geoPtr = DB::Get()->GetLink("GEO", "eos_vessel");
    Double_t fIVCylRadius = geoPtr->GetD("r_max");
    Double_t fIVCylHeight = geoPtr->GetD("size_z");  // This is HALF the height
    Double_t fIVCapRadius = geoPtr->GetD("top_radius") Double_t fIVCapHeight = geoPtr->GetD("top_height");
    Double_t fAcrylicThickness = 25.4;  // Acrylic is 1 inch thick
    Double_t fBarrelPMTRadius = 1049.6;
    Double_t fBarrelPMTHeight = 658.8;  // This is HALF the height
    Double_t fTopPMTRadius = -10.0;
    Double_t fBotPMTRadius = -10.0;
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
    std::string startRegion = GetRegion(eventPos);
    if (startRegion == "IV")
      fLightPathType

          // Calculate the refractive indices for the given wavlength
          fIVRefIndex = InterpolateTGraph(fIVRefIndexGraph, wavelength);
    fAcrylicRefIndex = InterpolateTGraph(fAcrylicRefIndexGraph, wavelength);
    fOVRefIndex = InterpolateTGraph(fOVRefIndexGraph, wavelength);

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
      CalculateStraightPath((fEndPos - fStartPos).Unit());
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
      pathResult = CalculateDistancesOV(fStartPos, fEndPos);
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
    fStartPos = startPos;
    fEndPos = endPos;

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

  Bool_t CalculateDistancesOV(const TVector3& eventPos, const TVector3& pmtPos) {}

  void CalculateStraightPath(const TVector3& eventPos, const TVector3& pmtPos) {
    fInitialLightVec = (pmtPos - eventPos).Unit();
  }

  TVector3 GetNormalVector(const TVector3& point) {
    // Check that the point is either within the IV or the acrylic
    if (GetRegion(point) == "OV" || GetRegion(point) == "Past PMT") {
      Log::Die("LightPathCalculator::GetNormalVector: Provided point is outside of the acrylic. Point is (" +
               to_string(point.X()) + ", " + to_string(point.Y()) + ", " + to_string(point.Z()) + ")");
    }

    TVector3 normalVector;
    if (TMath::Abs(point.Z()) < fBarrelPMTHeight) {
      normalVector.SetXYZ(point.X(), point.Y(), 0.0);
      normalVector = normalVector.Unit();
    } else {
      normalVector = point;
      normalVector = normalVector.Unit();
    }
    return normalVector;
  }

  std::string GetRegion(const TVector3& point) {
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

/*std::vector<double>::iterator idxIt = std::lower_bound(valOne.begin(), valOne.end(), wl);
idx = std::distance(valOne.begin(), idxIt);
// Linear interpolation
refIdxs[i] = (wl - valOne[idx]) / (valOne[idx+1] - valOne[idx]) * (valTwo[idx+1] - valTwo[idx]) + valTwo[idx];*/