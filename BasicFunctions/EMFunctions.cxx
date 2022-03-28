/// EMFunctions.cxx

#include "BasicFunctions/EMFunctions.h"
#include "BasicFunctions/Constants.h"

#include "TMath.h"

// Electric field at the field point, calculated from Lienard-Wiechert potentials
ROOT::Math::XYZVector rad::CalcEField(const ROOT::Math::XYZPoint fieldPoint,
				      const ROOT::Math::XYZPoint ePosition,
				      const ROOT::Math::XYZVector eVelocity,
				      const ROOT::Math::XYZVector eAcceleration)
{  
  const ROOT::Math::XYZVector beta = eVelocity * (1.0 / TMath::C());
  const ROOT::Math::XYZVector betaDot = eAcceleration * (1.0 / TMath::C());
  double premult = TMath::Qe() / (4.0 * EPSILON0 * TMath::Pi());
  const double r = TMath::Sqrt((fieldPoint - ePosition).Mag2());
  const ROOT::Math::XYZVector rHat = (fieldPoint - ePosition).Unit();
  ROOT::Math::XYZVector term1 = (rHat - beta)*(1 - beta.Dot(beta)) * (1.0/( pow(1-beta.Dot(rHat) , 3) * r*r));
  ROOT::Math::XYZVector term2 = rHat.Cross( (rHat - beta).Cross(betaDot) ) * (1.0/( TMath::C()*r*pow(1-beta.Dot(rHat), 3)));
  ROOT::Math::XYZVector field = (term1 + term2) * premult;
  return field;
}

ROOT::Math::XYZVector rad::CalcEField(const TVector3 fieldPoint, const TVector3 ePosition, const TVector3 eVelocity, const TVector3 eAcceleration)
{
  ROOT::Math::XYZPoint fp(fieldPoint.X(), fieldPoint.Y(), fieldPoint.Z());
  ROOT::Math::XYZPoint ep(ePosition.X(), ePosition.Y(), ePosition.Z());
  ROOT::Math::XYZVector ev(eVelocity.X(), eVelocity.Y(), eVelocity.Z());
  ROOT::Math::XYZVector ea(eAcceleration.X(), eAcceleration.Y(), eAcceleration.Z());
  ROOT::Math::XYZVector field = CalcEField(fp, ep, ev, ea);
  return field;
}

// Magnetic field at the field point, calculated from Lienard-Wiechert potentials
ROOT::Math::XYZVector rad::CalcBField(const ROOT::Math::XYZPoint fieldPoint,
				      const ROOT::Math::XYZPoint ePosition,
				      const ROOT::Math::XYZVector eVelocity,
				      const ROOT::Math::XYZVector eAcceleration)
{
  double premult = -1.0 * MU0 * TMath::Qe() / (4.0 * TMath::Pi());
  const ROOT::Math::XYZVector beta = eVelocity * (1.0 / TMath::C());
  const ROOT::Math::XYZVector betaDot = eAcceleration * (1.0 / TMath::C());
  const double r = TMath::Sqrt((fieldPoint - ePosition).Mag2());
  const ROOT::Math::XYZVector rHat = (fieldPoint - ePosition).Unit();
  ROOT::Math::XYZVector term1 = TMath::C() * rHat.Cross(beta) * (1 - beta.Dot(beta)) * (1.0 / ( r*r * pow(1.0 - beta.Dot(rHat), 3) ));
  ROOT::Math::XYZVector term2 = rHat.Cross(betaDot + rHat.Cross(beta.Cross(betaDot))) * (1.0 / (r * pow(1.0 - beta.Dot(rHat), 3) ));
  ROOT::Math::XYZVector field = premult * (term1 + term2);
  return field;
}

ROOT::Math::XYZVector rad::CalcBField(const TVector3 fieldPoint, const TVector3 ePosition,
				      const TVector3 eVelocity, const TVector3 eAcceleration)
{
  ROOT::Math::XYZPoint fp(fieldPoint.X(), fieldPoint.Y(), fieldPoint.Z());
  ROOT::Math::XYZPoint ep(ePosition.X(), ePosition.Y(), ePosition.Z());
  ROOT::Math::XYZVector ev(eVelocity.X(), eVelocity.Y(), eVelocity.Z());
  ROOT::Math::XYZVector ea(eAcceleration.X(), eAcceleration.Y(), eAcceleration.Z());
  ROOT::Math::XYZVector field = CalcBField(fp, ep, ev, ea);
  return field;
}

// Calculate Poynting vector
ROOT::Math::XYZVector rad::CalcPoyntingVec(const ROOT::Math::XYZPoint fieldPoint,
					   const ROOT::Math::XYZPoint ePosition,
					   const ROOT::Math::XYZVector eVelocity,
					   const ROOT::Math::XYZVector eAcceleration)
{
  ROOT::Math::XYZVector EField = CalcEField(fieldPoint, ePosition, eVelocity, eAcceleration);
  ROOT::Math::XYZVector BField = CalcBField(fieldPoint, ePosition, eVelocity, eAcceleration);
  ROOT::Math::XYZVector vec = EField.Cross(BField) * (1.0 / MU0);
  return vec;
}

ROOT::Math::XYZVector rad::CalcPoyntingVec(const TVector3 fieldPoint, const TVector3 ePosition,
					   const TVector3 eVelocity, const TVector3 eAcceleration)
{
  ROOT::Math::XYZPoint fp(fieldPoint.X(), fieldPoint.Y(), fieldPoint.Z());
  ROOT::Math::XYZPoint ep(ePosition.X(), ePosition.Y(), ePosition.Z());
  ROOT::Math::XYZVector ev(eVelocity.X(), eVelocity.Y(), eVelocity.Z());
  ROOT::Math::XYZVector ea(eAcceleration.X(), eAcceleration.Y(), eAcceleration.Z());
  ROOT::Math::XYZVector vec = CalcPoyntingVec(fp, ep, ev, ea);
  return vec;
}

ROOT::Math::XYZVector rad::CalcPoyntingVec(const ROOT::Math::XYZVector EField, const ROOT::Math::XYZVector BField)
{
  ROOT::Math::XYZVector vec = EField.Cross(BField) * (1.0 / MU0);
  return vec;
}
