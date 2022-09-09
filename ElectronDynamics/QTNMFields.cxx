// QTNMFields.cxx

#include "ElectronDynamics/QTNMFields.h"
#include "BasicFunctions/Constants.h"

#include <cmath>

#include "TVector3.h"
#include "TMath.h"
#include "TSpline.h"
#include "TGraph.h"

#include <boost/math/special_functions/ellint_1.hpp>
#include <boost/math/special_functions/ellint_2.hpp>
#include <boost/math/special_functions/heuman_lambda.hpp>

TVector3 rad::UniformField::evaluate_field_at_point(const TVector3 vec) {
  TVector3 BField(0, 0, fieldStrength);
  return BField;
}

rad::CoilField::CoilField(const double radius, const double current, const double z, const double mu) {
  coilRadius = radius;
  coilCurrent = current;
  coilZ = z;
  coilMu = mu;
}

double rad::CoilField::central_field() {
  return (coilCurrent * coilMu / coilRadius / 2.0);
}

double rad::CoilField::on_axis_field(const double z) {
  double field = coilMu * coilCurrent * pow(coilRadius, 2) / 2.0 / pow(coilRadius*coilRadius + pow(z - coilZ, 2), 1.5);
  return field;
}

TVector3 rad::CoilField::evaluate_field_at_point(const TVector3 vec) {
  double rad = sqrt(vec.X()*vec.X() + vec.Y()*vec.Y()); 

  // Is on axis
  if (rad / coilRadius < 1e-10) {
    TVector3 BField(0, 0, on_axis_field(vec.Z()));
    return BField;
  }

  // z relative to coil position
  double z_rel = vec.Z() - coilZ;

  double b_central = central_field();
  double rad_norm = rad / coilRadius;
  double z_norm = z_rel / coilRadius;
  double alpha = pow(1.0 + rad_norm, 2) + z_norm*z_norm;
  double root_alpha_pi = sqrt(alpha) * TMath::Pi();
  double beta = 4 * rad_norm / alpha;
  double int_k = boost::math::ellint_1(sqrt(beta));
  double int_e = boost::math::ellint_2(sqrt(beta));
  double gamma = alpha - 4 * rad_norm;

  double b_r = b_central * (int_e * ((1.0 + rad_norm*rad_norm + z_norm*z_norm) / gamma) - int_k) / root_alpha_pi *(z_rel / rad);
  double b_z = b_central * (int_e * ((1.0 - rad_norm*rad_norm - z_norm*z_norm) / gamma) + int_k) / root_alpha_pi;

  TVector3 BField(b_r*vec.X()/rad, b_r*vec.Y()/rad, b_z);
  return BField;
}

rad::BathtubField::BathtubField(const double radius, const double current, const double Z1, const double Z2, TVector3 background) {
  coil1 = CoilField(radius, current, Z1, MU0);
  coil2 = CoilField(radius, current, Z2, MU0);
  btBkg = background;
}

TVector3 rad::BathtubField::evaluate_field_at_point(const TVector3 vec) {
  TVector3 field1 = coil1.evaluate_field_at_point(vec);
  TVector3 field2 = coil2.evaluate_field_at_point(vec);
  TVector3 totalField = field1 + field2 + btBkg;
  return totalField;
}

double rad::SolenoidField::on_axis_field(const double z) {
  const double shiftedZ = z - zOff;
  const double xiPlus  = shiftedZ + l/2;
  const double xiMinus = shiftedZ - l/2;
  double field = (mu * n * i / 2) * ( (xiPlus/sqrt(xiPlus*xiPlus + r*r)) - (xiMinus/sqrt(xiMinus*xiMinus + r*r)) );
  return field;
}

TVector3 rad::SolenoidField::evaluate_field_at_point(const TVector3 vec) {
  double rad = sqrt(vec.X()*vec.X() + vec.Y()*vec.Y());

  // Check for on axis case
  if (rad / r < 1e-10) {
    TVector3 BField(0, 0, on_axis_field(vec.Z()));
    return BField;
  }

  double xiPlus  = (vec.Z() - zOff) + l/2;
  double xiMinus = (vec.Z() - zOff) - l/2;

  double premultR = (mu*n*i/TMath::Pi()) * TMath::Sqrt( r/rad );
  double premultZ = (mu*n*i/4);
  double kPlus  = sqrt( 4*rad*r / (xiPlus*xiPlus + pow(rad+r, 2)) );
  double kMinus = sqrt( 4*rad*r / (xiMinus*xiMinus + pow(rad+r, 2)) );

  double int_kPlus  = boost::math::ellint_1(kPlus);
  double int_kMinus = boost::math::ellint_1(kMinus);
  double int_ePlus  = boost::math::ellint_2(kPlus);
  double int_eMinus = boost::math::ellint_2(kMinus);
  
  double Br = ((2-kPlus*kPlus)/(2*kPlus)*int_kPlus - int_ePlus/kPlus) - ((2-kMinus*kMinus)/(2*kMinus)*int_kMinus - int_eMinus/kMinus);
  Br *= premultR;
  
  double phiPlus  = atan( abs(xiPlus / (r - rad)) );
  double phiMinus = atan( abs(xiMinus / (r - rad)) );

  double Bz = ( xiPlus*kPlus*int_kPlus/(TMath::Pi()*sqrt(r*rad)) + (r - rad)*xiPlus*boost::math::heuman_lambda(kPlus,phiPlus)/abs((r-rad)*xiPlus) ) - ( xiMinus*kMinus*int_kMinus/(TMath::Pi()*sqrt(r*rad)) + (r - rad)*xiMinus*boost::math::heuman_lambda(kMinus,phiMinus)/abs((r-rad)*xiMinus) );
  Bz *= premultZ;

  TVector3 BField(Br*vec.X()/rad, Br*vec.Y()/rad, Bz);  
  return BField;
}

TVector3 rad::InhomogeneousBackgroundField::evaluate_field_at_point(const TVector3 vec) {
  const double r = sqrt( vec.X()*vec.X() + vec.Y()*vec.Y() );
  // First calculate the field at the off-axis position r, at z = 0
  const double maxB = ((inhomRad*BCent)/pow(inhomRadRadius, 2)) * r * r + BCent;
  // Use this calculated field to get the axial variation
  const double squareConst = -inhomAx*maxB / pow(inhomAxPosition, 2);
  double field = squareConst * pow(vec.Z(), 2) + maxB;
  TVector3 BField(0, 0, field);
  return BField;
}

rad::InhomogeneousBathtubField::InhomogeneousBathtubField(const double radius, const double current, const double Z, const double centralField, const double fractionalInhomZ, const double fractionalInhomR)
{
  // Generate coil fields at +/- Z
  coil1 = CoilField(radius, current, -1.0*Z, MU0);
  coil2 = CoilField(radius, current, Z, MU0);
  // Generate the background field
  bkgField = InhomogeneousBackgroundField(centralField, fractionalInhomZ, Z, fractionalInhomR, radius);
}

TVector3 rad::InhomogeneousBathtubField::evaluate_field_at_point(const TVector3 vec) {
  TVector3 field1 = coil1.evaluate_field_at_point(vec);
  TVector3 field2 = coil2.evaluate_field_at_point(vec);
  TVector3 field3 = bkgField.evaluate_field_at_point(vec);
  TVector3 totalField = field1 + field2 + field3;
  return totalField;
}

rad::HarmonicField::HarmonicField(const double radius, const double current, const double background)
{
  coil = CoilField(radius, current, 0.0, MU0);
  btBkg = TVector3(0, 0, -background);
}

TVector3 rad::HarmonicField::evaluate_field_at_point(const TVector3 vec)
{
  TVector3 totalField = coil.evaluate_field_at_point(vec) + btBkg;
  return totalField;
}

rad::HTSMagnetUCL::HTSMagnetUCL()
{
  // Create the graphs and splines
  grFieldZ = new TGraph();
  grFieldZ->SetPoint(grFieldZ->GetN(), -0.1, 0.865);
  grFieldZ->SetPoint(grFieldZ->GetN(), -0.074, 1.0);
  grFieldZ->SetPoint(grFieldZ->GetN(), -0.03, 1.04057);
  grFieldZ->SetPoint(grFieldZ->GetN(), -0.025, 1.04069);
  grFieldZ->SetPoint(grFieldZ->GetN(), -0.02, 1.04073);
  grFieldZ->SetPoint(grFieldZ->GetN(), -0.01, 1.04075);
  grFieldZ->SetPoint(grFieldZ->GetN(), 0.0, 1.04075);
  grFieldZ->SetPoint(grFieldZ->GetN(), 0.01, 1.04075);
  grFieldZ->SetPoint(grFieldZ->GetN(), 0.02, 1.04073);
  grFieldZ->SetPoint(grFieldZ->GetN(), 0.025, 1.04069);
  grFieldZ->SetPoint(grFieldZ->GetN(), 0.03, 1.04057);
  grFieldZ->SetPoint(grFieldZ->GetN(), 0.074, 1.0);
  grFieldZ->SetPoint(grFieldZ->GetN(), 0.1, 0.865);
  spFieldZ = new TSpline3("", grFieldZ);

  grFieldR = new TGraph();
  grFieldR->SetPoint(grFieldR->GetN(), 0.000, 1.04075);
  grFieldR->SetPoint(grFieldR->GetN(), 0.005, 1.04075);
  grFieldR->SetPoint(grFieldR->GetN(), 0.010, 1.04075);
  grFieldR->SetPoint(grFieldR->GetN(), 0.015, 1.04075);
  grFieldR->SetPoint(grFieldR->GetN(), 0.020, 1.040751);
  grFieldR->SetPoint(grFieldR->GetN(), 0.027, 1.04076);
  grFieldR->SetPoint(grFieldR->GetN(), 0.030, 1.04077);
  grFieldR->SetPoint(grFieldR->GetN(), 0.040, 1.040835);
  spFieldR = new TSpline3("", grFieldR);  
}

rad::HTSMagnetUCL::~HTSMagnetUCL()
{
  delete grFieldZ;
  delete spFieldZ;
  delete grFieldR;
  delete spFieldR;
}

TVector3 rad::HTSMagnetUCL::evaluate_field_at_point(const TVector3 vec)
{
  double R = sqrt(vec.X()*vec.X() + vec.Y()*vec.Y());
  double scale = spFieldR->Eval(R) / spFieldR->Eval(0.0);
  double fieldZ = spFieldZ->Eval(vec.Z()) * scale;
  return TVector3(0.0, 0.0, fieldZ);
}

rad::HTSMagnetTrap::HTSMagnetTrap(double radius, double current)
{
  coil = CoilField(radius, current, 0.0, MU0);
}

TVector3 rad::HTSMagnetTrap::evaluate_field_at_point(const TVector3 vec)
{
  TVector3 coilField = coil.evaluate_field_at_point(vec);
  TVector3 bkgField  = bkg.evaluate_field_at_point(vec);
  TVector3 totalField = bkgField - coilField;
  return totalField;
}
