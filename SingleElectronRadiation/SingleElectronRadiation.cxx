#include <iostream>
#include "Eigen/Dense"
#include <math.h>
#include "SingleElectronRadiation/SingleElectronRadiation.h"

// ROOT includes
#include "TVector3.h"



constexpr double eCharge{ 1.602176634e-19 }; // Coulombs
constexpr double vacPerm{ 8.8541878128e-12 }; // Farads per metre
constexpr double pi{ 3.141592653589793238462643383279502884 };
constexpr double c{ 299792458 }; // speed of light, m/s
constexpr double eMass{ 9.1093837015e-31 }; // electron mass, kg


double rad::RelativisticTimeDelay(Eigen::Vector3d evaluationPosition,Eigen::Vector3d ePosition) {
    Eigen::Vector3d evaluationElectronSeparation{ (evaluationPosition-ePosition) };
    double separationMagnitude{ evaluationElectronSeparation.norm() };
    double timeDelay{ separationMagnitude/c };
    return timeDelay;
}

Eigen::Vector3d rad::RelFarEField(Eigen::Vector3d evaluationPosition,Eigen::Vector3d ePosition,Eigen::Vector3d eVelocity,Eigen::Vector3d eAcceleration) {
    double premultiplier{ eCharge/(4*pi*vacPerm*c) };
    Eigen::Vector3d evaluationElectronSeparation{ (evaluationPosition-ePosition) };
    double separationMagnitude{ evaluationElectronSeparation.norm() };
    Eigen::Vector3d normalizedEvaluationElectronSeparation{ evaluationElectronSeparation };
    normalizedEvaluationElectronSeparation.normalize();
    premultiplier = premultiplier * 1 / pow((1 - normalizedEvaluationElectronSeparation.dot(eVelocity) / c),3);

    Eigen::Vector3d farFieldPart{ normalizedEvaluationElectronSeparation.cross(
                                ( normalizedEvaluationElectronSeparation-eVelocity/c).cross(eAcceleration/c)
                                ) / separationMagnitude };
    Eigen::Vector3d relFarEField{ premultiplier*farFieldPart };

    return relFarEField;
}


Eigen::Vector3d rad::RelNearEField(Eigen::Vector3d evaluationPosition,Eigen::Vector3d ePosition,Eigen::Vector3d eVelocity,Eigen::Vector3d eAcceleration) {
    double premultiplier{ eCharge/(4*pi*vacPerm) };
    Eigen::Vector3d evaluationElectronSeparation{ (evaluationPosition-ePosition) };
    double separationMagnitude{ evaluationElectronSeparation.norm() };
    Eigen::Vector3d normalizedEvaluationElectronSeparation{ evaluationElectronSeparation };
    normalizedEvaluationElectronSeparation.normalize();
    premultiplier = premultiplier * 1 / pow((1 - normalizedEvaluationElectronSeparation.dot(eVelocity) / c),3);

    Eigen::Vector3d nearFieldPart{ (1-pow(eVelocity.norm()/c,2))*(normalizedEvaluationElectronSeparation-eVelocity/c) / pow(separationMagnitude,2) };
    Eigen::Vector3d relNearEField{ premultiplier*nearFieldPart };

    return relNearEField;
}

double rad::PoyntingVectorMagnitude(double electricFieldMagnitude) {
    double poyntingMagnitude{vacPerm*c*pow(electricFieldMagnitude,2)};
    return poyntingMagnitude;
}

double rad::HalfWaveDipoleEffectiveArea(double wavelength, double dipoleToEmissionAngle) {
    double cosFactor{ cos(cos(dipoleToEmissionAngle)*pi/2) / sin(dipoleToEmissionAngle) };
    double effectiveArea{ 1.64 * pow(cosFactor,2) * pow(wavelength,2) / (4*pi) };
    return effectiveArea;
}

double rad::HertzianDipoleEffectiveArea(double wavelength, double dipoleToEmissionAngle) {
    double effectiveArea{ 3/(8*pi)* pow(wavelength,2) * pow(sin(dipoleToEmissionAngle),2) };
    return effectiveArea;
}

double rad::HalfWaveDipoleEffectiveLength(double wavelength) {
    return wavelength/pi;
}


double rad::CalculateFrequency( double initialKEJ, double B) {
    double initialGamma = initialKEJ / (eMass * pow(c,2)) +1;
    double initialAngularFrequency = eCharge*B/(initialGamma*eMass);
    double initialFrequency = initialAngularFrequency / (2*pi);
    return initialFrequency;
}



// TVector3 versions of functions
TVector3 rad::RelFarEField(TVector3 evaluationPosition, TVector3 ePosition, TVector3 eVelocity, TVector3 eAcceleration) {
    double premultiplier{ eCharge/(4*pi*vacPerm*c) };
    TVector3 evaluationElectronSeparation{ (evaluationPosition-ePosition) };
    double separationMagnitude{ evaluationElectronSeparation.Mag() };
    TVector3 normalizedEvaluationElectronSeparation{ evaluationElectronSeparation * (1/evaluationElectronSeparation.Mag()) };
    premultiplier = premultiplier * 1 / pow((1 - normalizedEvaluationElectronSeparation.Dot(eVelocity) / c),3);

    TVector3 farFieldPart{ normalizedEvaluationElectronSeparation.Cross(
                                ( normalizedEvaluationElectronSeparation-eVelocity*(1/c)).Cross(eAcceleration*(1/c))
                                ) * (1/separationMagnitude) };
    TVector3 relFarEField{ premultiplier*farFieldPart };

    return relFarEField;
}

TVector3 rad::RelNearEField(TVector3 evaluationPosition, TVector3 ePosition, TVector3 eVelocity, TVector3 eAcceleration) {
    double premultiplier{ eCharge/(4*pi*vacPerm) };
    TVector3 evaluationElectronSeparation{ (evaluationPosition-ePosition) };
    double separationMagnitude{ evaluationElectronSeparation.Mag() };
    TVector3 normalizedEvaluationElectronSeparation{ evaluationElectronSeparation * (1/evaluationElectronSeparation.Mag()) };
    premultiplier = premultiplier * 1 / pow((1 - normalizedEvaluationElectronSeparation.Dot(eVelocity) * (1/c) ),3);

    TVector3 nearFieldPart{ (1-pow(eVelocity.Mag()/c,2)) * (normalizedEvaluationElectronSeparation-eVelocity*(1/c)) * (1/pow(separationMagnitude,2)) };
    TVector3 relNearEField{ premultiplier*nearFieldPart };

    return relNearEField;
}

