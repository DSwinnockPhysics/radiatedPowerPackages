#ifndef SINGLEELECTRONRADIATION_H
#define SINGLEELECTRONRADIATION_H

#include <iostream>
#include "Eigen/Dense"
#include <math.h>

// ROOT includes
#include "TVector3.h"


namespace rad
{
    double RelativisticTimeDelay(Eigen::Vector3d evaluationPosition,Eigen::Vector3d ePosition);
    Eigen::Vector3d RelFarEField(Eigen::Vector3d evaluationPosition,Eigen::Vector3d ePosition,Eigen::Vector3d eVelocity,Eigen::Vector3d eAcceleration);
    Eigen::Vector3d RelNearEField(Eigen::Vector3d evaluationPosition,Eigen::Vector3d ePosition,Eigen::Vector3d eVelocity,Eigen::Vector3d eAcceleration);
    double PoyntingVectorMagnitude(double electricFieldMagnitude);
    double HalfWaveDipoleEffectiveArea(double wavelength, double dipoleToEmissionAngle);
    double HertzianDipoleEffectiveArea(double wavelength, double dipoleToEmissionAngle);
    double HalfWaveDipoleEffectiveLength(double wavelength);
    double CalculateFrequency( double initialKEJ, double B);

    // TVector3 versions of functions
    TVector3 RelFarEField(TVector3 evaluationPosition, TVector3 ePosition, TVector3 eVelocity, TVector3 eAcceleration);
    TVector3 RelNearEField(TVector3 evaluationPosition, TVector3 ePosition, TVector3 eVelocity, TVector3 eAcceleration);
}

#endif