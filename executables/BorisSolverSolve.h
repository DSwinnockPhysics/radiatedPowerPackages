#ifndef BORISSOLVERSOLVE_H
#define BORISSOLVERSOLVE_H

#include "ElectronDynamics/BorisSolver.h"
#include "ElectronDynamics/MyCustomBorisSolver.h"
#include "TVector3.h"

namespace BorisSolverSolve
{
    void SolveBoris(rad::BorisSolver solver, double maxTime, double timeStep, TVector3 x0, TVector3 v0, int skipStepsInOutput, double pitchAngle, std::string outputFileName);
    std::vector<std::vector<double>> SolveBorisReturnVoltageAndTime(rad::BorisSolver solver, double maxTime, double timeStep, TVector3 x0, TVector3 v0, int skipStepsInOutput, double pitchAngle);
    std::vector<std::vector<double>> SolveBorisReturnVoltageAndTime(rad::CustomBorisSolver solver, double maxTime, double timeStep, TVector3 x0, TVector3 v0, int skipStepsInOutput, double pitchAngle);
    std::vector<std::vector<double>> SolveBorisReturnVoltageAndTimeWithCulling(rad::CustomBorisSolver solver, double maxTime, double timeStep, TVector3 x0, TVector3 v0, int skipStepsInOutput, double pitchAngle, TVector3 cullVector);


}


#endif