#include "BorisSolverSolve.h"


#include "ElectronDynamics/QTNMFields.h"
#include "ElectronDynamics/QTNMApproxFields.h"
#include "TVector3.h"
#include <chrono> 


// Constants
double eCharge{ 1.602176634e-19};
inline constexpr double pi{ 3.141592653589793238462643383279502884 };
inline constexpr double c{299792458};
constexpr double eMass{ 9.1093837015e-31 }; // electron mass, kg          
                
                
                
using namespace rad;
int main(int argc, char *argv[]){
    // Start time measurement
    auto start = std::chrono::high_resolution_clock::now();

    // set up the field
    const double centralField = 1.0; // Tesla
    double inhomAx  = 0.0;
    double inhomRad = 0.0;      
    double trapLength = 0.3; // m  
    const double RCoil = 0.025; // m
    const double trapDepth = 0.0049; // Tesla
    InhomogeneousBackgroundField* bkg  = new InhomogeneousBackgroundField(centralField, inhomAx, trapLength/2, inhomRad, RCoil);
    const double trapFieldOffset = centralField - bkg->evaluate_field_at_point(TVector3(0, 0, trapLength/2)).Mag();
    const double ICoil = 2.0 * (trapDepth + trapFieldOffset) * RCoil / MU0; // Amps
    HarmonicField* harmonicField = new HarmonicField(RCoil, ICoil, centralField);
    //McDonaldHarmonic* mcDHarmonic = new McDonaldHarmonic(RCoil, ICoil, centralField);

    // Set up antenna
    TVector3 antennaPosition( 0.02, 0, 0 );
    TVector3 antennaAlignment( 0, 1, 0 );

    // Electron dynamics
    const double TElec = 18600; // eV
    const double gamma = TElec * TMath::Qe() / (ME * TMath::C()*TMath::C()) + 1;
    const double betaSq = 1 - 1 / pow(gamma, 2);
    const double initialSpeed = sqrt(betaSq)*TMath::C();
    const double tau = 2 * R_E / (3 * TMath::C());
    double pitchAngle = 88;
    TVector3 posVec(0, 0.0022, 0);
    TVector3 vDirection(1, 0, 1/tan(pitchAngle/180*pi));
    vDirection = vDirection * (1 / vDirection.Mag()) ;
    TVector3 velVec(vDirection * initialSpeed);

    // Track the electrons and somehow combine their signals on the antenna together
    double maxTime{ 1e-6 };
    const double timeStep = 1e-12; 
    int skipStepsInOutput{ 5 };           
    
    //BorisSolver solver(harmonicField, -TMath::Qe(), eMass, tau);
    CustomBorisSolver solver(harmonicField, -TMath::Qe(), eMass, tau);
    std::cout << solver.charge << " " << solver.mass << " " << solver.tau << " " << std::endl;

    std::vector<std::vector<double>> signalTimeAndVoltage;
    signalTimeAndVoltage = BorisSolverSolve::SolveBorisReturnVoltageAndTime(solver, maxTime, timeStep, posVec, velVec, skipStepsInOutput, pitchAngle);
    
    unsigned int precision = 13;
    std::cout << std::fixed;
    std::ofstream signalVoltageFile;
    std::string outputFileName{"BackgroundVoltage_1Electron"};
    signalVoltageFile.open("Outputs/" + outputFileName + ".txt");
    for (int i=0; i<signalTimeAndVoltage[0].size(); ++i){
        signalVoltageFile << std::setprecision(12) << signalTimeAndVoltage[0][i] << "," << signalTimeAndVoltage[1][i] << std::endl;
    }



    // Get the end time:
    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(stop - start);
    std::cout << "Running time: " << duration.count() << std::endl;

}