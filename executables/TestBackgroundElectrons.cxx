// STL includes
#include <string>
#include <chrono>
#include <getopt.h>
#include <iostream>
#include "ElectronDynamics/QTNMFields.h"
#include <math.h> 
#include <random>

// My includes
#include "SingleElectronRadiation/SingleElectronRadiation.h"
#include "TBetaGenerator.hpp"
#include "executables/BorisSolverSolve.h"

// ROOT includes
#include "TRandom3.h"
#include "TVector3.h"


// Constants
double eCharge{ 1.602176634e-19};
inline constexpr double pi{ 3.141592653589793238462643383279502884 };
inline constexpr double c{299792458};
constexpr double eMass{ 9.1093837015e-31 }; // electron mass, kg



double GetFrequencyFromKEandB(double KE, double BMagnitude){
    return eCharge*BMagnitude / (eMass+KE/pow(c,2)) / (2*pi);
}

double CalculateGamma( double KEJ ){
    return KEJ/( eMass*pow(c,2)) +1;
}

double VelocityFromGamma(double gamma){
    return pow( pow(c,2)*(1-1/pow(gamma,2)),0.5 );
}

double CalcCyclotronRadius( double perpendicularVelocity, double BMagnitude, double gamma, double frequency ){
    return (eCharge*perpendicularVelocity*BMagnitude) / (gamma * eMass * pow( frequency*2*pi, 2 ));
}

double CalcAcceleration( double frequency, double cyclotronRadius ) {
    return pow(frequency*2*pi, 2) * cyclotronRadius;
}


// build generation functor
// operator calls eactly one 
// distribution function.
// Make different functor for 
// alternative distribution.
class betaGenerator {
    // data parameters
    private:
        bool order_;
        double munu_;
        double mNu_;
        double mixing_;

    public:
        // fix parameter at construction
        betaGenerator(bool o, double ml, double mn, double eta) :
            order_(o),
            munu_(ml),
            mNu_(mn),
            mixing_(eta) {} // with parameter constructor

        double operator() (double x) {
            return TBeta::dGammadE(order_, munu_, mNu_, mixing_, x);
        }
};





using namespace rad;
int main(int argc, char *argv[]){
    // Start time measurement
    auto start = std::chrono::high_resolution_clock::now();

    // Set up the electron energy distribution
    // This is the generator for random numbers
    std::random_device rd;
    std::mt19937 mt(rd());

    // beta decay parameter from G4 macro, for instance
    bool order = true;    // normal order
    double mnu = 1.0e-4;  // 0.1 eV neutrino mass [keV]
    double mN  = 0.0;     // no sterile neutrino
    double eta = 0.0;     // no mixing

    // distribution parameter
    int nw = 100; // nw - number of bins
    double lbound = 0.2; // lower energy bound [keV]
    double ubound = TBeta::endAt(mnu, 1); // max energy
    std::piecewise_linear_distribution<double> ed(nw, lbound, ubound, betaGenerator(order, mnu, mN, eta));



    // set up output location
    int opt;
    std::string outputDir = " ";
    std::string outputFileName = " ";
    const option long_opts[] = {
    {"outputDir", required_argument, nullptr, 'd'},
    {"outputFileName", required_argument, nullptr, 'f'}
    };
    while((opt = getopt_long(argc, argv, ":d:f:n:r:l:kh", long_opts, nullptr)) != -1) {
        switch(opt) {
        case 'd':
                outputDir = optarg;
                break;
        case 'f':
                outputFileName = optarg;
                break;
        }
    }
    if (outputDir == " ") {
        std::cout<<"Must specify directory file with -d"<<std::endl;
        exit(1);
    }
    if (outputFileName == " ") {
        std::cout << "Must specify file name with -f" << std::endl;
        exit(1);
    }



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

    // Set up antenna
    TVector3 antennaPosition( 0.02, 0, 0 );
    TVector3 antennaAlignment( 0, 1, 0 );

    // Set up voltage signal to see general level of voltage output
    double totalVoltage{ 0 };

    // Generate A bunch of random electron positions
    TRandom3* randomGenerator = new TRandom3(0);

    double randomNumber = randomGenerator->Uniform();
    std::cout << randomNumber << std::endl;

    int electronsInFrequencyRange{ 0 };
    int numElectronsToSimulate{ 100 };

    std::vector<std::vector<double>> backgroundTimeAndVoltage{};
    std::vector<double> backgroundOutputTimes{};
    std::vector<double> backgroundOutputVoltages{};


    for ( int electronNumber=0; electronNumber < numElectronsToSimulate; ++electronNumber ) {
        if (electronNumber%1000==0){
            std::cout << "Progress: " << electronNumber << std::endl;
        }
        
        double randomRadius = sqrt(randomGenerator->Uniform()) * 0.02;
        double randomAngle = randomGenerator->Uniform() * 2 * TMath::Pi(); 
        double randomXPos{ randomRadius * cos(randomAngle) };
        double randomYPos{ randomRadius * sin(randomAngle) };
        double randomZPos{ (randomGenerator->Uniform()-0.5) * 0.1 };

        TVector3 randomPosition( randomXPos, randomYPos, randomZPos );

        double maxKEJ{ 18600*eCharge };
        double randomKEJ{ ed(mt) * 1000 * eCharge };  // Random kinetic energy from distribution, given in keV so convert to J

        // Calculate the B field at the generated position
        TVector3 currentField(harmonicField->evaluate_field_at_point(randomPosition));


        double generatedFrequency{ GetFrequencyFromKEandB( randomKEJ, currentField.Mag()) };
        double bandwidthLow{ 26879829000 };
        double bandwidthHigh{ 26880829000 };

        bool calcAllElectronFields{ true };

        // originally only looked at electrons with main frequency in a given range, now just calculate all electrons instead
        if ( (generatedFrequency<bandwidthHigh) && (generatedFrequency>bandwidthLow ) || (calcAllElectronFields==true) ) {
            electronsInFrequencyRange+=1;
            TVector3 velocityDirection{ randomGenerator->Uniform()-0.5, randomGenerator->Uniform()-0.5, randomGenerator->Uniform()-0.5 };
            velocityDirection = velocityDirection *  (1 / velocityDirection.Mag());
            double gamma{ CalculateGamma(randomKEJ) };
            TVector3 randomlyOrientedVelocity{ VelocityFromGamma(gamma) * velocityDirection };
            double cosPitchAngle{ randomlyOrientedVelocity.Dot(currentField) / randomlyOrientedVelocity.Mag() / currentField.Mag() };
            double pitchAngle{ acos(cosPitchAngle) };
            double perpendicularVelocity{ sin(pitchAngle)*randomlyOrientedVelocity.Mag() };
            double currentElectronCyclotronRadius{ CalcCyclotronRadius(perpendicularVelocity, currentField.Mag(), gamma, generatedFrequency ) };
            double currentElectronAcceleration{ CalcAcceleration( generatedFrequency, currentElectronCyclotronRadius ) };
            TVector3 acceleration{ eCharge / eMass / gamma * randomlyOrientedVelocity.Cross(currentField)  };



            // Track the electrons and combine their signals on the antenna together
            double maxTime{ 1e-6 }; // seconds
            const double timeStep = 1e-12; // seconds
            int skipStepsInOutput{ 5 }; 
            const double tau = 2 * R_E / (3 * TMath::C());
            BorisSolver solver(harmonicField, -TMath::Qe(), eMass, tau);
            
            if (electronNumber==0) {
                // For the first electron, just store its voltage in a std::vector and save the times as well
                backgroundTimeAndVoltage = BorisSolverSolve::SolveBorisReturnVoltageAndTime(solver, maxTime, timeStep, randomPosition, randomlyOrientedVelocity, skipStepsInOutput, pitchAngle);
                backgroundOutputTimes = backgroundTimeAndVoltage[0];
                backgroundOutputVoltages = backgroundTimeAndVoltage[1];
            }
            else{
                // For subsequent electrons, add the voltages to get the total voltage
                backgroundTimeAndVoltage = BorisSolverSolve::SolveBorisReturnVoltageAndTime(solver, maxTime, timeStep, randomPosition, randomlyOrientedVelocity, skipStepsInOutput, pitchAngle);
                for (int i=0; i<backgroundOutputTimes.size(); ++i) {
                    backgroundOutputVoltages[i] += backgroundTimeAndVoltage[1][i];
                }
            } 

        }
        else {
            std::cout << "Electron NOT in frequency range" << std::endl;
        }

    }
    std::cout << "Total voltage: " << std::endl;
    std::cout << "Time: " << backgroundOutputTimes[0] << std::endl;
    std::cout << "Time: " << backgroundOutputTimes[1] << std::endl;
    std::cout << "Time: " << backgroundOutputTimes[2] << std::endl;
    std::cout << "Time: " << backgroundOutputTimes[3] << std::endl;
    std::cout << "Time: " << backgroundOutputTimes[4] << std::endl;
    std::cout << "Voltage: " << backgroundOutputVoltages[0] << std::endl;
    std::cout << "Voltage: " << backgroundOutputVoltages[1] << std::endl;
    std::cout << "Voltage: " << backgroundOutputVoltages[2] << std::endl;
    std::cout << "Voltage: " << backgroundOutputVoltages[3] << std::endl;
    std::cout << "Voltage: " << backgroundOutputVoltages[4] << std::endl;
    
    std::ofstream backgroundVoltageFile;
    backgroundVoltageFile.open("Outputs/" + outputFileName + ".txt");
    double RMSBackgroundVoltage{ 0 };
    for (int i=0; i<backgroundOutputTimes.size(); ++i) {
        backgroundOutputVoltages[i] += backgroundTimeAndVoltage[1][i];
        backgroundVoltageFile << std::setprecision(12) << backgroundOutputTimes[i] << "," << backgroundOutputVoltages[i] << std::endl; 
        RMSBackgroundVoltage += pow(backgroundOutputTimes[i],2);
    }  
    backgroundVoltageFile.close();

    RMSBackgroundVoltage /= backgroundOutputTimes.size();
    RMSBackgroundVoltage = pow(RMSBackgroundVoltage, 0.5);

    std::cout << std::scientific;
    std::cout << "Electrons in frequency range: " << electronsInFrequencyRange << std::endl;
    std::cout << "Electrons simulated: " << numElectronsToSimulate << std::endl;
    double fractionOfElectronsInFrequencyRange{ double(electronsInFrequencyRange) / numElectronsToSimulate };
    std::cout << "Fraction of electrons in frequency range: " << fractionOfElectronsInFrequencyRange << std::endl;
    std::cout << "Voltage induced from background electrons: " << totalVoltage << std::endl;
    std::cout << "RMS background electron voltage: " << RMSBackgroundVoltage << std::endl;
    std::cout << "RMS Voltage from background temperature: " << pow(4*1.380649e-23*8*73.2*1e6 , 0.5) <<std::endl;

    // Get the end time:
    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(stop - start);
    std::cout << "Running time: " << duration.count() << std::endl;
}