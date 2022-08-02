#include "ElectronDynamics/BorisSolver.h"

#include "ElectronDynamics/MyCustomBorisSolver.h"

#include "ElectronDynamics/QTNMFields.h"
#include "BasicFunctions/Constants.h"

#include "Antennas/HalfWaveDipole.h"
#include "SignalProcessing/Signal.h"
#include "SignalProcessing/InducedVoltage.h"

#include "Eigen/Dense"
#include "SingleElectronRadiation/SingleElectronRadiation.h"
#include "BorisSolverSolve.h"

// STL includes
#include <unistd.h>
#include <getopt.h>
#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <ctime>
#include <tuple>
#include <filesystem>

#include <iomanip>

// ROOT includes
#include "TFile.h"
#include "TTree.h"
#include "TRandom3.h"
#include "TMath.h"
#include "TVector3.h"
#include "TSystem.h"

// My includes for measuring time
#include <chrono>




using namespace rad;
inline constexpr double pi{ 3.141592653589793238462643383279502884 };
inline constexpr double c{299792458};
constexpr double eMass{ 9.1093837015e-31 }; // electron mass, kg

void BorisSolverSolve::SolveBoris(BorisSolver solver, double maxTime, double timeStep, TVector3 x0, TVector3 v0, int skipStepsInOutput, double pitchAngle, std::string outputFileName){
      std::cout << "Starting solver:" << std::endl;
      TVector3 posVec{x0};
      TVector3 velVec{v0};
      TVector3 accVec;

      double xPos;
      double yPos;
      double zPos;
      double xVel;
      double yVel;
      double zVel;
      double time;
      double xAcc;
      double yAcc;
      double zAcc;

      bool isTrapped;
      int isTrapped_val;
      int nRecordedSteps{0};

      //Write the results to a text file
      unsigned int precision = 13;
      std::cout << std::fixed;
      std::ofstream electronTrajectory;
      std::ofstream signalVoltageFile;
      electronTrajectory.open("Outputs/" + outputFileName + ".txt");
      signalVoltageFile.open("Outputs/" + outputFileName + "_Voltage.txt");

      int nTimeSteps = maxTime / timeStep;
      for (int iStep = 0; iStep < nTimeSteps; iStep++) {
            std::tuple<TVector3, TVector3> outputStep = solver.advance_step(timeStep, posVec, velVec);
            posVec = std::get<0>(outputStep);
            velVec = std::get<1>(outputStep);
            accVec = solver.acc(posVec, velVec);

            time = double(iStep+1) * timeStep;
            if (std::fmod(time, 5e-6) < timeStep) std::cout<<time<<" seconds generated"<<std::endl; 
            
            xPos = posVec.X();
            yPos = posVec.Y();
            zPos = posVec.Z();
            xVel = velVec.X();
            yVel = velVec.Y();
            zVel = velVec.Z();
            xAcc = accVec.X();
            yAcc = accVec.Y();
            zAcc = accVec.Z();    
            
            // Save the result after some steps have been skipped - means you can get better accuracy on the trajectory but still lower file size
            if (iStep % skipStepsInOutput == 0) {
                  electronTrajectory << std::setprecision(12) << time << "," << xPos << "," <<  yPos << "," <<  zPos << "," <<  xVel << "," <<  yVel << "," <<  zVel << "," << xAcc << "," << yAcc << "," << zAcc << std::endl;

                  double RVel = sqrt( velVec.X()*velVec.X() + velVec.Y()*velVec.Y() );

                  double thisPitchAngle = abs( atan(RVel / velVec.Z()) );
                  if (thisPitchAngle < pitchAngle) pitchAngle = thisPitchAngle;

                  Eigen::Vector3d antennaPosition( 0.02, 0, 0 );
                  Eigen::Vector3d antennaAlignment( 0, 1, 0 );

                  TVector3 currentField( solver.calc_b_field(posVec) );
                  double B = currentField.Mag();
                  double gamma = pow( 1-velVec.Mag2()/(c*c), 0.5 );
                  double initialKEJ = (gamma-1) * eMass * c * c;
                  double frequency( CalculateFrequency(initialKEJ, B) );
                  double wavelength( 299792458 / frequency);

                  Eigen::Vector3d EPosition(xPos, yPos, zPos);
                  Eigen::Vector3d EVelocity(xVel, yVel, zVel);
                  Eigen::Vector3d EAcceleration(xAcc, yAcc, zAcc);

                  Eigen::Vector3d EField = RelFarEField(antennaPosition, EPosition, EVelocity, EAcceleration)  + RelNearEField(antennaPosition, EPosition, EVelocity, EAcceleration);
                  Eigen::Vector3d vectorEffectiveLength = antennaAlignment * HalfWaveDipoleEffectiveLength(wavelength);
                  
                  double signalVoltage = vectorEffectiveLength.dot(EField);

                  signalVoltageFile << std::setprecision(12) << time << "," << signalVoltage << std::endl; 
            }
            nRecordedSteps++;
      }

      //close the text file
      signalVoltageFile.close();
      electronTrajectory.close();
      
}




// Note the antenna properties here are separate from the one above, this could be an issue at some point
std::vector<std::vector<double>> BorisSolverSolve::SolveBorisReturnVoltageAndTime(BorisSolver solver, double maxTime, double timeStep, TVector3 x0, TVector3 v0, int skipStepsInOutput, double pitchAngle){
      std::vector<std::vector<double>> timeAndVoltage{};
      std::vector<double> outputTimes{};
      std::vector<double> outputVoltages{};   

      //std::cout << "Starting solver:" << std::endl;
      TVector3 posVec{x0};
      TVector3 velVec{v0};
      TVector3 accVec;

      double xPos;
      double yPos;
      double zPos;
      double xVel;
      double yVel;
      double zVel;
      double time;
      double xAcc;
      double yAcc;
      double zAcc;

      bool isTrapped;
      int isTrapped_val;
      int nRecordedSteps{0};

      unsigned int precision = 13;
      std::cout << std::fixed;

      int nTimeSteps = maxTime / timeStep;
      for (int iStep = 0; iStep < nTimeSteps; iStep++) {
            std::tuple<TVector3, TVector3> outputStep = solver.advance_step(timeStep, posVec, velVec);
            posVec = std::get<0>(outputStep);
            velVec = std::get<1>(outputStep);
            accVec = solver.acc(posVec, velVec);

            time = double(iStep+1) * timeStep;
            if (std::fmod(time, 5e-6) < timeStep) std::cout<<time<<" seconds generated"<<std::endl; 
            
            xPos = posVec.X();
            yPos = posVec.Y();
            zPos = posVec.Z();
            xVel = velVec.X();
            yVel = velVec.Y();
            zVel = velVec.Z();
            xAcc = accVec.X();
            yAcc = accVec.Y();
            zAcc = accVec.Z(); 


            
            // Save the result after some steps have been skipped - means you can get better accuracy on the trajectory but still lower size
            if (iStep % skipStepsInOutput == 0) {

                  double RVel = sqrt( velVec.X()*velVec.X() + velVec.Y()*velVec.Y() );

                  double thisPitchAngle = abs( atan(RVel / velVec.Z()) );
                  if (thisPitchAngle < pitchAngle) pitchAngle = thisPitchAngle;

                  Eigen::Vector3d antennaPosition( 0.02, 0, 0 );
                  Eigen::Vector3d antennaAlignment( 0, 1, 0 );

                  TVector3 currentField( solver.calc_b_field(posVec) );
                  double B = currentField.Mag();
                  double gamma = pow( 1-velVec.Mag2()/(c*c), 0.5 );
                  double initialKEJ = (gamma-1) * eMass * c * c;
                  double frequency( CalculateFrequency(initialKEJ, B) );
                  double wavelength( 299792458 / frequency);

                  Eigen::Vector3d EPosition(xPos, yPos, zPos);
                  Eigen::Vector3d EVelocity(xVel, yVel, zVel);
                  Eigen::Vector3d EAcceleration(xAcc, yAcc, zAcc);

                  Eigen::Vector3d EField = RelFarEField(antennaPosition, EPosition, EVelocity, EAcceleration)  + RelNearEField(antennaPosition, EPosition, EVelocity, EAcceleration);
                  Eigen::Vector3d vectorEffectiveLength = antennaAlignment * HalfWaveDipoleEffectiveLength(wavelength);
                  
                  double signalVoltage = vectorEffectiveLength.dot(EField);


                  outputTimes.push_back(time);
                  outputVoltages.push_back(signalVoltage);
            }
            nRecordedSteps++;
      }

      std::cout << std::scientific; // return the output type to scientific 
      timeAndVoltage.push_back(outputTimes);
      timeAndVoltage.push_back(outputVoltages);
      return timeAndVoltage;
}


// This version of the above function uses the CustomBorisSolver instead of the normal one for speed
std::vector<std::vector<double>> BorisSolverSolve::SolveBorisReturnVoltageAndTime(CustomBorisSolver solver, double maxTime, double timeStep, TVector3 x0, TVector3 v0, int skipStepsInOutput, double pitchAngle){
      std::vector<std::vector<double>> timeAndVoltage{};
      std::vector<double> outputTimes{};
      std::vector<double> outputVoltages{};   

      //std::cout << "Starting solver:" << std::endl;
      TVector3 posVec{x0};
      TVector3 velVec{v0};
      TVector3 accVec;

      double xPos;
      double yPos;
      double zPos;
      double xVel;
      double yVel;
      double zVel;
      double time;
      double xAcc;
      double yAcc;
      double zAcc;

      bool isTrapped;
      int isTrapped_val;
      int nRecordedSteps{0};

      unsigned int precision = 13;
      std::cout << std::fixed;

      int nTimeSteps = maxTime / timeStep;
      for (int iStep = 0; iStep < nTimeSteps; iStep++) {
            std::tuple<TVector3, TVector3> outputStep = solver.custom_advance_step(timeStep, posVec, velVec);
            posVec = std::get<0>(outputStep);
            velVec = std::get<1>(outputStep);
            accVec = solver.acc(posVec, velVec);

            time = double(iStep+1) * timeStep;
            if (std::fmod(time, 5e-6) < timeStep) std::cout<<time<<" seconds generated"<<std::endl; 
            
            xPos = posVec.X();
            yPos = posVec.Y();
            zPos = posVec.Z();
            xVel = velVec.X();
            yVel = velVec.Y();
            zVel = velVec.Z();
            xAcc = accVec.X();
            yAcc = accVec.Y();
            zAcc = accVec.Z(); 


            
            // Save the result after some steps have been skipped - means you can get better accuracy on the trajectory but still lower size
            if (iStep % skipStepsInOutput == 0) {

                  double RVel = sqrt( velVec.X()*velVec.X() + velVec.Y()*velVec.Y() );

                  double thisPitchAngle = abs( atan(RVel / velVec.Z()) );
                  if (thisPitchAngle < pitchAngle) pitchAngle = thisPitchAngle;

                  Eigen::Vector3d antennaPosition( 0.02, 0, 0 );
                  Eigen::Vector3d antennaAlignment( 0, 1, 0 );

                  TVector3 currentField( solver.calc_b_field(posVec) );
                  double B = currentField.Mag();
                  double gamma = pow( 1-velVec.Mag2()/(c*c), 0.5 );
                  double initialKEJ = (gamma-1) * eMass * c * c;
                  double frequency( CalculateFrequency(initialKEJ, B) );
                  double wavelength( 299792458 / frequency);

                  Eigen::Vector3d EPosition(xPos, yPos, zPos);
                  Eigen::Vector3d EVelocity(xVel, yVel, zVel);
                  Eigen::Vector3d EAcceleration(xAcc, yAcc, zAcc);

                  Eigen::Vector3d EField = RelFarEField(antennaPosition, EPosition, EVelocity, EAcceleration)  + RelNearEField(antennaPosition, EPosition, EVelocity, EAcceleration);
                  Eigen::Vector3d vectorEffectiveLength = antennaAlignment * HalfWaveDipoleEffectiveLength(wavelength);
                  
                  double signalVoltage = vectorEffectiveLength.dot(EField);


                  outputTimes.push_back(time);
                  outputVoltages.push_back(signalVoltage);
            }
            nRecordedSteps++;
      }

      std::cout << std::scientific; // return the output type to scientific 
      timeAndVoltage.push_back(outputTimes);
      timeAndVoltage.push_back(outputVoltages);
      return timeAndVoltage;
}



std::vector<std::vector<double>> BorisSolverSolve::SolveBorisReturnVoltageAndTimeWithCulling(CustomBorisSolver solver, double maxTime, double timeStep, TVector3 x0, TVector3 v0, int skipStepsInOutput, double pitchAngle, TVector3 cullVector){
      std::vector<std::vector<double>> timeAndVoltage{};
      std::vector<double> outputTimes{};
      std::vector<double> outputVoltages{};   

      //std::cout << "Starting solver:" << std::endl;
      TVector3 posVec{x0};
      TVector3 velVec{v0};
      TVector3 accVec;

      double xPos;
      double yPos;
      double zPos;
      double xVel;
      double yVel;
      double zVel;
      double time;
      double xAcc;
      double yAcc;
      double zAcc;

      bool isTrapped;
      int isTrapped_val;
      int nRecordedSteps{0};

      unsigned int precision = 13;
      std::cout << std::fixed;

      int nTimeSteps = maxTime / timeStep;
      for (int iStep = 0; iStep < nTimeSteps; iStep++) {
            std::tuple<TVector3, TVector3> outputStep = solver.custom_advance_step(timeStep, posVec, velVec);
            posVec = std::get<0>(outputStep);
            velVec = std::get<1>(outputStep);
            accVec = solver.acc(posVec, velVec);

            time = double(iStep+1) * timeStep;
            if (std::fmod(time, 5e-6) < timeStep) std::cout<<time<<" seconds generated"<<std::endl; 
            
            xPos = posVec.X();
            yPos = posVec.Y();
            zPos = posVec.Z();
            xVel = velVec.X();
            yVel = velVec.Y();
            zVel = velVec.Z();
            xAcc = accVec.X();
            yAcc = accVec.Y();
            zAcc = accVec.Z(); 

            // Stop tracking the electron if it leaves the specified bounds
            if ((xPos < -cullVector[0]) || (xPos > cullVector[0])) {
                  std::cout << "Electron culled at position: " << xPos << ", " << yPos << ", " << zPos << std::endl;
                  break;
            }
            if ((yPos < -cullVector[1]) || (yPos > cullVector[1])) {
                  std::cout << "Electron culled at position: " << xPos << ", " << yPos << ", " << zPos << std::endl;
                  break;
            }
            if ((zPos < -cullVector[2]) || (zPos > cullVector[2])) {
                  std::cout << "Electron culled at position: " << xPos << ", " << yPos << ", " << zPos << std::endl;
                  break;
            }

            
            // Save the result after some steps have been skipped - means you can get better accuracy on the trajectory but still lower size
            if (iStep % skipStepsInOutput == 0) {

                  double RVel = sqrt( velVec.X()*velVec.X() + velVec.Y()*velVec.Y() );

                  // double thisPitchAngle = abs( atan(RVel / velVec.Z()) );
                  // if (thisPitchAngle < pitchAngle) pitchAngle = thisPitchAngle;

                  Eigen::Vector3d antennaPosition( 0.02, 0, 0 );
                  Eigen::Vector3d antennaAlignment( 0, 1, 0 );

                  TVector3 currentField( solver.calc_b_field(posVec) );
                  double B = currentField.Mag();
                  double gamma = pow( 1-velVec.Mag2()/(c*c), 0.5 );
                  double initialKEJ = (gamma-1) * eMass * c * c;
                  double frequency( CalculateFrequency(initialKEJ, B) );
                  double wavelength( 299792458 / frequency);

                  Eigen::Vector3d EPosition(xPos, yPos, zPos);
                  Eigen::Vector3d EVelocity(xVel, yVel, zVel);
                  Eigen::Vector3d EAcceleration(xAcc, yAcc, zAcc);

                  Eigen::Vector3d EField = RelFarEField(antennaPosition, EPosition, EVelocity, EAcceleration)  + RelNearEField(antennaPosition, EPosition, EVelocity, EAcceleration);
                  Eigen::Vector3d vectorEffectiveLength = antennaAlignment * HalfWaveDipoleEffectiveLength(wavelength);
                  
                  double signalVoltage = vectorEffectiveLength.dot(EField);


                  outputTimes.push_back(time);
                  outputVoltages.push_back(signalVoltage);
            }
            nRecordedSteps++;
      }

      std::cout << std::scientific; // return the output type to scientific 
      timeAndVoltage.push_back(outputTimes);
      timeAndVoltage.push_back(outputVoltages);
      return timeAndVoltage;
}








int TestProduceTrajectory(int argc, char *argv[]){
      // Start time measurement
      auto start = std::chrono::high_resolution_clock::now();

      // set up output location
      int opt;
      std::string outputDir = " ";

      const option long_opts[] = {
      {"outputDir", required_argument, nullptr, 'd'}
      };

      while((opt = getopt_long(argc, argv, ":d:n:r:l:kh", long_opts, nullptr)) != -1) {
            switch(opt) {
            case 'd':
                  outputDir = optarg;
                  break;
            }
      }

      if (outputDir == " ") {
            std::cout<<"Must specify directory file with -d"<<std::endl;
      exit(1);
      }

      //Name the output file here
      std::string outputFileName{ "88deg_Harmonic_18600eV_ZDisplaced15cm_1e-5s" };



      // set up the field
      const double centralField = 1.0; // Tesla
      double inhomAx  = 0.0;
      double inhomRad = 0.0;      
      double trapLength = 0.3; // m  
      const double RCoil = 0.025; // m
      const double trapDepth = 0.0049; // Tesla

      InhomogeneousBackgroundField* bkg  = new InhomogeneousBackgroundField(centralField, inhomAx, trapLength/2, inhomRad, RCoil);
      const double trapFieldOffset = centralField - bkg->evaluate_field_at_point(TVector3(0, 0, trapLength/2)).Mag();
      // This trap field offset just fixes the field back to the requested centralField after the inhomogeneities are added as they may move it a bit

      const double ICoil = 2.0 * (trapDepth + trapFieldOffset) * RCoil / MU0; // Amps
      InhomogeneousBathtubField* bathtubField = new InhomogeneousBathtubField(RCoil, ICoil, trapLength/2, centralField, inhomAx, inhomRad);

      HarmonicField* harmonicField = new HarmonicField(RCoil, ICoil, centralField);

      // Check magnetic field, output to file
      std::ofstream bFieldFile;
      std::ofstream bFieldFile3D;
      bFieldFile.open("Outputs/bFieldOutput.txt");
      bFieldFile3D.open("Outputs/bFieldOutput3D.txt");
      for ( int zStep = 0; zStep < 1001; zStep++ ) {
            TVector3 currentPosition(0.024, 0.0, -0.3 + zStep*0.6/1000);
            TVector3 currentField(harmonicField->evaluate_field_at_point(currentPosition));
            bFieldFile << currentPosition.X() << "," << currentPosition.Y() << "," <<  currentPosition.Z() << "," <<  currentField.X() << "," <<  currentField.Y() << "," <<  currentField.Z() << std::endl;
            for ( int xStep=0; xStep < 11; xStep++ ) {
                  for (int yStep=0; yStep < 11; yStep++ ) {
                        TVector3 currentPosition(-0.040+xStep*0.080/10, -0.040+yStep*0.080/10, -0.3 + zStep*0.6/1000);
                        TVector3 currentField(harmonicField->evaluate_field_at_point(currentPosition));
                        bFieldFile3D << currentPosition.X() << "," << currentPosition.Y() << "," <<  currentPosition.Z() << "," <<  currentField.X() << "," <<  currentField.Y() << "," <<  currentField.Z() << std::endl;
                  }
            }
      }
      bFieldFile.close();



      // Electron dynamics
      const double TElec = 18600; // eV
      const double gamma = TElec * TMath::Qe() / (ME * TMath::C()*TMath::C()) + 1;
      const double betaSq = 1 - 1 / pow(gamma, 2);
      const double initialSpeed = sqrt(betaSq)*TMath::C();
      const double tau = 2 * R_E / (3 * TMath::C());

      BorisSolver solver(harmonicField, -TMath::Qe(), ME, tau);


      const double maxSimTime = 1e-5;     // seconds
      const double timeStepSize = 1e-12; // seconds, bear in mind the function also skips some number of steps in the output
      const double stepsToSkip = 5;

      double pitchAngle = 88;



      TVector3 posVec(0, 0.0022, 0-0.15);
      TVector3 vDirection(1, 0, 1/tan(pitchAngle/180*pi));
      vDirection = vDirection * (1 / vDirection.Mag()) ;

      TVector3 velVec(vDirection * initialSpeed);
      TVector3 eAcc;

      double BMean{0.0};

      BorisSolverSolve::SolveBoris(solver, maxSimTime, timeStepSize, posVec, velVec, stepsToSkip, pitchAngle, outputFileName);

      // Get the end time:
      auto stop = std::chrono::high_resolution_clock::now();
      auto duration = std::chrono::duration_cast<std::chrono::seconds>(stop - start);
      std::cout << "Running time: " << duration.count() << std::endl;

      return 0;
}

