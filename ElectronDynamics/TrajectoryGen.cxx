/// TrajectoryGen.cxx

#include "ElectronDynamics/TrajectoryGen.h"
#include "ElectronDynamics/BorisSolver.h"
#include "ElectronDynamics/BaseField.h"

#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TVector3.h"
#include "TMath.h"

#include <tuple>
#include <cmath>

rad::ElectronTrajectoryGen::ElectronTrajectoryGen(TString outputFile, BaseField* field, TVector3 initPos, TVector3 initVel, double simStepSize, double simTime, double initialSimTime, double tau)
{
  outputFilePath = outputFile;
  // Check the file path can be opened in
  TFile* fout = new TFile(outputFilePath, "RECREATE");
  if (!fout) {
    // File path not opened correctly
    std::cout<<"File cannot be created. Exiting..."<<std::endl;
    delete fout;
    exit(1);
  }
  else {
    // File path is all good
    fout->Close();
    delete fout;
  }
  
  solver = BorisSolver(field, -TMath::Qe(), ME, tau);

  // Check that various input values make sense
  if (simStepSize <= 0) {
    std::cout<<"Invalid simulation step size ("<<simStepSize<<"). Exiting..."<<std::endl;
    exit(1);
  }
  if (simTime <= 0) {
    std::cout<<"Invalid simulation time ("<<simTime<<"). Exiting..."<<std::endl;
    exit(1);
  }
  
  stepSize  = simStepSize;
  startTime = initialSimTime;
  simulationTime = simTime;
  initialPosition = initPos;
  initialVelocity = initVel;
}

void rad::ElectronTrajectoryGen::GenerateTraj()
{
  // Open the output file
  TFile* fout = new TFile(outputFilePath, "RECREATE");
  TTree* tree = new TTree("tree", "tree");
  double time;
  double xPos, yPos, zPos;
  double xVel, yVel, zVel;
  double xAcc, yAcc, zAcc;
  tree->Branch("time", &time);
  tree->Branch("xPos", &xPos);
  tree->Branch("yPos", &yPos);
  tree->Branch("zPos", &zPos);
  tree->Branch("xVel", &xVel);
  tree->Branch("yVel", &yVel);
  tree->Branch("zVel", &zVel);
  tree->Branch("xAcc", &xAcc);
  tree->Branch("yAcc", &yAcc);
  tree->Branch("zAcc", &zAcc);

  // Set the initial state
  TVector3 eAcc = solver.acc(initialPosition, initialVelocity);
  time = startTime;
  xPos = initialPosition.X();
  yPos = initialPosition.Y();
  zPos = initialPosition.Z();
  xVel = initialVelocity.X();
  yVel = initialVelocity.Y();
  zVel = initialVelocity.Z();
  xAcc = eAcc.X();
  yAcc = eAcc.Y();
  zAcc = eAcc.Z();
  tree->Fill();

  TVector3 ePos = initialPosition;
  TVector3 eVel = initialVelocity;

  int nTimeSteps = simulationTime / stepSize;
  // Advance through the time steps
  for (int i = 1; i < nTimeSteps; i++) {
    time = startTime + double(i) * stepSize;
    std::tuple<TVector3, TVector3> outputStep = solver.advance_step(stepSize, ePos, eVel);

    if (std::fmod(time, 1e-6) < stepSize) {
      std::cout<<time<<" seconds of trajectory simulated..."<<std::endl;
    }
    
    ePos = std::get<0>(outputStep);
    eVel = std::get<1>(outputStep);
    eAcc = solver.acc(ePos, eVel);

    xPos = ePos.X();
    yPos = ePos.Y();
    zPos = ePos.Z();
    xVel = eVel.X();
    yVel = eVel.Y();
    zVel = eVel.Z();
    xAcc = eAcc.X();
    yAcc = eAcc.Y();
    zAcc = eAcc.Z();

    tree->Fill();
  }
  fout->cd();
  tree->Write();
  
  delete tree;
  fout->Close();
  delete fout;
}
