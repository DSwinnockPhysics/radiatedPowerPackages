/*
  writeTrajectory.cxx

  
*/

#include "ElectronDynamics/BorisSolver.h"
#include "ElectronDynamics/QTNMFields.h"
#include "BasicFunctions/Constants.h"

#include <unistd.h>
#include <iostream>
#include <cmath>
#include <tuple>
#include <string>

#include "TFile.h"
#include "TTree.h"
#include "TMath.h"

using namespace rad;

int main(int argc, char *argv[])
{
  int opt;

  std::string outputFile = " ";
  double simTime = -1;
  double simStepSize = -1;
  char* inputFile;
  double pitchAngle = 90;
  bool energyLoss = false;
  bool hasInputFile = false;
  
  while((opt = getopt(argc, argv, ":o:t:s:i:p:e")) != -1) {
    switch(opt) {
    case 'o':
      outputFile = optarg;
      std::cout<<"Output file is "<<outputFile<<std::endl;
      break;
    case 't':
      simTime = atof(optarg);
      std::cout<<"Simulation time is "<<simTime<<std::endl;
      break;
    case 's':
      simStepSize = atof(optarg);
      std::cout<<"Simulation step size is "<<simStepSize<<std::endl;
      break;
    case 'i':
      inputFile = optarg;
      hasInputFile = true;
      std::cout<<"Input file is "<<inputFile<<std::endl;
      break;
    case 'p':
      pitchAngle = atof(optarg);
      std::cout<<"Pitch angle is "<<pitchAngle<<std::endl;
      break;
    case 'e':
      energyLoss = true;
      std::cout<<"Energy loss turned on"<<std::endl;
      break;
    case ':':
      std::cout<<"Option needs a value"<<std::endl;
      break;
    case '?':
      std::cout<<"Unknown option: "<<optopt<<std::endl;
      break;
    }
  }

  // Check mandatory parameters
  if (outputFile == " ") {
    std::cout<<"Must specify output file with -o"<<std::endl;
    exit(1);
  }
  if (simTime == -1) {
    std::cout<<"Please specify a simulation time (in seconds) with -t"<<std::endl;
    exit(1);
  }
  if ((simStepSize == -1) && !hasInputFile) {
    std::cout<<"Either specify a time step size with -s or an input file with -i"<<std::endl;
    exit(1);
  }
  
  
  const double pitchAngleRad = pitchAngle * TMath::Pi() / 180;
  const double TElec = 18600; // eV
  const double gamma = TElec * TMath::Qe() / (ME * TMath::C()*TMath::C()) + 1;
  const double betaSq = 1 - 1 / pow(gamma, 2);
  const double V0 = sqrt(betaSq) * TMath::C();
  const TVector3 B0(0, 0, 1.0); // Background 1T field in Z direction
  double tau = 0.0;

  if (energyLoss) tau = 2 * R_E / (3 * TMath::C());

  TVector3 X0(0, -0.000460441, 0);
  TVector3 vInitial(V0 * sin(pitchAngleRad), 0, V0 * cos(pitchAngleRad));

  double simStartTime = 0;
  // Get initial state from input style
  if (hasInputFile) {
    TFile* fin = new TFile(inputFile, "READ");
    TTree* intree = (TTree*)fin->Get("tree");
    double startX, startY, startZ;
    double startXVel, startYVel, startZVel;
    intree->SetBranchAddress("time", &simStartTime);
    intree->SetBranchAddress("xPos", &startX);
    intree->SetBranchAddress("yPos", &startY);
    intree->SetBranchAddress("zPos", &startZ);
    intree->SetBranchAddress("xVel", &startXVel);
    intree->SetBranchAddress("yVel", &startYVel);
    intree->SetBranchAddress("zVel", &startZVel);

    intree->GetEntry(0);
    double time0 = simStartTime;
    intree->GetEntry(1);
    double time1 = simStartTime;
    simStepSize = time1 - time0;

    intree->GetEntry(intree->GetEntries()-1);
    X0.SetX(startX);
    X0.SetY(startY);
    X0.SetZ(startZ);
    vInitial.SetX(startXVel);
    vInitial.SetY(startYVel);
    vInitial.SetZ(startZVel);
    
    fin->Close();
    delete fin;
  }

  // Set up the QTNM bathtub trap with coils at +/- 25cm
  double zc1 = -0.25;
  double zc2 = 0.25;
  double Rcoil = 0.03;
  double I = 2.0 * 0.0049 * Rcoil / MU0;
  BathtubField* bathtub = new BathtubField(Rcoil, I, zc1, zc2, B0);
  TVector3 maxVec(0, 0, zc1);
  std::cout<<"Max field perturbation = "<<(bathtub->evaluate_field_at_point(maxVec) - B0).Mag()<<std::endl;

  // Set up the Boris solver
  BorisSolver solver(bathtub, -TMath::Qe(), ME, tau);

  // Open the output ROOT file
  TFile* fout = new TFile(outputFile.data(), "RECREATE");
  TTree* tree = new TTree("tree", "tree");

  double time;
  double xPos, yPos, zPos;
  double xVel, yVel, zVel;
  double xAcc, yAcc, zAcc;

  tree->Branch("time", &time, "time/D");
  tree->Branch("xPos", &xPos, "xPos/D");
  tree->Branch("yPos", &yPos, "yPos/D");
  tree->Branch("zPos", &zPos, "zPos/D");
  tree->Branch("xVel", &xVel, "xVel/D");
  tree->Branch("yVel", &yVel, "yVel/D");
  tree->Branch("zVel", &zVel, "zVel/D");
  tree->Branch("xAcc", &xAcc, "xAcc/D");
  tree->Branch("yAcc", &yAcc, "yAcc/D");
  tree->Branch("zAcc", &zAcc, "zAcc/D");

  int nTimeSteps = simTime / simStepSize;

  std::cout<<"Setting initial state"<<std::endl;
  TVector3 eAcc = solver.acc(X0, vInitial);
  time = simStartTime;
  xPos = X0.X();
  yPos = X0.Y();
  zPos = X0.Z();
  xVel = vInitial.X();
  yVel = vInitial.Y();
  zVel = vInitial.Z();
  xAcc = eAcc.X();
  yAcc = eAcc.Y();
  zAcc = eAcc.Z();

  // If there's no input file then write the initial state
  if (!hasInputFile) tree->Fill();

  TVector3 posVec = X0;
  TVector3 velVec = vInitial;  
  // Loop through the remaining steps and advance the dynamics
  for (int i = 1; i < nTimeSteps; i++) {
    time = simStartTime + double(i) * simStepSize;
    std::tuple<TVector3, TVector3> outputStep = solver.advance_step(simStepSize, posVec, velVec);
    posVec = std::get<0>(outputStep);
    velVec = std::get<1>(outputStep);
    std::cout<<posVec.X()<<", "<<posVec.Y()<<", "<<posVec.Z()<<std::endl;
    eAcc = solver.acc(posVec, velVec);
    
    xPos = posVec.X();
    yPos = posVec.Y();
    zPos = posVec.Z();
    xVel = velVec.X();
    yVel = velVec.Y();
    zVel = velVec.Z();
    xAcc = eAcc.X();
    yAcc = eAcc.Y();
    zAcc = eAcc.Z();

    tree->Fill();
  }

  fout->cd();
  tree->Write();
  
  fout->Close();
  delete fout;
  
  return 0;
}
