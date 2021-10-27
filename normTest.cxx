// normTest.cxx

// STL
#include <iostream>
#include <cmath>

// ROOT includes
#include "TFile.h"
#include "TTree.h"
#include "TVector3.h"
#include "TGraph.h"
#include "TAxis.h"
#include "TMath.h"

// ROOT FFTW Wrapper
#include "FFTtools.h"

// My files
#include "FieldClasses/FieldClasses.h"
#include "BasicFunctions/BasicFunctions.h"

using namespace rad;

int main()
{
  TFile* fin = new TFile("/home/sjones/work/qtnm/trajectories/90DegOnAxis.root", "read");  
  TTree* tree = (TTree*)fin->Get("tree");
  double time;
  double xPos, yPos, zPos;
  double xVel, yVel, zVel;
  double xAcc, yAcc, zAcc;
  tree->SetBranchAddress("time", &time);
  tree->SetBranchAddress("xPos", &xPos);
  tree->SetBranchAddress("yPos", &yPos);
  tree->SetBranchAddress("zPos", &zPos);
  tree->SetBranchAddress("xVel", &xVel);
  tree->SetBranchAddress("yVel", &yVel);
  tree->SetBranchAddress("zVel", &zVel);
  tree->SetBranchAddress("xAcc", &xAcc);
  tree->SetBranchAddress("yAcc", &yAcc);
  tree->SetBranchAddress("zAcc", &zAcc);

  const double maxTime = 1e-6;
  TGraph *grEy = new TGraph();

  TVector3 antennaPoint(0.02, 0, 0);
  
  // Loop through the entries and get the fields at each point
  for (int e = 0; e < tree->GetEntries(); e++) {
    tree->GetEntry(e);
    if (time > maxTime) break;
    TVector3 ePos(xPos, yPos, zPos);
    TVector3 eVel(xVel, yVel, zVel);
    TVector3 eAcc(xAcc, yAcc, zAcc);
    TVector3 EFieldCalc = CalcEField(antennaPoint, ePos, eVel, eAcc);
    TVector3 BFieldCalc = CalcBField(antennaPoint, ePos, eVel, eAcc);
    grEy->SetPoint(grEy->GetN(), time, EFieldCalc.Y());
  }
  TGraph* grEyPower = FFTtools::makePowerSpectrumPeriodogram(grEy);
  
  double vsquared = FFTtools::sumVoltageSquared(grEy, -1, -1)/grEy->GetN();
  double powerSum = FFTtools::sumPower(grEyPower);
  std::cout<<"vsquared, power = "<<vsquared<<", "<<powerSum<<std::endl; 
  
  return 0;
}
