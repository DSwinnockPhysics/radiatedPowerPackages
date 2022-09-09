// BasicFunctions.cxx

#include <cmath>

#include "BasicFunctions/BasicFunctions.h"
#include "BasicFunctions/Constants.h"

#include "TGraph.h"
#include "TAxis.h"
#include "TVector3.h"
#include "TMath.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TRandom3.h"

#include "FFTtools.h"
#include "FFTWComplex.h"

#include <iostream>
#include <vector>

void rad::setGraphAttr(TGraph *gr)
{
  gr->GetXaxis()->SetTitleSize(0.05);
  gr->GetYaxis()->SetTitleSize(0.05);
  gr->GetXaxis()->SetLabelSize(0.05);
  gr->GetYaxis()->SetLabelSize(0.05);
  gr->SetLineWidth(2);
}

void rad::SetHistAttr(TH1 *h)
{
  h->GetXaxis()->SetTitleSize(0.05);
  h->GetYaxis()->SetTitleSize(0.05);
  h->GetXaxis()->SetLabelSize(0.05);
  h->GetYaxis()->SetLabelSize(0.05);
  h->SetLineWidth(2);
}

void rad::SetHistAttr(TH2 *h)
{
  h->GetXaxis()->SetTitleSize(0.05);
  h->GetYaxis()->SetTitleSize(0.05);
  h->GetZaxis()->SetTitleSize(0.05);
  h->GetXaxis()->SetLabelSize(0.05);
  h->GetYaxis()->SetLabelSize(0.05);
  h->GetZaxis()->SetLabelSize(0.05);
}

double rad::CalcAeHertzianDipole(const double wavelength, const ROOT::Math::XYZVector dipoleDir,
				 const ROOT::Math::XYZPoint ePosition,
				 const ROOT::Math::XYZPoint antennaPoint)			 
{
  double Ae = 3 * pow(wavelength, 2) / (8*TMath::Pi());
  const double psi = TMath::ACos(((ePosition - antennaPoint).Unit()).Dot(dipoleDir));
  Ae *= pow(TMath::Sin(psi), 2);
  return Ae;
}

double rad::CalcAlHertzianDipole(const double wavelength, const ROOT::Math::XYZVector dipoleDir,
				 const ROOT::Math::XYZPoint ePosition,
				 const ROOT::Math::XYZPoint antennaPoint)
{
  double Al = wavelength * TMath::Sqrt(3 / (8*TMath::Pi()));
  const double psi = TMath::ACos(((ePosition - antennaPoint).Unit()).Dot(dipoleDir));
  Al *= TMath::Sin(psi);
  return Al;
}

double rad::CalcRetardedTime(const ROOT::Math::XYZPoint fieldPoint, const ROOT::Math::XYZPoint ePosition, const double labTime)
{
  double time = labTime - TMath::Sqrt((ePosition - fieldPoint).Mag2()) / TMath::C();
  return time;
}

double rad::CalcTimeFromRetardedTime(ROOT::Math::XYZPoint fieldPoint, ROOT::Math::XYZPoint ePosition,
				     double tRet)
{
  double time = tRet + TMath::Sqrt((ePosition - fieldPoint).Mag2()) / TMath::C();
  return time;
}

double rad::CalcTimeFromRetardedTime(TVector3 fieldPoint, TVector3 ePosition, double tRet)
{
  double time = tRet + ((ePosition - fieldPoint).Mag() / TMath::C());
  return time;
}

double rad::GetSpeedFromKE(double T, double particleMass)
{
  double gamma = T*TMath::Qe() / (ME*TMath::C()*TMath::C()) + 1;
  double betaSq = 1 - 1 / pow(gamma, 2);
  double speed = sqrt(betaSq)*TMath::C();
  return speed;
}

double rad::GetGyroradius(TVector3 velocity, TVector3 bField, double particleMass)
{
  double gamma{ 1 / sqrt( 1 - velocity.Dot(velocity)/(TMath::C()*TMath::C()) ) };
  TVector3 vPerp{ velocity - (velocity.Dot(bField.Unit())*bField) };
  double rg{ gamma*particleMass*vPerp.Mag() / (TMath::Qe()*bField.Mag()) };
  return rg;
}

// Very similar to the FFTtools implementation but without the scaling of the x axis the MHz
TGraph* rad::MakePowerSpectrumNorm(const TGraph* grWave)
{
  double *oldY = grWave->GetY();
  double *oldX = grWave->GetX();
  double deltaT = oldX[1] - oldX[0];
  int length = grWave->GetN();
  FFTWComplex *theFFT = FFTtools::doFFT(length, oldY);

  int newLength = (length/2) + 1;
  double *newY = new double[newLength];
  double *newX = new double[newLength];

  double deltaF = 1/(deltaT*length);

  double tempF = 0;
  for(int i=0;i<newLength;i++) {
    float power=pow(FFTtools::getAbs(theFFT[i]),2);
    if(i>0 && i<newLength-1) power*=2; //account for symmetry                                  
    power*=deltaT/(length); //For time-integral squared amplitude                               
    power/=deltaF;//Just to normalise bin-widths                                                   
    //Ends up the same as dt^2, need to integrate the power (multiply by df)                      
    //to get a meaningful number out.                                                            
    newX[i]=tempF;
    newY[i]=power;
    tempF+=deltaF;
  }

  TGraph *grPower = new TGraph(newLength,newX,newY);
  setGraphAttr(grPower);
  delete [] theFFT;
  delete [] newY;
  delete [] newX;
  return grPower;
}

TGraph* rad::MakePowerSpectrumPeriodogram(const TGraph* grWave)
{
  double *oldY = grWave->GetY();
  double *oldX = grWave->GetX();
  double deltaT = oldX[1] - oldX[0];
  int length = grWave->GetN();
  FFTWComplex *theFFT = FFTtools::doFFT(length, oldY);
  double lengthDub = (double)length;
  int newLength = (length/2) + 1;
  double *newY = new double[newLength];
  double *newX = new double[newLength];

  double deltaF = 1/(deltaT*length);

  double tempF = 0;
  for(int i=0;i<newLength;i++) {
    float power=pow(FFTtools::getAbs(theFFT[i]),2);
    if(i>0 && i<newLength-1) power*=2; //account for symmetry
    double scale = lengthDub*lengthDub;
    power /= scale;
    newX[i]=tempF;
    newY[i]=power;
    tempF+=deltaF;    
  }

  TGraph *grPower = new TGraph(newLength,newX,newY);
  setGraphAttr(grPower);
  grPower->GetXaxis()->SetTitle("Frequency [Hz]");
  delete [] theFFT;
  delete [] newY;
  delete [] newX;
  return grPower;
}

double rad::IntegratePowerNorm(const TGraph* grFFT, Int_t firstBin, Int_t lastBin)
{
  double integral = FFTtools::integratePower(grFFT, firstBin, lastBin);
  // Multiply by frequency bin width
  double deltaF = grFFT->GetPointX(1) - grFFT->GetPointX(0);
  integral *= deltaF;
  return integral;
}

// Re-implementation of the filter from FFTtools but without the conversion factors
TGraph* rad::BandPassFilter(const TGraph* grWave, const double minFreq, const double maxFreq) {
    double *oldY = grWave->GetY();
    double *oldX = grWave->GetX();
    double deltaT=oldX[1]-oldX[0];
    int length=grWave->GetN();
    FFTWComplex *theFFT = FFTtools::doFFT(length,oldY);

    int newLength=(length/2)+1;
    double deltaF=1/(deltaT*length); //Hz

    double tempF=0;
    for(int i=0;i<newLength;i++) {
      if(tempF<minFreq || tempF>maxFreq) {
	theFFT[i].re=0;
	theFFT[i].im=0;
      }      
      tempF+=deltaF;
    }

    double *filteredVals = FFTtools::doInvFFT(length,theFFT);

    TGraph *grFiltered = new TGraph(length,oldX,filteredVals);
    delete [] theFFT;
    delete [] filteredVals;
    return grFiltered;
}

double rad::RayleighPDF(const double x, const double sigma) {
  double sigmaSq = sigma*sigma;
  double prob = (x / sigmaSq) * exp(-1*x*x / (2*sigmaSq));
  return prob;
}

double rad::RayleighCDF(const double x, const double sigma) {
  double f = 1.0 - exp( -1*x*x / (2*sigma*sigma) );
  return f;
}

double rad::RayleighPDFFunc(double *x, double *par) {
  double xx = x[0];
  double retVal = RayleighPDF(xx, par[0]);
  return retVal;
}

double rad::RayleighCDFFunc(double *x, double *par) {
  double xx = x[0];
  double retVal = RayleighCDF(xx, par[0]);
  return retVal;
}

void rad::AddWhiteNoiseFrequencyDomainPowerNorm(TGraph* grIn, const double Teff, const int seed) {
  gRandom->SetSeed(seed);
  const double sampleRate = 2*grIn->GetPointX(grIn->GetN()-1);
  const double deltaT = 1.0 / sampleRate;
  const double deltaF = grIn->GetPointX(1) - grIn->GetPointX(0);
  const double sigma = TMath::Sqrt( TMath::K() * Teff * sampleRate );
  TF1* f1 = new TF1("f1", RayleighPDFFunc, 0, 4*sigma, 1);
  f1->SetParameter(0, sigma);

  // Now loop through the bins of the graph and add the noise
  for (int i = 0; i < grIn->GetN(); i++) {
    double noise = pow(f1->GetRandom()*sqrt(0.5), 2)*(1/deltaF)*(deltaT);
    grIn->SetPointY(i, grIn->GetPointY(i) + noise);
  }
  
  delete f1;
}

TVector3 rad::calculate_omega(const TVector3 BField, const double charge, const double energy, const double mass) {
  double gamma_m0 = mass + energy * TMath::Qe() / pow(TMath::C(), 2);
  return (charge * BField * (1.0 / gamma_m0));
}

TGraph* rad::DownmixInPhase(TGraph* grInput, const double freq) {
  TGraph* grOut = new TGraph();
  setGraphAttr(grOut);
  for (int i = 0; i < grInput->GetN(); i++) {
    double time = grInput->GetPointX(i);
    grOut->SetPoint(i, time, grInput->GetPointY(i)*TMath::Cos(2*TMath::Pi()*freq*time));
  }
  return grOut;
}

TGraph* rad::DownmixQuadrature(TGraph* grInput, const double freq) {
  TGraph* grOut = new TGraph();
  setGraphAttr(grOut);
  for (int i = 0; i < grInput->GetN(); i++) {
    double time = grInput->GetPointX(i);
    grOut->SetPoint(i, time, grInput->GetPointY(i)*TMath::Sin(2*TMath::Pi()*freq*time));
  }
  return grOut;
}

void rad::ScaleGraph(TGraph* grInput, const double scale) {
  for (int i = 0; i < grInput->GetN(); i++) {
    grInput->SetPointY(i, grInput->GetPointY(i) * scale);
  }
}

double rad::CalcCyclotronFreq(const double KE, const double B) {
  double freq = TMath::Qe()*B / (ME + (KE*TMath::Qe()/pow(TMath::C(), 2)));
  return freq / (2*TMath::Pi());
}

TGraph* rad::SumGraphs(std::vector<TGraph*> grInput)
{
  TGraph* grOut = new TGraph();
  setGraphAttr(grOut);
  
  // First of all check that our graphs have the same x (time) spacing
  double testSpacing = grInput[0]->GetPointX(1) - grInput[0]->GetPointX(0);
  for (int iGr = 0; iGr < grInput.size(); iGr++) {
    double thisSpacing = grInput[iGr]->GetPointX(1) - grInput[iGr]->GetPointX(0);
    // Return empty graph if not
    if ((thisSpacing-testSpacing)/testSpacing > 1e-10) {
      std::cout<<"Graphs do not have equivalent time spacing! -- returning empty graph."<<std::endl;
      return grOut;
    }
  }

  // The time series may be of different lengths
  // Need to determine the overlapping times between all the graphs
  double latestStart = -DBL_MAX;
  double earliestEnd = DBL_MAX;
  for (int iGr = 0; iGr < grInput.size(); iGr++) {
    if (grInput[iGr]->GetPointX(0) > latestStart)
      latestStart = grInput[iGr]->GetPointX(0);

    if (grInput[iGr]->GetPointX(grInput[iGr]->GetN()-1) < earliestEnd)
      earliestEnd = grInput[iGr]->GetPointX(grInput[iGr]->GetN()-1);
  }
  
  std::vector<int> startIndices;
  std::vector<int> endIndices;
  // Find the indices where these start and end times are satisfied
  for (int iGr = 0; iGr < grInput.size(); iGr++) {
    // Loop through points from start
    for (int iPnt = 0; iPnt < grInput[iGr]->GetN(); iPnt++) {
      if (grInput[iGr]->GetPointX(iPnt) == latestStart) {
	startIndices.push_back(iPnt);
	break;
      }
    }
    // Now loop through the points from the end to get the end 
    for (int iPnt = grInput[iGr]->GetN()-1; iPnt >= 0; iPnt--) {
      if (grInput[iGr]->GetPointX(iPnt) == earliestEnd) {
	endIndices.push_back(iPnt);
	break;
      }
    } // Loop through points
  }

  // Now sum the graphs between the appropriate ranges
  int nPointsToSum = endIndices[0] - startIndices[0];
  for (int iPnt = 0; iPnt < nPointsToSum; iPnt++) {
    double xVal = grInput[0]->GetPointX(startIndices[0]+iPnt);
    double yVal = 0;
    for (int iGr = 0; iGr < grInput.size(); iGr++) {
      yVal += grInput[iGr]->GetPointY(startIndices[iGr]+iPnt);
    }
    grOut->SetPoint(grOut->GetN(), xVal, yVal);
  }
  
  return grOut;
}

TGraph* rad::SampleWaveform(TGraph* grInput, const double sRate)
{
  TGraph* grOut = new TGraph();
  const double sampleSpacing = 1.0 / sRate;
  double sampleTime = grInput->GetPointX(0);

  for (int i = 0; i < grInput->GetN(); i++) {
    double time = grInput->GetPointX(i);
    if (time < sampleTime) continue;
    else if (i == 0) {
      double calcV = grInput->GetPointY(0);
      grOut->SetPoint(grOut->GetN(), sampleTime, calcV);
      sampleTime += sampleSpacing;
    }
    else {
      // Sample the distribution using linear interpolation
      double calcV = grInput->GetPointY(i-1) + (sampleTime - grInput->GetPointX(i-1)) * (grInput->GetPointY(i) - grInput->GetPointY(i-1)) / (time - grInput->GetPointX(i-1));
      grOut->SetPoint(grOut->GetN(), sampleTime, calcV);
      sampleTime += sampleSpacing;
    }
  }

  return grOut;
}

TH1D* rad::GraphToHistogram(TGraph* grInput)
{
  double binWidth = grInput->GetPointX(1) - grInput->GetPointX(0);
  double firstPoint = grInput->GetPointX(0);
  double lastPoint  = grInput->GetPointX(grInput->GetN()-1);
  TH1D* h = new TH1D("h", "", grInput->GetN(), firstPoint - binWidth/2, lastPoint + binWidth/2);
  SetHistAttr(h);
  for (int i = 0; i < grInput->GetN(); i++) {
    h->SetBinContent(i+1, grInput->GetPointY(i)); 
  }
  
  return h;
}

TGraph* rad::SignalProcessGraph(TGraph* grInput, const double downmixFreq, const double sampleRate)
{
  TGraph* grDM = DownmixInPhase(grInput, downmixFreq);
  TGraph* grS1 = SampleWaveform(grDM, 10*sampleRate);
  delete grDM;
  TGraph* grF  = BandPassFilter(grS1, 0.0, sampleRate/2.0);
  delete grS1;
  TGraph* grS2 = SampleWaveform(grF, sampleRate);
  delete grF;

  return grS2;
}

TGraph* rad::MakeFFTMagGraph(TGraph* grInput) {
  double *oldY = grInput->GetY();
  double *oldX = grInput->GetX();
  double deltaT = oldX[1] - oldX[0];
  int length = grInput->GetN();
  FFTWComplex *theFFT = FFTtools::doFFT(length, oldY);
  double lengthDub = (double)length;
  int newLength = (length/2) + 1;
  double *newY = new double[newLength];
  double *newX = new double[newLength];

  double deltaF = 1/(deltaT*length);

  double tempF = 0;
  for(int i = 0; i < newLength; i++) {
    float mag = FFTtools::getAbs(theFFT[i]);
    newX[i] = tempF;
    newY[i] = mag;
    tempF += deltaF;
  }

  TGraph* grMag = new TGraph(newLength, newX, newY);
  setGraphAttr(grMag);
  
  delete [] theFFT;
  delete [] newY;
  delete [] newX;
  return grMag;
}

double rad::HeavisideFunc(double x)
{
  if (x > 0.0)
    return 1.0;
  else
    return 0.0;
}
