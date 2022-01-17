// NoiseFunc.cxx

#include <cmath>
#include <time.h>

#include "TMath.h"
#include "TRandom3.h"

#include "SignalProcessing/NoiseFunc.h"

rad::GaussianNoise::GaussianNoise(double T, double R, int setSeed) {
  noiseTemp = T;
  resistance = R;
  numGen = new TRandom3(time(NULL));//(setSeed == 1234) ? new TRandom3(setSeed) : new TRandom3(time(NULL));
}

rad::GaussianNoise::GaussianNoise() {
  noiseTemp = 0;
  resistance = 0;
  sampleFreq = 0;
  sigma = 0;
  numGen = new TRandom3(1234);
}

rad::GaussianNoise::~GaussianNoise() {
  delete numGen;
}

rad::GaussianNoise::GaussianNoise(const GaussianNoise &g1) {
  noiseTemp = g1.noiseTemp;
  sampleFreq = g1.sampleFreq;
  resistance = g1.resistance;
  sigma = g1.sigma;
  numGen = new TRandom3(time(NULL));
}

void rad::GaussianNoise::SetSampleFreq(double fs) {
  sampleFreq = fs;
}

void rad::GaussianNoise::SetResistance(double r) {
  resistance = r;
}

void rad::GaussianNoise::SetSigma() {
  sigma = TMath::Sqrt( TMath::K() * noiseTemp * (sampleFreq/2) );
}

double rad::GaussianNoise::GetNoiseVoltage(bool IsComponent) {
  double premult = IsComponent ? TMath::Sqrt(resistance*0.5*0.5) : TMath::Sqrt(resistance);
  double volt = premult * numGen->Gaus(0.0, sigma);
  return volt;
}



