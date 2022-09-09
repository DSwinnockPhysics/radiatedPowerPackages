// Signal.cxx

#include "SignalProcessing/Signal.h"
#include "SignalProcessing/InducedVoltage.h"
#include "SignalProcessing/LocalOscillator.h"
#include "FieldClasses/FieldClasses.h"
#include "BasicFunctions/BasicFunctions.h"

#include <vector>
#include <iostream>
#include <cassert>
#include <cmath>

#include "TGraph.h"
#include "TH2.h"
#include "TRandom3.h"
#include "TAxis.h"
#include "TMath.h"

#include "FFTtools.h"

rad::Signal::~Signal() {
  delete grVITime;
  delete grVQTime;
}

rad::Signal::Signal()
{
  sampleRate = 0;
  grVITime = 0;
  grVQTime = 0;
}

rad::Signal::Signal(std::vector<FieldPoint> fp, LocalOscillator lo, double srate,
		    std::vector<GaussianNoise> noiseTerms, const bool kUseRetardedTime) {
  assert(fp.size() > 0);
  
  sampleRate = srate;
  
  // Make sure the noise terms are all set up correctly
  for (int iNoise = 0; iNoise < noiseTerms.size(); iNoise++) {
    (noiseTerms.at(iNoise)).SetSampleFreq(sampleRate);
    (noiseTerms.at(iNoise)).SetSigma();
  }

  // The actual output graphs
  grVITime = new TGraph();
  grVQTime = new TGraph();  
  setGraphAttr(grVITime);
  setGraphAttr(grVQTime);
  grVITime->GetXaxis()->SetTitle("Time [s]");
  grVQTime->GetXaxis()->SetTitle("Time [s]");
  
  std::vector<TGraph*> vecVITime;
  std::vector<TGraph*> vecVQTime;
  
  // For each antenna point generate the field and do the signal processing
  for (int point = 0; point < fp.size(); point++) {
    TGraph* grInputVoltageTemp = fp[point].GetAntennaLoadVoltageTimeDomain(kUseRetardedTime, -1, -1);
    
    std::cout<<"Performing the downmixing..."<<std::endl;
    TGraph* grVITimeUnfiltered = DownmixInPhase(grInputVoltageTemp, lo);
    TGraph* grVQTimeUnfiltered = DownmixQuadrature(grInputVoltageTemp, lo);
    delete grInputVoltageTemp;
    
    std::cout<<"First downsampling..."<<std::endl;
    TGraph* grVITimeFirstSample = SampleWaveform(grVITimeUnfiltered, 10*sampleRate);
    TGraph* grVQTimeFirstSample = SampleWaveform(grVQTimeUnfiltered, 10*sampleRate);

    delete grVITimeUnfiltered;
    delete grVQTimeUnfiltered;
    
    // Now we need to filter and then sample these signals
    std::cout<<"Filtering.."<<std::endl;
    TGraph* grVITimeUnsampled = BandPassFilter(grVITimeFirstSample, 0.0, sampleRate/2.0);
    TGraph* grVQTimeUnsampled = BandPassFilter(grVQTimeFirstSample, 0.0, sampleRate/2.0);

    delete grVITimeFirstSample;
    delete grVQTimeFirstSample;
    
    // Now do sampling
    // Use simple linear interpolation for the job
    std::cout<<"Sampling..."<<std::endl;
    TGraph* grVITimeTemp = SampleWaveform(grVITimeUnsampled);
    TGraph* grVQTimeTemp = SampleWaveform(grVQTimeUnsampled);

    delete grVITimeUnsampled;
    delete grVQTimeUnsampled;
    
    std::cout<<"Adding noise..."<<std::endl;
    AddGaussianNoise(grVITimeTemp, noiseTerms);
    AddGaussianNoise(grVQTimeTemp, noiseTerms);

    vecVITime.push_back(grVITimeTemp);
    vecVQTime.push_back(grVQTimeTemp);
  } // Generate individual waveforms

  std::cout<<"Summing fields"<<std::endl;
  const int nSampledTimePoints = vecVITime[0]->GetN();
  for (int t = 0; t < nSampledTimePoints; t++) {
    // Loop over the signals from each individual antenna and sum
    double sumVI = 0;
    double sumVQ = 0;
    for (int field = 0; field < vecVITime.size(); field++) {
      sumVI += vecVITime[field]->GetPointY(t);
      sumVQ += vecVQTime[field]->GetPointY(t);     
    }
    grVITime->SetPoint(t, vecVITime[0]->GetPointX(t), sumVI);
    grVQTime->SetPoint(t, vecVQTime[0]->GetPointX(t), sumVQ);
  }
  
  // Remaining cleanup
  // for (auto p : vecVITime) {
  //    delete p;
  // }
  // vecVITime.clear();
  
  // for (auto p : vecVQTime) {
  //    delete p;
  // }
  // vecVQTime.clear();
  
  std::cout<<"Constructor completed"<<std::endl;
  
}

rad::Signal::Signal(FieldPoint fp, LocalOscillator lo, double srate,
		    std::vector<GaussianNoise> noiseTerms, const bool kUseRetardedTime) {
  sampleRate = srate;
  
  // Make sure the noise terms are all set up correctly
  for (int iNoise = 0; iNoise < noiseTerms.size(); iNoise++) {
    (noiseTerms.at(iNoise)).SetSampleFreq(sampleRate);
    (noiseTerms.at(iNoise)).SetSigma();
  }

  // The actual output graphs
  grVITime = new TGraph();
  grVQTime = new TGraph();  
  setGraphAttr(grVITime);
  setGraphAttr(grVQTime);
  grVITime->GetXaxis()->SetTitle("Time [s]");
  grVQTime->GetXaxis()->SetTitle("Time [s]");

  TGraph* grInputVoltageTemp = fp.GetAntennaLoadVoltageTimeDomain(kUseRetardedTime, -1, -1);
    
  std::cout<<"Performing the downmixing..."<<std::endl;
  TGraph* grVITimeUnfiltered = DownmixInPhase(grInputVoltageTemp, lo);
  TGraph* grVQTimeUnfiltered = DownmixQuadrature(grInputVoltageTemp, lo);
  delete grInputVoltageTemp;
  
  std::cout<<"First downsampling..."<<std::endl;
  TGraph* grVITimeFirstSample = SampleWaveform(grVITimeUnfiltered, 10*sampleRate);
  TGraph* grVQTimeFirstSample = SampleWaveform(grVQTimeUnfiltered, 10*sampleRate);

  delete grVITimeUnfiltered;
  delete grVQTimeUnfiltered;
  
  // Now we need to filter and then sample these signals
  std::cout<<"Filtering.."<<std::endl;
  TGraph* grVITimeUnsampled = BandPassFilter(grVITimeFirstSample, 0.0, sampleRate/2.0);
  TGraph* grVQTimeUnsampled = BandPassFilter(grVQTimeFirstSample, 0.0, sampleRate/2.0);

  delete grVITimeFirstSample;
  delete grVQTimeFirstSample;
  
  // Now do sampling
  // Use simple linear interpolation for the job
  std::cout<<"Sampling..."<<std::endl;
  grVITime = SampleWaveform(grVITimeUnsampled);
  grVQTime = SampleWaveform(grVQTimeUnsampled);

  delete grVITimeUnsampled;
  delete grVQTimeUnsampled; 
  
  std::cout<<"Adding noise..."<<std::endl;
  AddGaussianNoise(grVITime, noiseTerms);
  AddGaussianNoise(grVQTime, noiseTerms);  
}

void rad::Signal::ProcessTimeChunk(InducedVoltage iv, LocalOscillator lo,
				   double thisChunk, double lastChunk,
				   std::vector<GaussianNoise> noiseTerms, 
				   double &firstSampleTime, double &firstSample10Time,
				   bool firstVoltage)
{
  iv.ResetVoltage();
  iv.GenerateVoltage(lastChunk, thisChunk);

  TGraph* grInputVoltageTemp = iv.GetVoltageGraph();

  std::cout<<"Performing the downmixing..."<<std::endl;
  TGraph* grVITimeUnfiltered = DownmixInPhase(grInputVoltageTemp, lo);
  TGraph* grVQTimeUnfiltered = DownmixQuadrature(grInputVoltageTemp, lo);
  delete grInputVoltageTemp;
  
  std::cout<<"First downsampling"<<std::endl;
  TGraph* grVITimeFirstSample = SampleWaveform(grVITimeUnfiltered, 10*sampleRate, firstSample10Time);
  TGraph* grVQTimeFirstSample = SampleWaveform(grVQTimeUnfiltered, 10*sampleRate, firstSample10Time);
  delete grVITimeUnfiltered;
  delete grVQTimeUnfiltered;
  firstSample10Time = grVITimeFirstSample->GetPointX(grVITimeFirstSample->GetN()-1) + 1/(10*sampleRate);

  // Check that the graph contains points
  if (grVITimeFirstSample->GetN() == 0) {
    std::cout<<"We have produced a graph with no sampled points! Finishing this signal chunk."<<std::endl;
    delete grVITimeFirstSample;
    delete grVQTimeFirstSample;    
  }
  else {  
    // Now we need to filter and then sample these signals
    std::cout<<"Filtering.."<<std::endl;
    TGraph* grVITimeUnsampled = BandPassFilter(grVITimeFirstSample, 0.0, sampleRate/2.0);
    TGraph* grVQTimeUnsampled = BandPassFilter(grVQTimeFirstSample, 0.0, sampleRate/2.0);

    delete grVITimeFirstSample;
    delete grVQTimeFirstSample;

    // Now do sampling  
    // Use simple linear interpolation for the job
    std::cout<<"Sampling..."<<std::endl;
    TGraph* grVITimeTemp = SampleWaveform(grVITimeUnsampled, sampleRate, firstSampleTime);
    TGraph* grVQTimeTemp = SampleWaveform(grVQTimeUnsampled, sampleRate, firstSampleTime);
    delete grVITimeUnsampled;
    delete grVQTimeUnsampled;
    firstSampleTime = grVITimeTemp->GetPointX(grVITimeTemp->GetN()-1) + 1/sampleRate;
  
    if (iv.GetLowerAntennaBandwidth() != -DBL_MAX || iv.GetUpperAntennaBandwidth() != DBL_MAX) {
      std::cout<<"Implementing antenna bandwidth..."<<std::endl;
      grVITimeTemp = BandPassFilter(grVITimeTemp, iv.GetLowerAntennaBandwidth()-lo.GetFrequency(), iv.GetUpperAntennaBandwidth()-lo.GetFrequency());
      grVQTimeTemp = BandPassFilter(grVQTimeTemp, iv.GetLowerAntennaBandwidth()-lo.GetFrequency(), iv.GetUpperAntennaBandwidth()-lo.GetFrequency());
    }

    // Now add the information from these temporary graphs to the larger ones
    if (firstVoltage) {
      // We are filling this graph or a new section of the graph for the first time
      for (int i = 0; i < grVITimeTemp->GetN(); i++) {
	grVITime->SetPoint(grVITime->GetN(), grVITimeTemp->GetPointX(i), grVITimeTemp->GetPointY(i));
	grVQTime->SetPoint(grVQTime->GetN(), grVQTimeTemp->GetPointX(i), grVQTimeTemp->GetPointY(i));
      }
    }
    else {
      double tempStartTime = grVITimeTemp->GetPointX(0);
      int startPnt = -1;
      // Loop through main graph to find start points
      for (int iMain = 0; iMain < grVITime->GetN(); iMain++) {
	if (abs(tempStartTime - grVITime->GetPointX(iMain)) < 1e-11) {
	  startPnt = iMain;
	  break;
	}
      }
    
      // Now add the points to the main graph
      for (int i = 0; i < grVITimeTemp->GetN(); i++) {
	// We are adding to existing voltages
	double viTmp = grVITime->GetPointY(startPnt + i);
	double vqTmp = grVQTime->GetPointY(startPnt + i);
	viTmp += grVITimeTemp->GetPointY(i);
	vqTmp += grVQTimeTemp->GetPointY(i);
	grVITime->SetPointY(startPnt + i, viTmp);
	grVQTime->SetPointY(startPnt + i, vqTmp);
      }
    }
      
    delete grVITimeTemp;
    delete grVQTimeTemp;
  } // Sampled graph has non-zero size
}

rad::Signal::Signal(InducedVoltage iv, LocalOscillator lo, double srate,
		    std::vector<GaussianNoise> noiseTerms, double maxTime) {
  sampleRate = srate;

  // Make sure the noise terms are all set up correctly
  for (int iNoise = 0; iNoise < noiseTerms.size(); iNoise++) {
    (noiseTerms.at(iNoise)).SetSampleFreq(sampleRate);
    (noiseTerms.at(iNoise)).SetSigma();
  }

  // The actual output graphs
  grVITime = new TGraph();
  grVQTime = new TGraph();  
  setGraphAttr(grVITime);
  setGraphAttr(grVQTime);
  grVITime->GetYaxis()->SetTitle("V_{I}");
  grVQTime->GetYaxis()->SetTitle("V_{Q}");
  grVITime->GetXaxis()->SetTitle("Time [s]");
  grVQTime->GetXaxis()->SetTitle("Time [s]");

  // Split the signal up into chunks to avoid memory issues
  if (maxTime < 0) maxTime = iv.GetFinalTime();
  
  double lastChunk = 0;
  double thisChunk = lastChunk + iv.GetChunkSize();
  if (thisChunk > maxTime) thisChunk = maxTime;

  double thisSample = 0;
  double this10Sample = 0;
  
  while (thisChunk <= maxTime && thisChunk != lastChunk) {
    ProcessTimeChunk(iv, lo, thisChunk, lastChunk, noiseTerms, thisSample, this10Sample);
    lastChunk = thisChunk;
    thisChunk += iv.GetChunkSize();
    if (thisChunk > maxTime) thisChunk = maxTime;
  }

  if (noiseTerms.size() > 0) {
    std::cout<<"Adding noise..."<<std::endl;
    AddGaussianNoise(grVITime, noiseTerms);
    AddGaussianNoise(grVQTime, noiseTerms);
  }
}

rad::Signal::Signal(std::vector<InducedVoltage> iv, LocalOscillator lo, double srate,
		    std::vector<GaussianNoise> noiseTerms, double maxTime)
{
  sampleRate = srate;

  // Make sure the noise terms are all set up correctly
  for (int iNoise = 0; iNoise < noiseTerms.size(); iNoise++) {
    (noiseTerms.at(iNoise)).SetSampleFreq(sampleRate);
    (noiseTerms.at(iNoise)).SetSigma();
  }

  // The actual output graphs
  grVITime = new TGraph();
  grVQTime = new TGraph();  
  setGraphAttr(grVITime);
  setGraphAttr(grVQTime);
  grVITime->GetYaxis()->SetTitle("V_{I}");
  grVQTime->GetYaxis()->SetTitle("V_{Q}");
  grVITime->GetXaxis()->SetTitle("Time [s]");
  grVQTime->GetXaxis()->SetTitle("Time [s]");

  // Process each voltage one at a time
  for (int iVolt = 0; iVolt < iv.size(); iVolt++) {
  
    if (maxTime < 0) maxTime = iv[iVolt].GetFinalTime();

    // Split the signal up into chunks to avoid memory issues
    double lastChunk = 0;
    double thisChunk = lastChunk + iv[iVolt].GetChunkSize();
    if (thisChunk > maxTime) thisChunk = maxTime;

    double thisSample = 0;
    double this10Sample = 0;

    while (thisChunk <= maxTime && thisChunk != lastChunk) {
      bool firstVoltage = (iVolt == 0);
      ProcessTimeChunk(iv[iVolt], lo, thisChunk, lastChunk, noiseTerms, thisSample, this10Sample, firstVoltage);
      lastChunk = thisChunk;
      thisChunk += iv[iVolt].GetChunkSize();
      if (thisChunk > maxTime) thisChunk = maxTime;
    }    
  } // Loop over InducedVoltage vector

  if (noiseTerms.size() > 0) {
    std::cout<<"Adding noise..."<<std::endl;
    AddGaussianNoise(grVITime, noiseTerms);
    AddGaussianNoise(grVQTime, noiseTerms);
  }
}

rad::Signal::Signal(const Signal &s1) {
  sampleRate = s1.sampleRate;
  grVITime = (TGraph*)s1.grVITime->Clone();
  grVQTime = (TGraph*)s1.grVQTime->Clone();
}

TGraph* rad::Signal::GetVITimeDomain(int firstPoint, int lastPoint) {
  if (firstPoint < 0) firstPoint = 0;
  if (lastPoint < 0 ) lastPoint = grVITime->GetN() - 1;
  
  TGraph* gr = new TGraph();
  setGraphAttr(gr);
  gr->GetXaxis()->SetTitle("Time [s]");
  gr->GetYaxis()->SetTitle("V_{I} [V]");

  for (int i = firstPoint; i <= lastPoint; i++) {
    gr->SetPoint(gr->GetN(), grVITime->GetPointX(i), grVITime->GetPointY(i));
  }
  return gr;
}

TGraph* rad::Signal::GetVQTimeDomain(int firstPoint, int lastPoint) {
  if (firstPoint < 0) firstPoint = 0;
  if (lastPoint < 0 ) lastPoint = grVQTime->GetN() - 1;
  
  TGraph* gr = new TGraph();
  setGraphAttr(gr);
  gr->GetXaxis()->SetTitle("Time [s]");
  gr->GetYaxis()->SetTitle("V_{Q} [V]");

  for (int i = firstPoint; i <= lastPoint; i++) {
    gr->SetPoint(gr->GetN(), grVQTime->GetPointX(i), grVQTime->GetPointY(i));
  }
  return gr;
}

TGraph* rad::Signal::DownmixInPhase(TGraph* grInput, LocalOscillator lo) {
  TGraph* grOut = new TGraph();
  for (int i = 0; i < grInput->GetN(); i++) {
    grOut->SetPoint(i, grInput->GetPointX(i), grInput->GetPointY(i)*lo.GetInPhaseComponent(grInput->GetPointX(i)));
  }
  return grOut;
}

TGraph* rad::Signal::DownmixQuadrature(TGraph* grInput, LocalOscillator lo) {
  TGraph* grOut = new TGraph();
  for (int i = 0; i < grInput->GetN(); i++) {
    grOut->SetPoint(i, grInput->GetPointX(i), grInput->GetPointY(i)*lo.GetQuadratureComponent(grInput->GetPointX(i)));
  }
  return grOut;
}

TGraph* rad::Signal::SampleWaveform(TGraph* grInput) {
  TGraph* grOut = new TGraph();
  double sampleSpacing = 1.0 / sampleRate;
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

TGraph* rad::Signal::SampleWaveform(TGraph* grInput, const double sRate) {
  TGraph* grOut = new TGraph();
  double sampleSpacing = 1.0 / sRate;
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

TGraph* rad::Signal::SampleWaveform(TGraph* grInput, const double sRate, const double firstSampleTime) {
  TGraph* grOut = new TGraph();
  double sampleSpacing = 1.0 / sRate;
  double sampleTime = firstSampleTime;
 
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

void rad::Signal::AddGaussianNoise(TGraph* grInput, std::vector<GaussianNoise> noiseTerms,
				   bool IsComponent) {
  double sampleFreqCalc = 1 / (grInput->GetPointX(1) - grInput->GetPointX(0));
  for (int noise = 0; noise < noiseTerms.size(); noise++) {
    (noiseTerms.at(noise)).SetSampleFreq(sampleFreqCalc);
    (noiseTerms.at(noise)).SetSigma();
  }
  
  for (int i = 0; i < grInput->GetN(); i++) {
    double voltage = grInput->GetPointY(i);
    for (int noise = 0; noise < noiseTerms.size(); noise++) {
      voltage += (noiseTerms.at(noise)).GetNoiseVoltage(IsComponent);  
    }
    grInput->SetPointY(i, voltage);
  }
}

TGraph* rad::Signal::GetVIPowerNorm(const double loadResistance, int firstPoint, int lastPoint) {
  TGraph* grTime = GetVITimeDomain(firstPoint, lastPoint);
  TGraph* grOut = MakePowerSpectrumNorm(grTime);
  delete grTime;
  ScaleGraph(grOut, 1/loadResistance);
  setGraphAttr(grOut);
  grOut->GetXaxis()->SetTitle("Frequency [Hz]");
  grOut->GetYaxis()->SetTitle("#frac{V_{I}^{2}}{R} #times (#Deltat)^{2} [W s^{2}]");
  return grOut;
}

TGraph* rad::Signal::GetVQPowerNorm(const double loadResistance, int firstPoint, int lastPoint) {
  TGraph* grTime = GetVQTimeDomain(firstPoint, lastPoint);
  TGraph* grOut = MakePowerSpectrumNorm(grTime);
  delete grTime;
  ScaleGraph(grOut, 1/loadResistance);
  setGraphAttr(grOut);
  grOut->GetXaxis()->SetTitle("Frequency [Hz]");
  grOut->GetYaxis()->SetTitle("#frac{V_{Q}^{2}}{R} #times (#Deltat)^{2} [W s^{2}]");
  return grOut;
}

TGraph* rad::Signal::GetVIPowerPeriodogram(const double loadResistance, int firstPoint, int lastPoint) {
  TGraph* grTime = GetVITimeDomain(firstPoint, lastPoint);
  TGraph* grOut = MakePowerSpectrumPeriodogram(grTime);
  delete grTime;
  setGraphAttr(grOut);
  grOut->GetXaxis()->SetTitle("Frequency [Hz]");
  grOut->GetYaxis()->SetTitle("Power [W]");
  ScaleGraph(grOut, 1/loadResistance);
  return grOut;
}

TGraph* rad::Signal::GetVQPowerPeriodogram(const double loadResistance, int firstPoint, int lastPoint) {
  TGraph* grTime = GetVQTimeDomain(firstPoint, lastPoint);
  TGraph* grOut = MakePowerSpectrumPeriodogram(grTime);
  delete grTime;
  setGraphAttr(grOut);
  grOut->GetXaxis()->SetTitle("Frequency [Hz]");
  grOut->GetYaxis()->SetTitle("Power [W]");
  ScaleGraph(grOut, 1/loadResistance);
  return grOut;
}

TH2D* rad::Signal::GetVISpectrogram(const double loadResistance, const int NSamplesPerTimeBin) {
  // Get the number of time bins given the specified NSamplesPerTimeBin
  const int nTimeBins = int(floor(double(grVITime->GetN()) / double(NSamplesPerTimeBin)));
  const double firstTime = grVITime->GetPointX(0);
  const double lastTime  = grVITime->GetPointX(nTimeBins * NSamplesPerTimeBin - 1);
  const int nFreqBins = (NSamplesPerTimeBin/2)+1;
  const double deltaF = 1.0 / ((1.0/sampleRate) * NSamplesPerTimeBin);
  const double lastFreq = nFreqBins * deltaF;

  TH2D* h2 = new TH2D("h2", "Spectrogram V_{I}; Time [s]; Frequency [Hz]; Power [W]", nTimeBins, firstTime, lastTime, nFreqBins, 0, lastFreq);
  SetHistAttr(h2);
  // Loop through the time bins and generate the power spectrum at each point
  for (int t = 1; t <= nTimeBins; t++) {
    TGraph* grPower = GetVIPowerPeriodogram(loadResistance, (t-1)*NSamplesPerTimeBin, t*NSamplesPerTimeBin-1);
    // Add this data to the histogram
    for (int i = 0; i < grPower->GetN(); i++) {
      h2->SetBinContent(t, i+1, grPower->GetPointY(i));
    }
    delete grPower;
  }
  
  return h2;
}

TH2D* rad::Signal::GetVQSpectrogram(const double loadResistance, const int NSamplesPerTimeBin) {
  // Get the number of time bins given the specified NSamplesPerTimeBin
  const int nTimeBins = int(floor(double(grVQTime->GetN()) / double(NSamplesPerTimeBin)));
  const double firstTime = grVQTime->GetPointX(0);
  const double lastTime  = grVQTime->GetPointX(nTimeBins * NSamplesPerTimeBin - 1);
  const int nFreqBins = (NSamplesPerTimeBin/2)+1;
  const double deltaF = 1.0 / ((1.0/sampleRate) * NSamplesPerTimeBin);
  const double lastFreq = nFreqBins * deltaF;

  TH2D* h2 = new TH2D("h2", "Spectrogram V_{I}; Time [s]; Frequency [Hz]; Power [W]", nTimeBins, firstTime, lastTime, nFreqBins, 0, lastFreq);
  SetHistAttr(h2);
  // Loop through the time bins and generate the power spectrum at each point
  for (int t = 1; t <= nTimeBins; t++) {
    TGraph* grPower = GetVQPowerPeriodogram(loadResistance, (t-1)*NSamplesPerTimeBin, t*NSamplesPerTimeBin-1);
    // Add this data to the histogram
    for (int i = 0; i < grPower->GetN(); i++) {
      h2->SetBinContent(t, i+1, grPower->GetPointY(i));
    }
    delete grPower;
  }
  
  return h2;
}

TH2D* rad::Signal::GetVISparseSpectrogram(const double loadResistance, const int NSamplesPerTimeBin, const double ThresholdPower) {
  TH2D* hSpec = GetVISpectrogram(loadResistance, NSamplesPerTimeBin);
  TH2D* h2 = new TH2D("h2", "Sparse spectrogram V_{I}; Time [s]; Frequency [Hz]", hSpec->GetNbinsX(), hSpec->GetXaxis()->GetBinLowEdge(1), hSpec->GetXaxis()->GetBinUpEdge(hSpec->GetNbinsX()), hSpec->GetNbinsY(), hSpec->GetYaxis()->GetBinLowEdge(1), hSpec->GetYaxis()->GetBinUpEdge(hSpec->GetNbinsY()));

  for (int x = 1; x <= hSpec->GetNbinsX(); x++) {
    for (int y = 1; y <= hSpec->GetNbinsY(); y++) {
      double binPower = hSpec->GetBinContent(x, y);
      (binPower > ThresholdPower) ? h2->SetBinContent(x, y, 1) : h2->SetBinContent(x, y, 0);  
    }
  }
  delete hSpec;
  
  return h2;
}

TH2D* rad::Signal::GetVQSparseSpectrogram(const double loadResistance, const int NSamplesPerTimeBin, const double ThresholdPower) {
  TH2D* hSpec = GetVQSpectrogram(loadResistance, NSamplesPerTimeBin);
  TH2D* h2 = new TH2D("h2", "Sparse spectrogram V_{Q}; Time [s]; Frequency [Hz]", hSpec->GetNbinsX(), hSpec->GetXaxis()->GetBinLowEdge(1), hSpec->GetXaxis()->GetBinUpEdge(hSpec->GetNbinsX()), hSpec->GetNbinsY(), hSpec->GetYaxis()->GetBinLowEdge(1), hSpec->GetYaxis()->GetBinUpEdge(hSpec->GetNbinsY()));

  for (int x = 1; x <= hSpec->GetNbinsX(); x++) {
    for (int y = 1; y <= hSpec->GetNbinsY(); y++) {
      double binPower = hSpec->GetBinContent(x, y);
      (binPower > ThresholdPower) ? h2->SetBinContent(x, y, 1) : h2->SetBinContent(x, y, 0);  
    }
  }
  delete hSpec;
  
  return h2;
}

TH2D* rad::Signal::GetVISpectrogramNorm(const double loadResistance, const int NSamplesPerTimeBin) {
  // Get the number of time bins given the specified NSamplesPerTimeBin
  const int nTimeBins = int(floor(double(grVITime->GetN()) / double(NSamplesPerTimeBin)));
  const double firstTime = grVITime->GetPointX(0);
  const double lastTime  = grVITime->GetPointX(nTimeBins * NSamplesPerTimeBin - 1);
  const int nFreqBins = (NSamplesPerTimeBin/2)+1;
  const double deltaF = 1.0 / ((1.0/sampleRate) * NSamplesPerTimeBin);
  const double lastFreq = nFreqBins * deltaF;

  TH2D* h2 = new TH2D("h2", "Spectrogram V_{I}; Time [s]; Frequency [Hz]; #frac{V_{I}^{2}}{R} #times (#Delta t)^{2} [W s^{2}]", nTimeBins, firstTime, lastTime, nFreqBins, 0, lastFreq);
  SetHistAttr(h2);
  // Loop through the time bins and generate the power spectrum at each point
  for (int t = 1; t <= nTimeBins; t++) {
    TGraph* grPower = GetVIPowerNorm(loadResistance, (t-1)*NSamplesPerTimeBin, t*NSamplesPerTimeBin-1);
    // Add this data to the histogram
    for (int i = 0; i < grPower->GetN(); i++) {
      h2->SetBinContent(t, i+1, grPower->GetPointY(i));
    }
    delete grPower;
  }
  
  return h2;
}

TH2D* rad::Signal::GetVQSpectrogramNorm(const double loadResistance, const int NSamplesPerTimeBin) {
  // Get the number of time bins given the specified NSamplesPerTimeBin
  const int nTimeBins = int(floor(double(grVQTime->GetN()) / double(NSamplesPerTimeBin)));
  const double firstTime = grVQTime->GetPointX(0);
  const double lastTime  = grVQTime->GetPointX(nTimeBins * NSamplesPerTimeBin - 1);
  const int nFreqBins = (NSamplesPerTimeBin/2)+1;
  const double deltaF = 1.0 / ((1.0/sampleRate) * NSamplesPerTimeBin);
  const double lastFreq = nFreqBins * deltaF;

  TH2D* h2 = new TH2D("h2", "Spectrogram V_{Q}; Time [s]; Frequency [Hz]; #frac{V_{Q}^{2}}{R} #times (#Delta t)^{2} [W s^{2}]", nTimeBins, firstTime, lastTime, nFreqBins, 0, lastFreq);
  SetHistAttr(h2);
  // Loop through the time bins and generate the power spectrum at each point
  for (int t = 1; t <= nTimeBins; t++) {
    TGraph* grPower = GetVQPowerNorm(loadResistance, (t-1)*NSamplesPerTimeBin, t*NSamplesPerTimeBin-1);
    // Add this data to the histogram
    for (int i = 0; i < grPower->GetN(); i++) {
      h2->SetBinContent(t, i+1, grPower->GetPointY(i));
    }
    delete grPower;
  }
  
  return h2;
}

TGraph* rad::Signal::GetDechirpedSignalTimeDomain(const double alpha, int firstPoint, int lastPoint)
{
  TGraph* vi = GetVITimeDomain(firstPoint, lastPoint);
  TGraph* vq = GetVQTimeDomain(firstPoint, lastPoint);

  TGraph* gr = new TGraph();
  setGraphAttr(gr);
  gr->GetXaxis()->SetTitle("Time [s]");
  gr->GetYaxis()->SetTitle("Voltage [V]");

  for (int i = 0; i < vi->GetN(); i++) {
    double time = vi->GetPointX(i);
    double dechirpRe = TMath::Cos( -alpha*time*time/2 );
    double dechirpIm = TMath::Sin( -alpha*time*time/2 );
    double dechirpedSignal = vi->GetPointY(i)*dechirpRe - vq->GetPointY(i)*dechirpIm;
    gr->SetPoint(gr->GetN(), time, dechirpedSignal);
  }
  delete vi;
  delete vq;
  
  return gr;
}

TH2D* rad::Signal::GetDechirpedSpectrogram(const double loadResistance, const int NSamplesPerTimeBin, const double alpha)
{  
  // Get the number of time bins given the specified NSamplesPerTimeBin
  TGraph* grTimeTest = GetVITimeDomain();
  const int nTimeBins = int(floor(double(grTimeTest->GetN()) / double(NSamplesPerTimeBin)));
  const double firstTime = grTimeTest->GetPointX(0);
  const double lastTime  = grTimeTest->GetPointX(nTimeBins * NSamplesPerTimeBin - 1);
  const int nFreqBins = (NSamplesPerTimeBin/2)+1;
  const double deltaF = 1.0 / ((1.0/sampleRate) * NSamplesPerTimeBin);
  const double lastFreq = nFreqBins * deltaF;
  delete grTimeTest;
  
  TH2D* h2 = new TH2D("h2", "Dechirped spectrogram; Time [s]; Frequency [Hz]; Power [W]", nTimeBins, firstTime, lastTime, nFreqBins, 0, lastFreq);
  SetHistAttr(h2);

  // Loop through the time bins and generate the power spectrum at each point
  for (int t = 1; t <= nTimeBins; t++) {
    // First get the dechirped time domain signal
    TGraph* grTime = GetDechirpedSignalTimeDomain(alpha, (t-1)*NSamplesPerTimeBin, t*NSamplesPerTimeBin-1);
    TGraph* grPower = MakePowerSpectrumPeriodogram(grTime);
    delete grTime;
    ScaleGraph(grPower, 1.0 / loadResistance);
    
    // Add this data to the histogram
    for (int i = 0; i < grPower->GetN(); i++) {
      h2->SetBinContent(t, i+1, grPower->GetPointY(i));
    }
    delete grPower;
  }
  
  return h2;
}
