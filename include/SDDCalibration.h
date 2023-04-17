#ifndef SDDCALIBRATION_H
#define SDDCALIBRATION_H

#include <string>
#include <TFile.h>
#include <TH1.h>
#include <TCanvas.h>

class SDDCalibration {
public:
  SDDCalibration(const std::string &filename);
  ~SDDCalibration();

  TH1D    *GetHistogram(int bus, int sdd);
  void    DrawHistogram(int bus, int sdd, const std::string &outputFileName);

private:
  TFile *root_file;

  static const int rebin_factor = 4;
  static const int adc_min = 1500;
  static const int adc_max = 4500;

  // Struct to hold pre-calibration values for different lines
  struct LineCalib {
    double energy;
    double intensity;
    double position;
  };

  void    SearchPeaks(int bus, int sdd, TH1D *&hist);
};

#endif
