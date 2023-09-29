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

  void          DrawHistogram(int bus, int sdd, const std::string &output_file);  
  void          PreFitPeaks(std::string input_file, std::string output_file);
  void          Calibrate(
                    std::string &input_file, std::string &input_par, 
                    int &bus_num, int &sdd_num);
  void           PrintParams(const std::string &input_par);

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

  TH1D                    *GetHistogram(int bus, int sdd);
  std::vector<double>     SearchPeaks(
                              int bus, int sdd, TH1D *&hist, 
                              const std::string &output_file);
  std::vector<int>        GetPeakPositions(
                              int bus_num, int sdd_num, 
                              const std::string &file_name);
  std::vector<double>     GetParams(
                              int bus_num, int sdd_num, 
                              const std::string &input_par);
  void                    FitSpectrum(TH1D *hist);
};

#endif
