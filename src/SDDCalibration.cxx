#include "../include/SDDCalibration.h"

#include <iostream>

#include "TROOT.h"
#include "TH1D.h"

SDDCalibration::SDDCalibration(const std::string &filename) : root_file(nullptr)
{
  root_file = TFile::Open(filename.c_str());
  if (!root_file || root_file->IsZombie()) {
    std::cerr << "Error opening ROOT file: " << filename << std::endl;
    exit(EXIT_FAILURE);
  }
}

SDDCalibration::~SDDCalibration()
{
  if (root_file) {
    root_file->Close();
    delete root_file;
  }
}

TH1D *SDDCalibration::GetHistogram(int bus, int sdd)
{
  std::string hist_name = 
      "h_adc_bus" + std::to_string(bus) + "_sdd" + std::to_string(sdd);

  TH1D *hist = dynamic_cast<TH1D *>(root_file->Get(hist_name.c_str()));
  if (!hist) {
    std::cerr << "Error getting histogram: " << hist_name << std::endl;
    exit(EXIT_FAILURE);
  }

  return hist;
}
