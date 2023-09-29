#include "../include/SDDCalibration.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include "TROOT.h"
#include "TFile.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TString.h"
#include "TF1.h"
#include "TMinuit.h"
#include "TStyle.h"
#include "TPaveStats.h"
#include "TGraphErrors.h"

#include "XrayLines.h"

std::vector<double> SDDCalibration::GetParams(
    int bus_num, int sdd_num, const std::string &input_par)
{
  std::ifstream input_file(input_par);
  std::vector<double> param_tmp;

  if (!input_file.is_open()) {
      std::cerr << "Error: Unable to open the file: " << input_par << std::endl;
      return param_tmp; // Return an empty vector on error
  }

  std::string line;
  while (std::getline(input_file, line)) {
    int bus, sdd;
    std::stringstream ss(line);
    ss >> bus >> sdd;
    if (bus == bus_num && sdd == sdd_num) {
      double param_val;
      while (ss >> param_val) {
        double val = static_cast<double>(param_val);
        param_tmp.push_back(val);
      }
      return param_tmp;
    }
  }
  input_file.close();
  return {};
}

void SDDCalibration::PrintParams(const std::string &input_par)
{
  for (int bus_num = 1; bus_num < 2; ++bus_num) {
    for (int sdd_num = 1; sdd_num < 65; ++sdd_num) {
      std::vector<double> param_tmp = GetParams(bus_num, sdd_num, input_par);
      std::cout << "BUS: " << bus_num << ", SDD: " << sdd_num << std::endl;
      std::cout << "Parameters: ";
      for (auto &param : param_tmp) {
        std::cout << param << " ";
      }
      std::cout << std::endl;
    }
  }
}