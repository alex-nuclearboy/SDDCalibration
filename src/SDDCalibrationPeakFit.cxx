#include "../include/SDDCalibration.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <memory>
#include <iomanip>

#include "TROOT.h"
#include "TH1D.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TString.h"

std::vector<int> SDDCalibration::GetPeakPositions(
    int bus_num, int sdd_num, const std::string &file_name) {
  std::ifstream in_file(file_name);
  std::vector<int> peak_pos;

  std::string line;
  while (std::getline(in_file, line)) {
    int bus, sdd;
    std::stringstream ss(line);
    ss >> bus >> sdd;
    if (bus == bus_num && sdd == sdd_num) {
      float peak_val;
      while (ss >> peak_val) {
        int val = static_cast<int>(peak_val);
        peak_pos.push_back(val);
      }
      return peak_pos;
    }
  }
  return {};
}

void SDDCalibration::PreFitPeaks(
    std::string input_file, std::string output_file) 
{
  gROOT->SetBatch(kTRUE); // Disable canvas drawing

  // Open output file for storing the fitted histograms
  TFile* out_img_file = new TFile(("output/precalibration/rootfiles/" 
                                  + output_file + ".root").c_str(), "RECREATE");
  
  // Create output file to store parameters for calibration
  std::ofstream out_par_file;
  out_par_file.open(("output/precalibration/parameters/init_par_" + output_file 
                    + ".txt").c_str()); 

  // Loop over all bus and SDD combinations
  for (int bus = 2; bus < 3; ++bus) {
    for (int sdd = 1; sdd < 21; ++sdd) {
      // Get peak positions for the current bus and SDD
      std::vector<int> peak_positions = GetPeakPositions(bus, sdd, input_file);
      if (peak_positions.empty()) {
        std::cerr << "Error: could not find bus " << bus << " and sdd "
                  << sdd << " in file " << input_file << std::endl;
        continue;
      }

      // Get histogram for the current bus and SDD
      TH1D* hist = GetHistogram(bus, sdd);
      hist->Rebin(rebin_factor);
      hist->SetTitle(Form("bus %d, sdd %d", bus, sdd));      
      hist->GetXaxis()->SetRangeUser(adc_min, adc_max);
      hist->GetYaxis()->SetTitle(Form("counts / %d channels", rebin_factor));

      // Construct the fitting function
      const int num_peaks = peak_positions.size();
      const int num_params = num_peaks * 3 + 3;
      std::string func_str = "";
      for (int i = 0; i < num_peaks; ++i) {
        func_str += Form("[%.0f]*exp(-0.5*((x-[%.0f])/[%.0f])**2) + ", 
                          3.0 * i, 3.0 * i + 1, 3.0 * i + 2);
      }
      func_str += Form("[%.0f] + exp((x - [%.0f])*[%.0f])", 
                        num_peaks * 3.0, num_peaks * 3.0 + 1, num_peaks * 3.0 + 2);
      TF1* fit_func = new TF1("fit_func", func_str.c_str(), adc_min, adc_max);

      // Set initial parameters and parameter limits for the fit function
      for (int i = 0; i < num_peaks; ++i) {        
        double bin_content = hist->GetBinContent(hist->FindBin(peak_positions[i]));
        fit_func->SetParameter(i * 3, bin_content);
        fit_func->SetParLimits(i * 3, 0, 1.5 * bin_content);
        fit_func->SetParameter(i * 3 + 1, peak_positions[i]);
        fit_func->SetParLimits(i * 3 + 1, peak_positions[i] - 5, peak_positions[i] + 5);
        fit_func->SetParameter(i * 3 + 2, 25);
        fit_func->SetParLimits(i * 3 + 2, 5, 40);
        fit_func->SetParName(i * 3, Form("peak_amp_%d", i + 1));
        fit_func->SetParName(i * 3 + 1, Form("peak_pos_%d", i + 1));
        fit_func->SetParName(i * 3 + 2, Form("peak_width_%d", i + 1));
      }
      fit_func->SetParameter(num_peaks * 3, hist->GetBinContent(hist->FindBin(4000)));
      fit_func->SetParameter(num_peaks * 3 + 1, 4000);
      fit_func->SetParameter(num_peaks * 3 + 2, 0.01);
      fit_func->SetParName(num_peaks * 3, "bkg1");
      fit_func->SetParName(num_peaks * 3 + 1, "bkg2");
      fit_func->SetParName(num_peaks * 3 + 2, "bkg3");

      fit_func->SetNpx(1000);
      
      // Perform the fit
      const int fit_result = hist->Fit(fit_func, "R");
      
      // Check if the fit was successful
      if (fit_result != 0) {
        std::cerr << "Fit failed for bus " << bus << " and sdd " << sdd << std::endl;
        continue;
      }

      out_par_file << std::setw(1) << bus << std::setw(4) << sdd;
      for (int i = 0; i < num_peaks; ++i) {
        double peak_pos = fit_func->GetParameter(i * 3 + 1);
        out_par_file << std::setw(10) << peak_pos;
      }
      out_par_file << std::setw(10) << fit_func->GetParameter(num_peaks * 3);
      out_par_file << std::setw(10) << fit_func->GetParameter(num_peaks * 3 + 1);
      out_par_file << std::setw(12) << fit_func->GetParameter(num_peaks * 3 + 2);      
      out_par_file << std::endl;
        
      // Plot the fitted histogram and function
      TCanvas* canvas = new TCanvas(Form("adc_bus%d_sdd%d", bus, sdd), "", 800, 600);
      hist->Draw();
      
      // Create a new function for each Gaussian and the background
      auto gaussFunc = [&](int peak_idx, double* x, double* par) {
        double amplitude = par[peak_idx * 3];
        double mean = par[peak_idx * 3 + 1];
        double sigma = par[peak_idx * 3 + 2];
        return amplitude * exp(-0.5 * pow((x[0] - mean) / sigma, 2));
      };
      
      auto bkgFunc = [&](double* x, double* par) {
        double bkg1 = par[num_peaks * 3];
        double bkg2 = par[num_peaks * 3 + 1];
        double bkg3 = par[num_peaks * 3 + 2];
        return bkg1 + exp((x[0] - bkg2) * bkg3);
      };
      
      // Draw individual Gaussians and background function
      for (int i = 0; i < num_peaks; ++i) {
        TF1* gauss_func = new TF1(
            Form("gauss_func_%d", i),
            [&, i](double* x, double* par) { return gaussFunc(i, x, par); },
        peak_positions[i] - 100, peak_positions[i] + 100, num_params);
        for (int j = 0; j < num_params; ++j) {
          gauss_func->SetParameter(j, fit_func->GetParameter(j));
          gauss_func->FixParameter(j, fit_func->GetParameter(j));
        }
        gauss_func->SetLineColor(kGreen);
        gauss_func->SetLineWidth(1);
        gauss_func->SetNpx(1000);
        gauss_func->Draw("same");
      }
      
      TF1* bkg_func = new TF1(
          "bkg_func",
          [&](double* x, double* par) { return bkgFunc(x, par); },
      adc_min, adc_max, num_params);
      for (int j = 0; j < num_params; ++j) {
        bkg_func->SetParameter(j, fit_func->GetParameter(j));
        bkg_func->FixParameter(j, fit_func->GetParameter(j));
      }
      bkg_func->SetLineColor(kMagenta);
      bkg_func->SetLineWidth(1);
      bkg_func->SetNpx(1000);
      bkg_func->Draw("same");// Save the fitted histogram to the output file

      canvas->Write();
      delete canvas;
      
      // Clean up memory
      for (int i = 0; i < num_peaks; ++i) {
        delete gROOT->FindObject(Form("gauss_func_%d", i));
      }
      delete bkg_func;
    }
  }
  out_par_file.close();
  out_img_file->Close();
  delete out_img_file;
}

