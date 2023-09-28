#include "../include/SDDCalibration.h"

#include <iostream>
#include <fstream>
#include <iomanip>

#include "TROOT.h"
#include "TH1D.h"
#include "TSpectrum.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TLine.h"
#include "TApplication.h"

#include "../include/XrayLines.h"

TApplication *gApplication = nullptr;

std::vector<double> SDDCalibration::SearchPeaks(
    int bus, int sdd, TH1D *&hist, const std::string &output_file)
{
  hist = GetHistogram(bus, sdd);

  // Rebin histogram
  hist->Rebin(rebin_factor);

  // Change X-axis range
  hist->GetXaxis()->SetRangeUser(1500, 4500);
  hist->GetYaxis()->SetTitle(Form("counts / %d channels", rebin_factor));

  const int npeaks = 10; // Desired number of peaks
  TSpectrum spectrum(40);
  float threshold = 0.04;
  int nfound = 0;
  int tries = 0;

  while (nfound < npeaks && threshold >= 0.0001) {
    nfound = spectrum.Search(hist, 3, "", threshold);
    threshold *= 0.1;
    tries++;
  }

  std::cout << "found " << nfound << " peaks" << std::endl;
  std::cout << "threshold " << threshold << std::endl;
 
  if (nfound < npeaks) {
    std::cerr << "Error: Could not find " << npeaks 
              << " peaks in histogram h_adc_bus" + std::to_string(bus) \
              + "_sdd" + std::to_string(sdd) << std::endl;
    // Return an empty vector to indicate failure
    return std::vector<double>();
  }

  std::cout << "Peaks found in " << tries << " attempts." << std::endl;
  const double *peaks = spectrum.GetPositionX();

  // Store the peak positions in a vector
  std::vector<double> peak_positions;
  for (int i = 0; i < npeaks; ++i) {
    peak_positions.push_back(peaks[i]);
  }

  // Find positions for TiKa and CuKa peaks
  double pos_TiKa = 0.0;
  double pos_CuKa = 0.0;
  double pos_TiKa_dist = 0.0;
  double pos_CuKa_dist = 0.0;
  for (int i = 0; i < nfound; ++i) {
    const double pos = peaks[i];
    const double bin_content = hist->GetBinContent(hist->GetXaxis()->FindBin(pos));
    if (pos > 1750.0 && pos < 2300.0 && bin_content > pos_TiKa_dist) {
      pos_TiKa = pos;
      pos_TiKa_dist = bin_content;
    }
    if (pos > 2800.0 && pos < 3500.0 && bin_content > pos_CuKa_dist) {
      pos_CuKa = pos;
      pos_CuKa_dist = bin_content;
    }
  }

  // Pre-calibration values for different lines
  LineCalib line_TiKa1 = {TiKa1, TiKa1_RI, pos_TiKa};
  LineCalib line_TiKa2 = {TiKa2, TiKa2_RI, 0.0};
  LineCalib line_TiKb1 = {TiKb1, TiKb1_RI, 0.0};
  LineCalib line_CuKa1 = {CuKa1, CuKa1_RI, pos_CuKa};
  LineCalib line_CuKa2 = {CuKa2, CuKa2_RI, 0.0};  
  LineCalib line_CuKb1 = {CuKb1, CuKb1_RI, 0.0};

  const double TiKa = (line_TiKa1.energy * line_TiKa1.intensity + line_TiKa2.energy * line_TiKa2.intensity) / (line_TiKa1.intensity + line_TiKa2.intensity);
  const double CuKa = (line_CuKa1.energy * line_CuKa1.intensity + line_CuKa2.energy * line_CuKa2.intensity) / (line_CuKa1.intensity + line_CuKa2.intensity);
  const double slope = (CuKa - TiKa) / (line_CuKa1.position - line_TiKa1.position);
  const double intercept = -1.0 * line_TiKa1.position * slope + TiKa;
  const double p_TiKb = (line_TiKb1.energy - intercept) / slope;
  const double p_CuKb = (line_CuKb1.energy - intercept) / slope;

  // Find positions for TiKb and CuKb peaks
  double pos_TiKb = 0.0;
  double pos_CuKb = 0.0;
  double pos_TiKb_dist = 1e9;
  double pos_CuKb_dist = 1e9;
  for (int i = 0; i < npeaks; ++i) {
    const double pos = peaks[i];
    const double diff1 = std::abs(pos - p_TiKb);
    const double diff2 = std::abs(pos - p_CuKb);
    if (diff1 < pos_TiKb_dist && pos != line_TiKa1.position && pos != line_CuKa1.position) {
      pos_TiKb = pos;
      pos_TiKb_dist = diff1;
    }
    if (diff2 < pos_CuKb_dist && pos != line_TiKa1.position && pos != line_CuKa1.position) {
      pos_CuKb = pos;
      pos_CuKb_dist = diff2;
    }
  }

  line_TiKb1.position = pos_TiKb;
  line_CuKb1.position = pos_CuKb;

  // Print the peak positions, with the main peaks and closest peaks first
  std::cout << "Peak positions:" << std::endl;
  std::cout << std::setw(1) << bus << std::setw(4) << sdd << std::setw(6) 
            << line_TiKa1.position << std::setw(6) << line_CuKa1.position
            << std::setw(6) << line_TiKb1.position << std::setw(6) 
            << line_CuKb1.position;
  for (int i = 0; i < npeaks; ++i) {
    const double pos = peaks[i];
    if ((pos != line_TiKa1.position) && (pos != line_CuKa1.position) 
        && (pos != line_TiKb1.position) && (pos != line_CuKb1.position)) {
      std::cout << std::setw(6) << pos;
    }
  }
  std::cout << std::endl << std::endl;
  for (int i = 0; i < nfound; ++i) {
    const double pos = peaks[i];
    std::cout << pos << " ";
  }
  std::cout << std::endl;

  // Open the output file
  std::string output_file_name = \
      std::string("output/precalibration/peaks/pf_") + output_file + std::string("_bus") \
      + std::to_string(bus) + std::string("_sdd") + std::to_string(sdd) \
      + std::string(".txt");
  std::ofstream output_txt_file(output_file_name);
  // Write the peak positions to the output file
  output_txt_file << "Desired 10 peaks: " << std::endl;
  output_txt_file << bus << "\t" << sdd << "\t" << line_TiKa1.position 
              << "\t" << line_CuKa1.position << "\t" << line_TiKb1.position 
              << "\t" << line_CuKb1.position;
  for (int i = 0; i < npeaks; ++i) {
    const double pos = peaks[i];
    if (pos != line_TiKa1.position && pos != line_CuKa1.position && pos != line_TiKb1.position && pos != line_CuKb1.position) {
      output_txt_file << " " << pos;
    }
  }
  output_txt_file << std::endl << std::endl;

  output_txt_file << "All " << nfound << " peaks (bus " << bus << ", sdd " 
                  << sdd << "): " << std::endl;
 for (int i = 0; i < nfound; ++i) {
    const double pos = peaks[i];
    output_txt_file << pos << "  ";
  }
  std::cout << std::endl;

  std::cout << "Peak positions saved to " << output_file_name << std::endl;

  return peak_positions;

}

void SDDCalibration::DrawHistogram(
    int bus, int sdd, const std::string &output_file)
{
  TH1D *hist = nullptr;
  std::vector<double> peak_positions = SearchPeaks(bus, sdd, hist, output_file);

  if (hist) {
    std::string output_image_file = \
      std::string("output/precalibration/peaks/plots/pf_") + output_file + std::string("_bus") \
      + std::to_string(bus) + std::string("_sdd") + std::to_string(sdd) \
      + std::string(".png");

    TCanvas *canvas = new TCanvas("canvas", "canvas", 800, 600);
    hist->GetYaxis()->SetTitle(Form("counts / %d channels", rebin_factor));
    hist->Draw();

    // Draw lines on the canvas for each peak found by TSpectrum
    for (const double &pos : peak_positions) {
      const int bin = hist->GetXaxis()->FindBin(pos);
      const double count = hist->GetBinContent(bin);
      TLine *line = new TLine(pos, 0, pos, count);
      line->SetLineColor(kGreen);
      line->SetLineWidth(1);
      line->Draw();
    }

    canvas->SaveAs(output_image_file.c_str());

    // Keep the canvas open until it is explicitly closed
    canvas->Connect("Closed()", "TApplication", gApplication, "Terminate()");
    canvas->ToggleEditor();
  }
}
