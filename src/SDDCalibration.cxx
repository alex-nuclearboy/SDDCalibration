#include "../include/SDDCalibration.h"

#include <iostream>

#include "TSpectrum.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TApplication.h"

#include "../include/XrayLines.h"

TApplication *gApplication = nullptr;

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

void SDDCalibration::DrawHistogram(
    int bus, int sdd, const std::string &output_file_name)
{
  TH1D *hist = nullptr;
  SearchPeaks(bus, sdd, hist);

  if (hist) {
    TCanvas *canvas = new TCanvas("canvas", "canvas", 800, 600);
    hist->GetYaxis()->SetTitle(Form("counts / %d channels", rebin_factor));
    hist->Draw();
    canvas->SaveAs(output_file_name.c_str());

    // Keep the canvas open until it is explicitly closed
    canvas->Connect("Closed()", "TApplication", gApplication, "Terminate()");
    canvas->ToggleEditor();
  }
}

void SDDCalibration::SearchPeaks(int bus, int sdd, TH1D *&hist)
{
  hist = GetHistogram(bus, sdd);

  // Rebin histogram
  hist->Rebin(4);

  // Change X-axis range
  hist->GetXaxis()->SetRangeUser(1500, 4500);

  const int npeaks = 14;
  TSpectrum spectrum(npeaks);
  float threshold = 0.04;
  int nfound = 0;
  int tries = 0;

  while (nfound < npeaks && threshold >= 0.0001) {
    nfound = spectrum.Search(hist, 3, "", threshold);
    threshold *= 0.1;
    tries++;
  }

  if (nfound != npeaks) {
    std::cerr << "Error: Could not find " << npeaks 
              << " peaks in histogram h_adc_bus" + std::to_string(bus) \
              + "_sdd" + std::to_string(sdd) << std::endl;
    return;
  }

  std::cout << "Peaks found in " << tries << " attempts." << std::endl;
  const double *peaks = spectrum.GetPositionX();

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
  std::cout << bus << " " << sdd << " " << line_TiKa1.position << " " << line_CuKa1.position
            << " " << line_TiKb1.position << " " << line_CuKb1.position;
  for (int i = 0; i < npeaks; ++i) {
    const double pos = peaks[i];
    if (pos != line_TiKa1.position && pos != line_CuKa1.position && pos != line_TiKb1.position && pos != line_CuKb1.position) {
      std::cout << " " << pos;
    }
  }
  std::cout << std::endl;

}
