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

void fcn(int &npar, double *gin, double &f, double *par, int iflag);
double function(double x, double *par);
double line(double *x, double *par)
{
  return par[0] + par[1] * x[0];
}

struct InitParameters {
    double value;
    double low_bnd;
    double upp_bnd;
    double step_size;
};

int num_par, num_bkg_par, num_def_par, num_peaks, bin_width, ndf;
double *par, *par_err;
double *bin_content, *bin_center; 
double chi_square, chi_square_red;
TH1D *fit_func, *fit_gauss, *fit_gauss_Ka2, *fit_tail, *fit_bkg;
double qv, mv, err_mv, err_qv;

TString status_fit = {"0"};


void SDDCalibration::Calibrate(
  std::string &input_file, std::string &input_par, int &bus_num, int &sdd_num)
{
  // Get the peak positions
  const std::vector<double> param_tmp = GetParams(bus_num, sdd_num, input_par);
  if (param_tmp.empty()) {
    std::cerr << "Error: could not find bus " << bus_num << " and sdd "
              << sdd_num << " in file" << input_par << std::endl;
  }

  // Open input file and retrieve histogram
  TFile *root_file = new TFile(input_file.c_str());
  if (!root_file->IsOpen()) {
    std::cerr << "Error: could not open ROOT file " << input_file << std::endl;
    return;
  }

  // Get the histogram with the given bus and SDD numbers
  const std::string hist_name = "h_adc_bus" + std::to_string(bus_num) \
                                + "_sdd" + std::to_string(sdd_num);
  TH1D *input_hist = dynamic_cast<TH1D *>(root_file->Get(hist_name.c_str()));
  if (!input_hist) {
    std::cerr << "Error getting histogram from ROOT file!" << std::endl;
    return;
  }
  if (input_hist->GetEntries() < 1000) {
    std::cout << "NOT ENOUGH STATISTICS (" << input_hist->GetEntries()
              << "). SKIPPING..." << std::endl;
    return;
  }

  // Rebin and set range of histogram
  input_hist->Rebin(4);
  input_hist->GetXaxis()->SetRangeUser(1500, 4500);

  num_bkg_par = 3;
  num_def_par = 4;
  num_peaks   = 14;
  num_par    = 2 * num_peaks + num_bkg_par + num_def_par;
  init  = new double[num_par];
  par   = new double[num_par];
  par_err  = new double[num_par];

  int idx = 5;
  for (int i = 0; i < param_tmp.size(); ++i) {
    init[idx] = param_tmp[i];
    if (idx < num_par - num_bkg_par - 1) {
      idx = idx + 2;
    }
    else {
      idx++;
    }
  }

  TH1D *hist = (TH1D *) input_hist->Clone("hist");
  ///////////
  fitfunc = new TH1D("fitfunc", "", hist->GetNbinsX(), hist->GetXaxis()->GetXmin(), hist->GetXaxis()->GetXmax());
  fit_gaussKa1 = new TH1D("fit_gaussKa1", "", hist->GetNbinsX(), hist->GetXaxis()->GetXmin(), hist->GetXaxis()->GetXmax());
  fit_gaussKa2 = new TH1D("fit_gaussKa2", "", hist->GetNbinsX(), hist->GetXaxis()->GetXmin(), hist->GetXaxis()->GetXmax());
  fit_tail = new TH1D("fit_tail", "", hist->GetNbinsX(), hist->GetXaxis()->GetXmin(), hist->GetXaxis()->GetXmax());
  fit_bkg = new TH1D("fit_bkg", "", hist->GetNbinsX(), hist->GetXaxis()->GetXmin(), hist->GetXaxis()->GetXmax());

  fit(input_hist);

  TCanvas *ca = new TCanvas("", "", 1900, 1000);
  ca->Divide(2, 2, 0.001, 0.001);
  ca->cd(1);
  hist->GetXaxis()->SetRangeUser(1500, 4500);
  gPad->SetLogy();
  hist->Draw("");
  hist->GetXaxis()->SetTitle("ADC [channel]") ;
  hist->GetYaxis()->SetTitle("counts / 4 channels");

  fitfunc->SetLineColor(2);
  fitfunc->SetLineWidth(1);
  fitfunc->Draw("same");
  fit_gaussKa1->SetLineColor(4);
  fit_gaussKa1->Draw("same");
  fit_gaussKa2->SetLineColor(4);
  fit_gaussKa2->Draw("same");
  fit_tail->SetLineColor(3);
  fit_tail->Draw("same");
  fit_bkg->SetLineColor(6);
  fit_bkg->Draw("same");
  ca->Update();

  int np = 2;
  double p1 = par[num_def_par + 1]; //TiKa1
  double ep1 = par_err[num_def_par + 1];
  double p2 = par[num_def_par + 3]; //CuKa1
  double ep2 = par_err[num_def_par + 3];

  double chn[2] = {p1, p2};
  double echn[2] = {ep1, ep2};
  double E[2] = {TiKa1, CuKa1};

  //calib_par
  TGraphErrors *gr1 = new TGraphErrors(2, chn, E, echn, 0);
  TF1 *fit_line = new TF1("fit_line", InterpolateLine, 0, 15000, 2);
  fit_line->SetParNames("Offset", "Slope");
  gr1->SetTitle("");
  fit_line->SetParameters(1, 1);
  fit_line->SetLineColor(kGreen);
  fit_line->SetLineWidth(1);
  gr1->Fit("fit_line");
  gr1->GetXaxis()->SetTitle("ADC[chn]") ;
  gr1->GetYaxis()->SetTitle("E[eV]");
  gr1->SetMarkerColor(kBlack);
  gr1->SetMarkerStyle(20);
  gr1->SetMarkerSize(0.8);


  qv = fit_line->GetParameter(0);
  err_qv = fit_line->GetParError(0);
  mv = fit_line->GetParameter(1);
  err_mv = fit_line->GetParError(1);

  ca->cd(2);
  gr1->Draw("AP");
  gStyle->SetOptFit(0111);
  TPaveStats *st = (TPaveStats *)gr1->FindObject("stats");
  st->SetX1NDC(0.15);
  st->SetX2NDC(.45);


  //Resolutions
  double FWHM[num_peaks];
  double err_FWHM[num_peaks];
  double Energy[num_peaks];

  for (int i = 0; i < num_peaks; i++) {
    FWHM[i] = 2.35 * sqrt(par[1] * SiW * (par[num_def_par + 1 + i * 2] * mv) + (par[0] * par[0] / pow(2.35, 2)));
    err_FWHM[i] = (par_err[1] * SiW * (par[num_def_par + 1 + i * 2] * mv) + par[1] * SiW * (par_err[num_def_par + 1 + i * 2] * mv) + par[1] * SiW * (par[num_def_par + 1 + i * 2] * err_mv) + (2 * par[0] * par_err[0]) / pow(2.35, 2)) / (2 * FWHM[i]);
    Energy[i] = par[num_def_par + 1 + i * 2] * mv + qv;
  }


  TGraphErrors *gr_resolution = new TGraphErrors(num_peaks, Energy, FWHM, 0, err_FWHM);
  gr_resolution->GetXaxis()->SetTitle("E[eV]") ;
  gr_resolution->GetYaxis()->SetTitle("FWHM[eV]");
  gr_resolution->SetTitle("FWHM");
  gr_resolution->SetMarkerColor(kBlue);
  gr_resolution->SetMarkerStyle(20);
  gr_resolution->SetMarkerSize(0.8);
  ca->cd(4);
  gr_resolution->Draw("AP");


  //PullPlot
  TH1D *h_pull = new TH1D("h_pull", "", hist->GetNbinsX(), hist->GetXaxis()->GetXmin(), hist->GetXaxis()->GetXmax());
  for (int i = hist->FindBin(1500); i < hist->FindBin(4500); i++) {
    double w = sqrt((double)hist->GetBinContent(i));
    double bincont = hist->GetBinContent(i) - (double)fitfunc->GetBinContent(i);
    h_pull->SetBinContent(i, bincont / w);
  }
  h_pull->GetXaxis()->SetRangeUser(1500, 4500);
  gStyle->SetOptStat(0);
  h_pull->GetXaxis()->SetTitle("ADC[chn]");
  h_pull->GetYaxis()->SetTitle("Residuals/Sigma");
  h_pull->SetMarkerStyle(20);
  h_pull->SetMarkerSize(0.5);
  ca->cd(3);
  h_pull->Draw("P");

  TFile *outfile = new TFile(Form("rootfiles/20230405_bus%d_sdd_%d.root", bus_num, sdd_num), "recreate");
  outfile->cd();
  ca->Write();
  ca->Close();
  outfile->Close();

  TCanvas *canvas = new TCanvas();
  hist->Draw();
  fitfunc->Draw("same");
  fit_gaussKa1->Draw("same");
  fit_gaussKa2->Draw("same");
  fit_tail->Draw("same");
  fit_bkg->Draw("same");
  canvas->SaveAs(Form("plots/20230405_bus%d_sdd_%d.pdf", bus_num, sdd_num));

  TCanvas *canvas1 = new TCanvas();
  gPad->SetLogy();
  hist->Draw();
  fitfunc->Draw("same");
  fit_gaussKa1->Draw("same");
  fit_gaussKa2->Draw("same");
  fit_tail->Draw("same");
  fit_bkg->Draw("same");
  canvas1->SaveAs(Form("plots/20230405_bus%d_sdd_%d_log.pdf", bus_num, sdd_num));

  delete fitfunc;
  delete fit_gaussKa1;
  delete fit_gaussKa2;
  delete fit_tail;
  delete fit_bkg;

  std::ofstream file1;
  file1.open(Form("output/bus_%d_sdd%d_calib_par.txt", bus_num, sdd_num));
  file1 << "bus" << "\t" << "sdd" << "\t" << "chi_square_red" << "\t" << "gain" << "\t" << "err_gain" << "offset" << "\t" << "err_offset" << "\t" << "StausFit" << std::endl;
  file1 << bus_num << "\t" << sdd_num << "\t" << chi_square_red << "\t" << mv << "\t" << err_mv << "\t" << qv << "\t" << err_qv << "\t" << status_fit;
  file1 << "\n";
  file1.close();

  std::ofstream file2;
  file2.open(Form("output/bus_%d_sdd%d_fit_par_results.txt", bus_num, sdd_num));
  file2 << bus_num << "\t" << sdd_num << "\t";
  for (int i = 0; i < num_par; i++) {
    if (i > 3 && i < num_par - num_bkg_par && i % 2 == 0) {
      file2 << par[i] << "\t" << par_err[i] << "\t";
    }
    else {
      file2 << par[i] << "\t" << par_err[i] << "\t";
    }
  }
  file2 << "\n";
  file2.close();


  std::ofstream file3;
  file3.open(Form("output/bus_%d_sdd%d_resolutions.txt", bus_num, sdd_num));
  file3 << "bus" << "\t" << "sdd" << "\t" << "E_peak" << "\t" << "FWHM[eV]" << "\t" << "Err_FWHM[eV]" << std::endl;
  for (int i = 0; i < num_peaks; i++) {
    file3 << bus_num << "\t" << sdd_num << "\t" << Energy[i] << "\t" << FWHM[i] << "\t" << err_FWHM[i] << "\n";
  }
  file3.close();


}

std::vector<double> SDDCalibration::GetParams(
    nt bus_num, int sdd_num, const std::string &input_par)
{
  std::ifstream input_file(input_par);
  std::vector<double> param_tmp;

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

void SDDCalibration::FitSpectrum(TH1D *hist)
{
  // Define an array of initial parameters struct
  InitParameters init_par[num_par];

  // Set parameter values for noise
  init_par[0].value = 10.0;
  init_par[0].upp_bnd = init_par[0].value * 10.0;
  init_par[0].low_bnd = init_par[0].value * 0.5;
  init_par[0].step_size = (init_par[0].upp_bnd - init_par[0].low_bnd) * 0.01;

  // Set parameter values for Fano Factor
  init_par[1].value = 0.116;
  init_par[1].upp_bnd = init_par[1].value * 2.0;
  init_par[1].low_bnd = init_par[1].value * 0.5;
  init_par[1].step_size = (init_par[1].upp_bnd - init_par[1].low_bnd) * 0.01;

  // Set parameter values for the amplitude ratio of the tail function
  init_par[2].value = 0.03;
  init_par[2].upp_bnd = init_par[2].value * 3.0;
  init_par[2].low_bnd = init_par[2].value * 0.3;
  init_par[2].step_size = (init_par[2].upp_bnd - init_par[2].low_bnd) * 0.01;

  // Set parameter values for slope of the tail
  init_par[3].value = 5.0;
  init_par[3].upp_bnd = init_par[3].value * 3.0;
  init_par[3].low_bnd = init_par[3].value * 0.2;
  init_par[3].step_size = (init_par[3].upp_bnd - init_par[3].low_bnd) * 0.01;

  // Set parameter values for amplitudes and means of peaks
  for (int k = num_def_par; k < num_par - num_bkg_par; k = k + 2) {
    init_par[k].value = hist->GetBinContent(hist->FindBin(init_par[k + 1].value));
    init_par[k].upp_bnd = init_par[k].value * 1.5;
    init_par[k].low_bnd = init_par[k].value * 0.1;
    init_par[k].step_size = (init_par[k].upp_bnd - init_par[k].low_bnd) * 0.01;
  }
  for (int k = num_def_par + 1; k < num_par - num_bkg_par; k = k + 2) {
    init_par[k].upp_bnd = init_par[k].value + 20.0;
    init_par[k].low_bnd = init_par[k].value - 20.0;
    init_par[k].step_size = (init_par[k].upp_bnd - init_par[k].low_bnd) * 0.01;
  }
  // Set parameter values for background
  for (int k = num_par - num_bkg_par; k < num_par; ++k) {
    init_par[k].upp_bnd = init_par[k].value * 1.5;
    init_par[k].low_bnd = init_par[k].value * 0.5;
    init_par[k].step_size = (init_par[k].upp_bnd - init_par[k].low_bnd) * 0.1;
    init_par[k + 1].upp_bnd = init_par[k + 1].value * (init_par[k + 1].value > 0 ? 3.0 : 0.1);
    init_par[k + 1].low_bnd = init_par[k + 1].value * (init_par[k + 1].value > 0 ? 0.1 : 3.0);
    init_par[k + 1].step_size = (init_par[k + 1].upp_bnd - init_par[k + 1].low_bnd) * 0.1;
  }

  num_bins = hist->GetNbinsX();
  bin_content = std::make_unique<double[]>(num_bins);
  bin_center = std::make_unique<double[]>(num_bins);

  // Loop through all bins in the histogram
  for (int bin = 0; bin < num_bins; ++bin) {
    bin_content[bin] = hist->GetBinContent(bin + 1);
    bin_center[bin] = hist->GetBinCenter(bin + 1);
    if (bin_center[bin] < adc_min || bin_center[bin] > adc_max) {
      bin_content[bin] = 0.;
    }
  }
  
  // Initialize TMinuit and set parameter names and initial values  
  TMinuit *minuit = new TMinuit(num_par);

  double arglist[10];   // argument list
  int ierflg = 0;       // flag indicating whether the function was successful
  
  minuit->SetFCN(fcn);  // set the function to be minimized

  // Set initial parameter values, step sizes, and boundaries
  minuit->mnparm(0, "Noise", init_par[0].value, init_par[0].step_size, 
                init_par[0].low_bnd, init_par[0].upp_bnd, ierflg);
  minuit->mnparm(1, "FanoFactor", init_par[1].value, init_par[1].step_size, 
                init_par[1].low_bnd, init_par[1].upp_bnd, ierflg);
  minuit->mnparm(2, "RatioTail", init_par[2].value, init_par[2].step_size, 
                init_par[2].low_bnd, init_par[2].upp_bnd, ierflg);
  minuit->mnparm(3, "SlopeTail", init_par[3].value, init_par[3].step_size, 
                init_par[3].low_bnd, init_par[3].upp_bnd, ierflg);

  for (int j = num_def_par; j < num_par - num_bkg_par; ++j) {    // peaks parameters
    if (j % 2 == 0) {
      minuit->mnparm(j, "Amp_" + std::to_string((j + 2 - num_def_par) / 2),
                     init_par[j].value, init_par[j].step_size, 
                     init_par[j].low_bnd, init_par[j].upp_bnd, ierflg);
    }
    if (j % 2 != 0) {
      minuit->mnparm(j, "Mean_" + std::to_string((j + 1 - num_def_par) / 2),
                     init_par[j].value, init_par[j].step_size, 
                     init_par[j].low_bnd, init_par[j].upp_bnd, ierflg);
    }
  }

  for (int j = num_par - num_bkg_par; j < num_par; ++j) {   // background parameters
    minuit->mnparm(j, "Bkg_p" + std::to_string(j - 2 * num_peaks - num_def_par + 1),
                   init_par[j].value, init_par[j].step_size, 
                   init_par[j].low_bnd, init_par[j].upp_bnd, ierflg);
  }

  // Set optimisation options
  arglist[0] = 1.;    // standard error definition (1 sigma)
  minuit->mnexcm("SET ERR", arglist, 1, ierflg);
  arglist[0] = 1.;    // standard fit strategy
  minuit->mnexcm("SET STRATEGY", arglist, 1, ierflg);
  arglist[0] = 0.;    // turn off printout during fit
  minuit->mnexcm("SET PRINTOUT", arglist, 1, ierflg);

  arglist[0] = 200000000;   // maximum number of iterations
  arglist[1] = 0.00000001;  // tolerance for convergence

  // Run the minimisation algorithm (MIGRAD)
  minuit->mnexcm("MIGRAD", arglist, 2, ierflg);

  // Retrieve the parameter values and errors from the Minuit minimiser
  for (int p = 0; p < num_par; p++) {
    minuit->GetParameter(p, par[p], par_err[p]);
  }

  status_fit = minuit->fCstatu;

  ndf = 0;  // number of degrees of freedom
  chi_square = 0;

  double noise = par[0];
  double FF = par[1];
  double ratio_tail = par[2];
  double slope_tail = par[3];
  double amp_TiKa1  = par[num_def_par];
  double amp_CuKa1 = par[num_def_par + 2];
  double mean_TiKa1 = par[num_def_par + 1];
  double mean_CuKa1 = par[num_def_par + 3];
  double slope = (mean_CuKa1 - mean_TiKa1) / ((CuKa1 - TiKa1));
  double sigma = 0;
  double gauss_func = 0, gauss_func_Ka2 = 0;
  double tail_func = 0;
  double bkg_func = 0;
  double total_fit_func = 0;

  double x;

  for (int bin = 0; bin < num_bins; ++bin) {
    x = bin_center[bin];
    
    for (int g = num_def_par; g < num_par - num_bkg_par; g += 2) {
      sigma = TMath::Sqrt(FF * SiW * par[g + 1] * slope +
                          TMath::Power(noise / 2.35, 2));
      gauss_func += par[g] * TMath::Exp(-0.5 * 
                                        TMath::Power((x - par[g + 1]) / sigma, 2));
      tail_func += ratio_tail * par[g] *
                   TMath::Exp(((x - par[g + 1]) / (sigma * slope_tail)) +
                              (1 / (2 * TMath::Power(slope_tail, 2)))) *
                   TMath::Erfc(((x - par[g + 1]) / (TMath::Sqrt(2) * sigma)) +
                               1 / (TMath::Sqrt(2) * slope_tail));
    }
    
    double mean_TiKa2 = mean_TiKa1 - ((TiKa1 - TiKa2) * slope);
    double sigma_TiKa2 = TMath::Sqrt(FF * SiW * mean_TiKa2 * slope +
                                     TMath::Power(noise / 2.35, 2));
    double mean_CuKa2 = mean_CuKa1 - ((CuKa1 - CuKa2) * slope);
    double sigma_CuKa2 = TMath::Sqrt(FF * SiW * mean_CuKa2 * slope + 
                                     TMath::Power(noise / 2.35, 2));
                                     
    double gauss_TiKa2 = TiKa2_RI / TiKa1_RI * amp_TiKa1 *
                         TMath::Exp(-0.5 * TMath::Power((x - mean_TiKa2) / sigma_TiKa2, 2));                         
    double tail_TiKa2 = TiKa2_RI / TiKa1_RI * ratio_tail * amp_TiKa1 *
                        TMath::Exp(((x - mean_TiKa2) / (sigma_TiKa2 * slope_tail)) +
                                   (1 / (2 * TMath::Power(slope_tail, 2)))) *
                        TMath::Erfc(((x - mean_TiKa2) / (TMath::Sqrt(2) * sigma_TiKa2)) +
                                    1 / (TMath::Sqrt(2) * slope_tail));
    double gauss_CuKa2 = CuKa2_RI / CuKa1_RI * amp_CuKa1 *
                         TMath::Exp(-0.5 * TMath::Power((x - mean_CuKa2) / sigma_CuKa2, 2));
    double tail_CuKa2 = CuKa2_RI / CuKa1_RI * ratio_tail * amp_CuKa1 *
                        TMath::Exp(((x - mean_CuKa2) / (sigma_CuKa2 * slope_tail)) +
                                   (1 / (2 * TMath::Power(slope_tail, 2)))) * 
                        TMath::Erfc(((x - mean_CuKa2) / (TMath::Sqrt(2) * sigma_CuKa2)) +
                                    1 / (TMath::Sqrt(2) * slope_tail));
    
    tail_func += (tail_TiKa2 + tail_CuKa2);
    gauss_func_Ka2 = gauss_TiKa2 + gauss_CuKa2

    bkg_func = par[2 * num_peaks + num_def_par] + 
               TMath::Exp((x - par[2 * num_peaks + num_def_par + 1]) * par[2 * num_peaks + num_def_par + 2]);

    total_fit_func = gauss_func + gauss_func_Ka2 + tail_func + bkg_func;

    if (bin_content[bin] != 0) {
      ndf = ndf + 1;
      chi_square += TMath::Power(total_fit_func - bin_content[bin], 2) / bin_content[bin];
      fit_func->SetBinContent(bin + 1, total_fit_func);
      fit_gauss->SetBinContent(bin + 1, gauss_func);
      fit_gauss_Ka2->SetBinContent(bin + 1, gauss_func_Ka2);
      fit_tail->SetBinContent(bin + 1, tail_func);
      fit_bkg->SetBinContent(bin + 1, bkg_func);
    }
  }

  chi_square_red = chi_square / (ndf + 1);    // reduced chi-square
}

// Calculate the chi-square value for the given parameters
void fcn(int &npar, double *gin, double &f, double *par, int iflag)
{
  double x;

  // chi_square
  double chisq = 0;
  double delta = 0;
  ndf = 0;
  for (int i = 0; i < num_bins; i++) {
    if (bin_content[i] != 0) {
      x = bin_center[i];
      delta = (bin_content[i] - function(x, par)) / sqrt(bin_content[i]);
      chisq += delta * delta;
    }
  }
  f = chisq;
}

double function(double x, double *par)
{
  double noise = par[0];
  double FF = par[1];
  double ratio_tail = par[2];
  double slope_tail = par[3];
  double amp_TiKa1  = par[num_def_par];
  double amp_CuKa1 = par[num_def_par + 2];
  double mean_TiKa1 = par[num_def_par + 1];
  double mean_CuKa1 = par[num_def_par + 3];
  double slope = (mean_CuKa1 - mean_TiKa1) / ((CuKa1 - TiKa1));
  double sigma = 0;
  double gauss_func = 0, gauss_func_Ka2 = 0;
  double tail_func = 0;
  double bkg_func = 0;

  for (int g = num_def_par; g < num_par - num_bkg_par; g += 2) {
    sigma = TMath::Sqrt(FF * SiW * par[g + 1] * slope +
                        TMath::Power(noise / 2.35, 2));
    gauss_func += par[g] * TMath::Exp(-0.5 * 
                                      TMath::Power((x - par[g + 1]) / sigma, 2));
    tail_func += ratio_tail * par[g] *
                 TMath::Exp(((x - par[g + 1]) / (sigma * slope_tail)) +
                            (1 / (2 * TMath::Power(slope_tail, 2)))) *
                 TMath::Erfc(((x - par[g + 1]) / (TMath::Sqrt(2) * sigma)) +
                             1 / (TMath::Sqrt(2) * slope_tail));
  }
    
  double mean_TiKa2 = mean_TiKa1 - ((TiKa1 - TiKa2) * slope);
  double sigma_TiKa2 = TMath::Sqrt(FF * SiW * mean_TiKa2 * slope +
                                   TMath::Power(noise / 2.35, 2));
  double mean_CuKa2 = mean_CuKa1 - ((CuKa1 - CuKa2) * slope);
  double sigma_CuKa2 = TMath::Sqrt(FF * SiW * mean_CuKa2 * slope + 
                                   TMath::Power(noise / 2.35, 2));
                                     
  double gauss_TiKa2 = TiKa2_RI / TiKa1_RI * amp_TiKa1 *
                       TMath::Exp(-0.5 * TMath::Power((x - mean_TiKa2) / sigma_TiKa2, 2));                         
  double tail_TiKa2 = TiKa2_RI / TiKa1_RI * ratio_tail * amp_TiKa1 *
                      TMath::Exp(((x - mean_TiKa2) / (sigma_TiKa2 * slope_tail)) +
                                 (1 / (2 * TMath::Power(slope_tail, 2)))) *
                      TMath::Erfc(((x - mean_TiKa2) / (TMath::Sqrt(2) * sigma_TiKa2)) +
                                  1 / (TMath::Sqrt(2) * slope_tail));
  double gauss_CuKa2 = CuKa2_RI / CuKa1_RI * amp_CuKa1 *
                       TMath::Exp(-0.5 * TMath::Power((x - mean_CuKa2) / sigma_CuKa2, 2));
  double tail_CuKa2 = CuKa2_RI / CuKa1_RI * ratio_tail * amp_CuKa1 *
                      TMath::Exp(((x - mean_CuKa2) / (sigma_CuKa2 * slope_tail)) +
                                 (1 / (2 * TMath::Power(slope_tail, 2)))) * 
                      TMath::Erfc(((x - mean_CuKa2) / (TMath::Sqrt(2) * sigma_CuKa2)) +
                                  1 / (TMath::Sqrt(2) * slope_tail));
  gauss_func_Ka2 = gauss_TiKa2 + gauss_CuKa2;
  tail_func += (tail_TiKa2 + tail_CuKa2);
  
  bkg_func = par[2 * num_peaks + num_def_par] + 
             TMath::Exp((x - par[2 * num_peaks + num_def_par + 1]) * par[2 * num_peaks + num_def_par + 2]);
             
  return gauss_func + gauss_func_Ka2 + tail_func + bkg_func;
}

int main(int argc, char *argv[])
{
  int b = std::stoi(argv[1]);
  int s = std::stoi(argv[2]);

  std::string file_dir =
    "/home/aleks/SIDDHARTA2/SpectrumAnalyser/output/rootfiles/";

  std::string file_name = file_dir +
                          "histos_20230405_2143_0406_0101_25kv_100ua_cal.root";
  std::string input_par = "../calib_parameters/input_param_bus2.dat";

  calib(file_name, input_par, b, s);

  return 0;
}

