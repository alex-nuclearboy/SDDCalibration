#include "SDDCalibration.h"
#include <iostream>

int main(int argc, char *argv[]) 
{
  if (argc != 3) {
    std::cerr << "Usage: " << argv[0] 
              << " <input_root_file> <input_peaks_file> " << std::endl;
    return 1;
  }

  std::cout << "Starting the program..." << std::endl;
  std::cout << std::endl;

  std::string input_root_file_name = argv[1];
  std::string input_root_file = 
      "/home/aleks/SIDDHARTA2/SpectrumAnalyser/output/rootfiles/" + \
      input_root_file_name;
  std::string input_peaks_file_name = argv[2];
  std::string input_peaks_file = 
      "output/precalibration/parameters/" + input_peaks_file_name;  
  std::string output_file_name = input_root_file_name.substr(7, 23);

  SDDCalibration ff(input_root_file);  
  ff.PreFitPeaks(input_peaks_file, output_file_name);

  std::cout << "End." << std::endl;

  return 0;
}
