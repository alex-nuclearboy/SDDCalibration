#include "./include/SDDCalibration.h"
#include <iostream>

int main(int argc, char *argv[]) 
{
  if (argc != 3) {
    std::cerr << "Usage: " << argv[0] 
              << " <input_root_file_name> <input_peaks_file> " << std::endl;
    return 1;
  }

  std::cout << "Starting the program..." << std::endl;
  std::cout << std::endl;

  std::string input_peaks_file_name = argv[2];
  std::string input_peaks_file = 
      "output/precalibration/parameters/" + \
      input_peaks_file_name;
  const std::string data_directory  = "${XRAYSPECTRA}";
  std::string input_root_file_name = argv[1];
  const std::string input_root_file = data_directory + "/" + input_root_file_name;
  // std::string output_file_name = input_peaks_file_name.substr(7, 23);

  SDDCalibration ff(input_root_file);  
  ff.PrintParams(input_peaks_file);

  std::cout << "End." << std::endl;

  return 0;
}
