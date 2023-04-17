#include "SDDCalibration.h"
#include <iostream>

#include "TApplication.h"

int main(int argc, char *argv[]) 
{
  if (argc != 4) {
    std::cerr << "Usage: " << argv[0] 
              << " <input_root_file> <bus_number> <sdd_number>" << std::endl;
    return 1;
  }

  TApplication app("myapp", 0, 0);

  std::cout << "Starting the Peak Finder program..." << std::endl;
  std::cout << std::endl;

  std::string input_root_file_name = argv[1];
  std::string input_root_file = 
      "/home/aleks/SIDDHARTA2/SpectrumAnalyser/output/rootfiles/" + \
      input_root_file_name;
  int bus_number = std::stoi(argv[2]);
  int sdd_number = std::stoi(argv[3]);
  std::string output_file_name = input_root_file_name.substr(7, 23);
  std::string output_image_file = \
      std::string("plots/pf_") + output_file_name + std::string("_bus") \
      + std::string(argv[2]) + std::string("_sdd") + std::string(argv[3]) \
      + std::string(".png");

  SDDCalibration pf(input_root_file);  
  pf.DrawHistogram(bus_number, sdd_number, output_image_file);

  std::cout << "End." << std::endl;

  app.Run();

  return 0;
}
