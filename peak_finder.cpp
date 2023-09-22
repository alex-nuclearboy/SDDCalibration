#include "SDDCalibration.h"
#include <iostream>

#include "TApplication.h"

std::string CheckExtension(const std::string& filename);

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
  const std::string data_directory  = "${XRAYSPECTRA}";
  const std::string input_file_name = CheckExtension(argv[1]);
  const std::string input_file_path = data_directory + "/" + input_file_name;
  const int bus_number = std::stoi(argv[2]);
  const int sdd_number = std::stoi(argv[3]);
  const std::string output_file_name = input_file_name.substr(7, 23);
  
  SDDCalibration pf(input_file_path);  
  pf.DrawHistogram(bus_number, sdd_number, output_file_name);
  
  std::cout << "End." << std::endl;

  app.Run();

  return 0;
}

std::string CheckExtension(const std::string& filename) {
  if (filename.substr(filename.size() - 5) != ".root") {
    return filename + ".root";
  }
  return filename;
}
