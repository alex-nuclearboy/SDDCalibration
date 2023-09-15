# Energy Calibration of SDDs

This repository contains the code used in the energy calibration of Silicon Drift Detectors (SDDs) in the <a href="https://iopscience.iop.org/article/10.1088/1402-4896/ac7fc0">SIDDHARTA-2 experiment</a> for high-precision X-ray spectroscopy measurements of kaonic atoms. 

The SDDs are placed around a cryogenic target cell made of high purity aluminum structure and 75 µm thick Kapton walls. 
The <a href="https://iopscience.iop.org/article/10.1088/1402-4896/ac95da">energy calibration</a> of the SDDs is performed using two X-ray tubes installed on two sides of the vacuum chamber and a multi-element target made of high purity titanium and copper strips placed on the target cell walls. The X-ray tubes induce the fluorescence emission of the target elements, and the characteristic Kα and Kβ transitions are detected by the SDDs. 

The energy response function of the detector is predominantly a Gaussian curve for every fluorescence X-ray peak. However, the response also has a low-energy component due to incomplete charge collection and electron-hole recombination. 
Thus, the total peak fit function is formed by two contributions: Gauss function (the main contribution to the peak shape) and tail function (an exponential function to reproduce the incomplete charge collection).

Only the Ti Kα and Cu Kα peaks are used to calibrate the detectors, as they have the highest signal-to-background ratio. 


This software package includes functionalities to find peaks on SDD spectra, conduct preliminary calibrations, and perform main fitting procedures using MINUIT. 

## Requirements

    ROOT: A data analysis framework
    CMake: Cross-platform build-system generator

## Installation & Compilation

Clone the repository

    git clone https://github.com/alex-nuclearboy/SDDCalibration.git

Change to the project directory

    cd SDDCalibration

Create a build directory and change to it

    mkdir SDD-build
    cd SDD-build

Generate Makefile and compile

    cmake ..
    make

## Usage

### Peak Finding

The automatic procedure to search for peaks on the spectra is as follows:

    ./peakfinder <ROOT file's name> <bus number> <sdd number>

This will list all the peaks found in the spectra. Out of these, select the 12 most intense peaks and record them in a table in the following order: Ti Kα, Cu Kα, Ti Kβ, Cu Kβ, Fe Kα, and the remaining peaks.

### Preliminary Calibration

To perform a preliminary (quick) calibration of the spectra, use the following command:

    ./precalibration <ROOT file's name> <peaks file's name>

This will produce a table with initial parameters needed for the main calibration procedure.

### Main Fitting Procedure

The main fitting procedure uses MINUIT and is initiated by calling the MIGRAD function:

    ./calibration <ROOT file's name> <init parameters file's name>
