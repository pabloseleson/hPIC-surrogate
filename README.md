
CODE DESCRIPTION 
=================

The hPIC-surrogate MATLAB code creates surrogate ion energy-angle distribution 
(IEAD) models based on hPIC data and performs sensitivity analysis.

CODE DEPENDENCIES
=================

The code requires the following:

(1) MATLAB files

(2) Tasmanian: see installation instructions below

(3) hPIC data: see downloading instructions below


MATLAB FILES
============

The function "main.m" is the main function that needs to be run. 

The function has two hardcoded parameters: the physical parameter space dimension 
(dim = 2) and the sparse grid level (sglevel = 7). These parameters can be 
changed to account for other cases, with the following options:
- dim = 2, 4
- sglevel = 7, 10, 13

The function requires definition of the following directories:
- Directory_tsg: directory corresponding to the Tasmanian Workfiles folder
- Directory_train: directory containing the training data 
                   (two directories need to be defined - one for 2D and one for 4D)
- Directory_test: directory containing the testing data 
                   (two directories need to be defined - one for 2D and one for 4D)                          


TASMANIAN INSTALLATION
======================

The following are steps for installation of Tasmanian and linking MATLAB to it.

STEP 1: Download the source code (TASMANIAN-7.7.tar.gz) file from: 

https://github.com/ORNL/TASMANIAN/releases/tag/v7.7

and move the file to a desired directory, e.g., "~/Documents/MATLAB/"

STEP 2: Untar source code file by running

tar -xf TASMANIAN-7.7.tar.gz

or, if the file is automatically unzipped,  

tar -xf TASMANIAN-7.7.tar

STEP 3: Run the following commands 

cd TASMANIAN-7.7

mkdir build

cd build

cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=../INSTALL_DIR -DTasmanian_MATLAB_WORK_FOLDER=../Workfiles/ ..

Note: You would need to have cmake installed (see https://cmake.org/download/ 
      for downloading instructions) and add it to the path. To add cmake to the 
      path you may use the following command:

PATH="/Applications/CMake.app/Contents/bin":"$PATH"

STEP 4: Run the following commands 

cmake --build .

cmake --build . --target install

STEP 5: Add Tasmanian to the MATLAB path 

Type the following in MATLAB (a command to add the path is output when running cmake --build . --target install)

addpath('~/Documents/MATLAB/TASMANIAN-7.7/INSTALL_DIR/share/Tasmanian/matlab/')


hPIC DATA 
=========

The hPIC data can be freely accessed online through the Constellation portal 
(https://doi.ccs.ornl.gov/) of the Oak Ridge National Laboratory Leadership 
Computing Facility under the following doi: 10.13139/ORNLNCCS/1846780 

