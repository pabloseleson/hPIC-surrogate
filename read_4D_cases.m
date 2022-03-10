
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The function read_4D_cases reads the cases to read data from in 4D
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Input
% -----
% - cases_file: name of file to read cases from

% Output
% ------
% - filename: array of file names
% - case_number: array of case numbers
% - rep_number: array of repetition numbers
% - Log10_Te_Ti: array of log_{10} of temperature ratios
% - Te: array of electron temperatures
% - Ti: array of ion temperatures
% - Psi: array of magnetic field angles
% - B: array of magnetic field magnitudes
% - Log10_n_values: array of log_{10} of density values
% - n_values: array of density values

% Author: Pablo Seleson
% ------

% Last Modified: February 1, 2022
% -------------

function [filename,case_number,rep_number,Log10_Te_Ti,Te,Ti,Psi,B,Log10_n_values,n_values] = read_4D_cases(cases_file)

    % Open and read file
    fid = fopen(cases_file,'r');
    indata = textscan(fid, '%s %f %f %f %f %f %f %f %f %f', 'HeaderLines',1); % Number of lines to skip
    fclose(fid);

    % Assign data to arrays
    filename = indata{1};
    case_number = indata{2};
    rep_number = indata{3};
    Log10_Te_Ti = indata{4};
    Te = indata{5};
    Ti = indata{6};
    Psi = indata{7};
    B = indata{8};
    Log10_n_values = indata{9};
    n_values = indata{10};
    
end