
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The function moments_E_theta computes the mean of E and theta and the 
% covariance matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Input
% -----
% - theta: array of angle values
% - E: array of energy values

% Output
% ------
% - thetabar: mean of angle values
% - Ebar: mean of energy values
% - Theta: covariance matrix

% Author: Pablo Seleson
% ------

% Last Modified: February 1, 2022
% -------------

function [thetabar,Ebar,Theta] = moments_E_theta(theta,E)

    % Compute number of data
    N = length(theta);
    
    % Compute average quantities
    thetabar = sum(theta)/N;
    Ebar = sum(E)/N;
  
    % Compute covariance matrix
    Theta = zeros(2,2);

    Theta(1,1) = sum((theta - thetabar).^2)/N;
    Theta(2,2) = sum((E - Ebar).^2)/N;
    
    Theta(1,2) = sum((theta - thetabar).*(E - Ebar))/N;
    Theta(2,1) = Theta(1,2);

end