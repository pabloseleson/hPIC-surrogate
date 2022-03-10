
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The function transform_E_theta transforms the values of E and theta
% for given moments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Input
% -----
% - theta: array of angle values
% - E: array of energy values
% - thetabar: mean of angle values
% - Ebar: mean of energy values
% - Theta: covariance matrix

% Output
% ------
% - thetastar: array of transformed angle values
% - Estar: array of transformed energy values

% Author: Pablo Seleson
% ------

% Last Modified: February 1, 2022
% -------------

function [thetastar,Estar] = transform_E_theta(theta,E,thetabar,Ebar,Theta)

    v = [theta E];
    vbar = [thetabar Ebar];
    
    % Compute transformed distribution
    w = (sqrtm(Theta))\(v-vbar)';
    w = w';
    
    % Assign transformed velocities components
    thetastar = w(:,1);
    Estar = w(:,2);

end