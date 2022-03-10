
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The function evaluate_IEAD_from_IEADstar computes IEAD values using 
% an IEAD* (IEAD in the transformed coordinates) distribution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Input
% -----
% - theta: array with angle values
% - E: array with energy values
% - thetabar: mean angle
% - Ebar: mean energy
% - Theta: covariance matrix
% - thetastar: array with theta* values
% - Estar: array with E* values
% - IEADstar: array with IEAD* values

% Output
% ------
% - IEAD: array with IEAD values

% Author: Pablo Seleson
% ------

% Last Modified: February 25, 2022
% -------------

function [IEAD] = evaluate_IEAD_from_IEADstar(theta,E,thetabar,Ebar,Theta,thetastar,Estar,IEADstar)

    % Dimension of data
    dim = length(theta);

    % Transform theta and E via the moments
    [wx,wy] = transform_E_theta(theta,E,thetabar,Ebar,Theta);
    
    % Find closest bins in the transformed coordinates
    [~,I] = min(abs(wx'-thetastar) + abs(wy'-Estar));
    I = I';

    % Initialize array
    IEAD = zeros(dim,1);
    
    % Area conversion factor
    factor = 1/(sqrt(abs(det(Theta))));
    
    % Assign data from the closest bin in the transformed coordinates
    for n = 1:dim
        IEAD(n) = factor * IEADstar(I(n));
    end
    