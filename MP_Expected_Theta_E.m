
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The function MP_Expected_Theta_E computes the most probable and expected
% angle and energy
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Input
% -----
% - xNU: x-coordinates of nodes 
% - yNU: y-coordinates of nodes 
% - xNUpolyarray: x-coordinates of element edges 
% - yNUpolyarray: y-coordinates of element edges 
% - IEAD_NU: ion energy-angle distribution at the nodes
% - Te: electron temperature

% Output
% ------
% - MP_theta: most probably angle
% - MP_E: most probable energy
% - expected_theta: angle mean
% - expected_E: energy mean
% - std_theta: angle standard deviation
% - std_E: energy standard deviation

% Author: Pablo Seleson
% ------

% Last Modified: February 8, 2022
% -------------

function [MP_theta, MP_E, expected_theta, expected_E, std_theta, std_E] = MP_Expected_Theta_E(xNU,yNU,xNUpolyarray,yNUpolyarray,IEAD_NU,Te)

    % -------------------------------------------
    %  Most probable and expected angle & energy
    % -------------------------------------------

    % Compute area of polygons: Shoelace formula
    %
    % Link: https://en.wikipedia.org/wiki/Shoelace_formula
    %
    % xx = xNUpolyarray(:,n); yy = yNUpolyarray(:,n);
    %
    % Area = 0.5*(sum(xx.*[yy(2:end); yy(1)]) - sum(yy.*[xx(2:end); xx(1)]));
    %
    % polyarea(xx,yy);

    % Compute most probable (angle, energy)
    [~, I] = max(IEAD_NU);
    MP_theta = xNU(I);
    MP_E = yNU(I)*Te; % Convert value from E/T_e -> E

    % Compute expected (angle, energy)
    expected_theta = 0;
    expected_E = 0;
    Ncount = 0;
    
    % Define 1D arrays
    xNU_array = xNU(:);
    yNU_array = yNU(:);

    for n = 1:length(IEAD_NU)

        % Omit negative energy values
        if yNU_array(n) < 0

        else

            % Read bin vertices
            xx = xNUpolyarray(:,n);
            yy = yNUpolyarray(:,n);

            % Compute bin area
            bin_area = polyarea(xx,yy);

            % Bin count
            count = IEAD_NU(n)*bin_area;

            % Update arrays for expected angle and energy
            expected_theta = expected_theta + xNU_array(n)*count;
            expected_E  = expected_E + yNU_array(n)*count;

            % Update total count
            Ncount = Ncount + count;

        end
    end

    % Expected angle
    expected_theta = expected_theta/Ncount;

    % Expected energy
    expected_E = (expected_E/Ncount)*Te; % Convert E/T_e -> E

    % Compute standard deviation for (angle, energy)
    std_theta = 0;
    std_E = 0;
    Ncount = 0;

    for n = 1:length(IEAD_NU)

        % Omit negative energy values
        if yNU_array(n) < 0

        else

            % Read bin vertices
            xx = xNUpolyarray(:,n);
            yy = yNUpolyarray(:,n);

            % Compute bin area
            bin_area = polyarea(xx,yy);

            % Bin count
            count = IEAD_NU(n)*bin_area;

            % Update arrays for angle and energy standard deviations
            std_theta = std_theta + (xNU_array(n) - expected_theta)^2*count;
            std_E = std_E + (yNU_array(n) - expected_E/Te)^2*count;

            % Update total count
            Ncount = Ncount + count;

        end
    end

    % Angle standard deviation
    std_theta = sqrt(std_theta/(Ncount)); 

    % Energy standard deviation
    std_E = sqrt(std_E/(Ncount))*Te;

end
