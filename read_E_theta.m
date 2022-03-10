
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The function read_E_theta reads particle velocity data and computes 
% the corresponding kinetic energy and angle values
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Input
% -----
% - data_file: name of data file 

% Output
% ------
% - E_values: array of kinetic energies
% - theta_values: array of angles

% Note:
% ----
% The function assumes the ion mass = proton mass (hardcoded)

% We assume cylindrical coordinates with vx being the "z-axis". 
% The radial component corresponds to the radius of the cylinder: 
%     vr = sqrt{vy^2 + vz^2}
% The angle relative to the "z-axis" is:
%     theta = tan^{-1}(vr/|vx|) 
% note vx is negative.

% Author: Pablo Seleson
% ------

% Last Modified: February 1, 2022
% -------------

function [E_values,theta_values] = read_E_theta(data_file)

        % Conversion factor: Joules -> eV
        Joules2eV = 6.242e+18;
        
        % Proton mass
        mp =  1.67262192369E-27; % [in kg]

        % --------------------------------------------------
        %            Read particle velocities
        % --------------------------------------------------

        % Read data from file
        data_particle_velocities = importdata(data_file);

        % Assign data to velocity arrays
        vx = data_particle_velocities(:,1);
        vy = data_particle_velocities(:,2);
        vz = data_particle_velocities(:,3);

        % --------------------------------------------------
        %              Compute E and theta
        % --------------------------------------------------

        % Ion mass
        m = mp;  

        % Kinetic energy
        E_values = Joules2eV*0.5*m*(vx.^2 + vy.^2 + vz.^2);

        % Compute radial velocity values (in cylindrical coordinates)
        vr = sqrt(vy.^2 + vz.^2);

        % Angle values
        theta_values = atan(vr./abs(vx));

end