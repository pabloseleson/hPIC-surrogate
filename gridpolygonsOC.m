
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The function gridpolygonsOC creates a grid & polygons in the original 
% coordinates
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Output
% ------
% XNODES: 1D array of theta values
% YNODES: 1D array of E/Te values
% dimX: Number of theta values  
% dimY: Number of E/Te values 
% x: 2D array of theta grid values
% y: 2D array of E/Te grid values
% xpolyarray: 4 x (dimX*dimY) array of theta values of polygons vertices
% ypolyarray: 4 x (dimX*dimY) array of E/Te values of polygons vertices
% element: 2D array of element numbers

% Author: Pablo Seleson
% ------

% Last Modified: February 1, 2022
% -------------

function [XNODES, YNODES, dimX, dimY, x, y, xpolyarray, ypolyarray, element] = gridpolygonsOC

    % Number of points in each dimension
    Ntheta = 90;
    NE = 240;
    
    % Array of theta values
    XEDGES = linspace(0, pi*90/180, Ntheta+1);
    XNODES = (XEDGES(1:end-1) + XEDGES(2:end))/2;
    dimX = length(XNODES);
    
    % Array of E/Te values
    YEDGES = linspace(0, 24, NE+1); % Max energy is 24*Te
    YNODES = (YEDGES(1:end-1) + YEDGES(2:end))/2;
    dimY = length(YNODES);
    
    % Create grid
    x = zeros(dimX,dimY);
    y = zeros(dimX,dimY);
    
    % Initialize polygons arrays
    element = zeros(dimX,dimY);
    xpolyarray =  zeros(4,dimX*dimY);
    ypolyarray =  zeros(4,dimX*dimY);
    
    % Initialize element counter
    count = 0;
    
    % Run over theta and E/Te values
    for j = 1:dimY
        for i = 1:dimX
            
            % --------
            %   Grid
            % --------

            x(i,j) = XNODES(i);
            y(i,j) = YNODES(j);
            
            % --------
            % Polygons
            % --------
            count = count + 1;
            
            % xleft
            if i == 1
                xleft =  XEDGES(1);
            else
                xleft =  0.5*(XNODES(i)+XNODES(i-1)); 
            end
            
            % xright
            if i == dimX
                xright = XEDGES(end);
            else
                xright =  0.5*(XNODES(i)+XNODES(i+1)); 
            end

            % ydown
            if j == 1
                ydown =  YEDGES(1);
            else
                ydown =  0.5*(YNODES(j)+YNODES(j-1)); 
            end
            
            % yup
            if j == dimY
                yup =  YEDGES(end);
            else
                yup =  0.5*(YNODES(j)+YNODES(j+1)); 
            end
            
            element(i,j) = count;
            xpolyarray(:,count) = [xleft; xright; xright; xleft];
            ypolyarray(:,count) = [ydown; ydown; yup; yup];

        end
    end
    
end   