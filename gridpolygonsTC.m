
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The function gridpolygonsTC creates a grid & polygons in the transformed 
% coordinates
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Output
% ------
% XNODESstar: 1D array of theta* values
% YNODESstar: 1D array of E* values
% dimXstar: Number of theta* values 
% dimYstar: Number of E* values 
% xstar: 2D array of theta* grid values
% ystar: 2D array of E* grid values
% xstarpolyarray: 4 x (dimXstar*dimYstar) array of theta* values of polygons vertices
% ystarpolyarray: 4 x (dimXstar*dimYstar) array of E* values of polygons vertices
% elementstar: 2D array of element numbers

% Author: Pablo Seleson
% ------

% Last Modified: February 1, 2022
% -------------

function [XNODESstar, YNODESstar, dimXstar, dimYstar, xstar, ystar, xstarpolyarray, ystarpolyarray, elementstar] = gridpolygonsTC
    
    N = 25; % Number of points in [-1,1];
    xL = 5; % Number of standard deviations
    
    % Array of theta* values
    [XNODESstar] = normaldistgrid(N,xL);
    dimXstar = length(XNODESstar);
    
    % Array of E* values
    YNODESstar = XNODESstar;
    dimYstar = length(YNODESstar); 
    
    % Initialize grid arrays
    xstar = zeros(dimXstar,dimYstar);
    ystar = zeros(dimXstar,dimYstar);

    % Initialize polygons arrays
    elementstar = zeros(dimXstar,dimYstar);
    xstarpolyarray =  zeros(4,dimXstar*dimYstar);
    ystarpolyarray =  zeros(4,dimXstar*dimYstar);
    
    % Initialize element counter
    count = 0;
    
    % Run over theta* and E* values
    for j = 1:dimYstar
        for i = 1:dimXstar
            
            % --------
            %   Grid
            % --------

            xstar(i,j) = XNODESstar(i);
            ystar(i,j) = YNODESstar(j);
            
            % --------
            % Polygons
            % --------
            count = count + 1;
            
            % xleft
            if i == 1
                xleft =  -5;
            else
                xleft =  0.5*(XNODESstar(i)+XNODESstar(i-1)); 
            end
            
            % xright
            if i == dimXstar
                xright =  5;
            else
                xright =  0.5*(XNODESstar(i)+XNODESstar(i+1)); 
            end

            % ydown
            if j == 1
                ydown =  -5;
            else
                ydown =  0.5*(YNODESstar(j)+YNODESstar(j-1)); 
            end
            
            % yup
            if j == dimYstar
                yup =  5;
            else
                yup =  0.5*(YNODESstar(j)+YNODESstar(j+1)); 
            end
            
            elementstar(i,j) = count;
            xstarpolyarray(:,count) = [xleft; xright; xright; xleft];
            ystarpolyarray(:,count) = [ydown; ydown; yup; yup];

        end
    end
    
end