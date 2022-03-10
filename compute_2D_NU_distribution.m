
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The function compute_2D_NU_distribution computes a 2D distribution
% based on a (regular) non-uniform grid
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Input
% -----
% - wx: x-component of data
% - wy: y-component of data
% - XNODES: x-coordinates of nodal points in distribution 
% - YNODES: y-coordinates of nodal points in distribution  
% - xmin: minimum distribution domain in x-direction
% - xmax: maximum distribution domain in x-direction
% - ymin: minimum distribution domain in y-direction
% - ymax: maximum distribution domain in y-direction

% Output
% ------
% - x: array with x-components for distribution
% - y: array with y-components for distribution
% - f: distribution 
% - parea: array with bin areas

% Author: Pablo Seleson
% ------

% Last Modified: March 9, 2022
% -------------

function [x,y,f,parea] = compute_2D_NU_distribution(wx,wy,XNODES,YNODES,xmin,xmax,ymin,ymax)

    % Filter data outside domain limits
    tol = 1e-14;
    wxtemp = wx(wx > xmin-tol & wx < xmax+tol & wy > ymin-tol & wy < ymax+tol );
    wytemp = wy(wx > xmin-tol & wx < xmax+tol & wy > ymin-tol & wy < ymax+tol );
    wx = wxtemp;
    wy = wytemp;

    % Grid dimensions
    dimX = length(XNODES);
    dimY = length(YNODES); 
    
    % Grid resolutions
    dXNODES = zeros(dimX,1);
    dYNODES = zeros(dimY,1);
    
    % Compute element lengths per dimension
    for i = 1:dimX
        if i == 1
            dXNODES(i) = XNODES(i+1) - XNODES(i);
        elseif i == dimX
            dXNODES(i) = XNODES(i) - XNODES(i-1);
        else
            dXNODES(i) = 0.5*(XNODES(i+1) - XNODES(i)) + 0.5*(XNODES(i) - XNODES(i-1));
        end
    end
    
    for j = 1:dimY
        if j == 1
            dYNODES(j) = YNODES(j+1) - YNODES(j);
        elseif j == dimY
            dYNODES(j) = YNODES(j) - YNODES(j-1);
        else
            dYNODES(j) = 0.5*(YNODES(j+1) - YNODES(j)) + 0.5*(YNODES(j) - YNODES(j-1));
        end
    end
    
    % Data dimension
    dim_data = length(wx);
    
    % Find index in wx for distribution
    [~,Xminind] = min(abs(wx'-XNODES'));

    % Find index in wy for distribution
    [~,Yminind] = min(abs(wy'-YNODES'));

    % Initialize distribution
    x = zeros(dimX,dimY);
    y = zeros(dimX,dimY);
    f = zeros(dimX,dimY);
    parea = zeros(dimX,dimY);

    % Create grid
    for i = 1:dimX
        for j = 1:dimY
            x(i,j) = XNODES(i);
            y(i,j) = YNODES(j);
        end
    end
    
    % Create distribution
    for n = 1:dim_data
        i = Xminind(n);
        j = Yminind(n);        
        f(i,j) = f(i,j) + 1;
    end

    % Compute density value
    for i = 1:dimX
        for j = 1:dimY
            dX = dXNODES(i);
            dY = dYNODES(j);
            parea(i,j) = dX*dY;
            f(i,j) = f(i,j)/parea(i,j);
        end
    end

end
