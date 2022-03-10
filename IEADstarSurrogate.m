
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The function IEADstarSurrogate fits the IEAD* (IEAD computed in the 
% transformed coordinates) using data constructed with the moments 
% surrogate model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Input
% -----
% sglevel: data sparse grid level
% dim: dimension of parameter space
%      2D: dim = 2
%      4D: dim = 4

% Author: Pablo Seleson
% ------

% Last Modified: February 2, 2022
% -------------

function IEADstarSurrogate(sglevel,dim)
    
    % ====================================================================
    %                       Load IEAD Data
    % ====================================================================
    
    IEAD_data_filename = ['IEAD_data_' num2str(dim) 'D/Level_' num2str(sglevel) '_train.mat'];
    
    load(IEAD_data_filename)
    
    % Number of pixels
    dimXstar = length(XNODESstar);
    dimYstar = length(YNODESstar);
    dim_bins_star = dimXstar * dimYstar;
    
    % ====================================================================
    %           Create Basis Functions: Lagrange Polynomials
    % ====================================================================
    
    iDepth_basis = sglevel-2;
    sType = 'iptotal';
    iOut = dim_bins_star;

    gridname = ['LS_' num2str(dim) 'D_Grid_IEADstar_level_' num2str(sglevel)];

    if dim == 2
        % 2D case
        [lGrid_IEADstar, ~] = tsgMakeGlobal(gridname,2,iOut,'clenshaw-curtis',sType,iDepth_basis,[min(points(:,1)) max(points(:,1)); min(points(:,2)) max(points(:,2))]);
    else
        % 4D case
        [lGrid_IEADstar, ~] = tsgMakeGlobal(gridname,4,iOut,'clenshaw-curtis',sType,iDepth_basis,[min(points(:,1)) max(points(:,1)); min(points(:,2)) max(points(:,2)); min(points(:,3)) max(points(:,3)); min(points(:,4)) max(points(:,4))]);
    end
    
    % ====================================================================
    %                      Find Coefficients
    % ====================================================================
    
    % Evaluate basis functions at grid points
    % Note: B has dimensions [Npoints x Nbasis]
    % ----
    B = tsgEvaluateHierarchy(lGrid_IEADstar,points);

    % Solve coefficients
    % Note: rhs has dimensions[Npoints x iOut]
    % ----
    % Two cases:
    % (1) If A is an N-by-N matrix and B is a column vector with N components, 
    %     or a matrix with several such columns, then X = A\B is the solution 
    %     to the equation A*X = B.
    % (2) If A is an M-by-N matrix with M < or > N and B is a column vector 
    %     with M components, or a matrix with several such columns, then 
    %     X = A\B is the solution in the least-squares sense to the under- 
    %     or over-determined system of equations A*X = B.
 
    rhs = [IEADstar_array];
    a = B\rhs; % This solves B*a = rhs 
    
    % Load coefficients to the grid object
    tsgLoadHCoefficients(lGrid_IEADstar, a);

end
  
     
