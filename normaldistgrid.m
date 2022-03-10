
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The function normaldistgrid generates a grid with N points uniformly
% distributed in (-1,1) and smoothly coarsened outside that range following
% a scaled normal distribution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Input
% -----
% - N: Number of points in (-1,1)
% - xL: Grid limits (-xL,xL)

% Discussion
% ----------
% Consider a normal distribution f(x). We would like the following
% relation to hold:
%
% \int_{x_n}^{x_{n+1}} f(x) dx = (1/N) \int_{-1}^{1} f(x) dx
%
% Let the cumulative distribution be:
%
% \Phi(x) = \int_{-inft}^{x} f(y) dy
%
% Then, we have 
%
% \Phi(x_{n+1}) - \Phi(x_{n}) = (1/N) (\Phi(1) - \Phi(-1))
%
% or
%
% \Phi(x_{n+1}) = \Phi(x_{n}) + (1/N) (\Phi(1) - \Phi(-1))
%
% Link: https://www.mathworks.com/help/stats/normal-distribution.html
%
% Note (see link) that
% 
% \Phi(x) = (1/2) ( 1 - erf(-x/sqrt(2)) )
%
% Then
%
% erf(-x/sqrt(2)) = 1 - 2\Phi(x) => x = -sqrt(2) erfinv( 1 - 2\Phi(x) )

% Author: Pablo Seleson
% ------

% Last Modified: February 1, 2022
% -------------

function [xgrid] = normaldistgrid(N,xL)
    
    % Normal distribution
    syms x
    
    % Note: standard normal cumulative distribution function:
    %       Phi = @(x) 0.5*(1-erf(-x/sqrt(2)));
    
    rf = 3/5; % Scaling factor
    Phi = @(x) 0.5*(1-erf(-rf*x/sqrt(2))); % Scaled by rf

    % Generate mesh in [-1,1]
    dx = 2/N;
    xmesh = -1:dx:1;
    
    % RHS
    RHS = (Phi(1) - Phi(-1))/N;
    
    % Initialize x_n
    xn = 1;
    
    % Initialize grid
    rgrid = xn;
    
    while xn < xL
        
       % Note: We need xnp1 s.t. Phi(xnp1) = Phi(xn) + RHS
       % 
       % This is equivalent to (without scaling): 
       % 
       %     xnp1 = -sqrt(2) erfinv( 1 - 2\Phi(xnp1) ) = -sqrt(2) erfinv( 1 - 2(Phi(xn) + RHS))
       
        xnp1 = -(sqrt(2)/rf)*erfinv(1 - 2*(Phi(xn) + RHS)); % Includes scaling factor rf
        
        % Add entry to grid
        rgrid = [rgrid xnp1];
        
        % Update xn value
        xn = xnp1;
        
    end
    
    % Remove first entry (1) and last entry (NaN)
    rgrid = rgrid(2:end-1);
    
    % Define left grid
    lgrid = fliplr(-rgrid);
    
    % Concatenate grids
    xgrid = [lgrid xmesh rgrid];

end
        