
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The function transform_gridpolygonsTC transforms nodes and elements of
% a nonuniform grid in the transformed coordinates to the original
% coordinates 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Input
% -----
% - xNU: x-coordinates of nodes (in transformed coordinates)
% - yNU: y-coordinates of nodes (in transformed coordinates)
% - xNUpolyarray: x-coordinates of element edges (in transformed coordinates)
% - yNUpolyarray: y-coordinates of element edges (in transformed coordinates)
% - thetabar: angle mean
% - Ebar: energy mean
% - Theta: covariance tensor

% Output
% ------
% - xNU: x-coordinates of nodes (in original coordinates)
% - yNU: y-coordinates of nodes (in original coordinates)
% - xNUpolyarray: x-coordinates of element edges (in original coordinates)
% - yNUpolyarray: y-coordinates of element edges (in original coordinates)

% Author: Pablo Seleson
% ------

% Last Modified: February 8, 2022
% -------------

function [xNU, yNU, xNUpolyarray, yNUpolyarray] = transform_gridpolygonsTC(xNU, yNU, xNUpolyarray, yNUpolyarray, thetabar, Ebar, Theta)

    % Perform transformation of nodes
    w = [xNU(:) yNU(:)];
    vbar = [thetabar Ebar];
    v = sqrtm(Theta)*w' + vbar';

    xNU = v(1,:);
    yNU = v(2,:);

    % Perform transformation of elements
    xNUpolyarray_temp = zeros(size(xNUpolyarray));
    yNUpolyarray_temp = zeros(size(yNUpolyarray));

    for j = 1:4
        w = [xNUpolyarray(j,:)' yNUpolyarray(j,:)'];
        v = sqrtm(Theta)*w' + vbar';

        xNUpolyarray_temp(j,:) = v(1,:);
        yNUpolyarray_temp(j,:) = v(2,:);
    end

    xNUpolyarray = xNUpolyarray_temp;
    yNUpolyarray = yNUpolyarray_temp;

end
