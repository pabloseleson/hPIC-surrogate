
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The function GSA_Moments_2D computes the angle and energy moments profiles 
% (mean and standard deviation) using the 2D IEAD surrogate model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Input
% -----
% sglevel: data sparse grid level
% Npoints: number of points per dimension 

% Author: Pablo Seleson
% ------

% Last Modified: February 10, 2022
% -------------

function GSA_Moments_2D(sglevel, Npoints)

    % Check sglevel input
    if sglevel~= 7 && sglevel~= 10 && sglevel~= 13
        error('sglevel should be 7, 10, or 13.')
    end

    % ====================================================================
    %         Create Evaluation Grid in Physical Parameter Space
    % ====================================================================

    % Grid limits
    min_Log10_Te_Ti = -1.25;
    max_Log10_Te_Ti = 1.25;

    min_Psi = 0;
    max_Psi = 87;

    % Generate grid
    x_grid = linspace(min_Log10_Te_Ti,max_Log10_Te_Ti,Npoints);
    y_grid = linspace(min_Psi,max_Psi,Npoints);
    [Xgrid, Ygrid] = meshgrid(x_grid,y_grid);

    % Read physical parameter values
    Log10_Te_Ti = Xgrid(:);
    Psi = Ygrid(:);

    % Number of total points
    npoints = length(Log10_Te_Ti);

    % ====================================================================
    %                    Evaluate IEAD* model
    % ====================================================================

    % Reload 2D IEAD* fit grid
    gridname = ['LS_2D_Grid_IEADstar_level_' num2str(sglevel)];
    [lGrid_IEADstar] = tsgReloadGrid(gridname);

    % Evaluate IEAD* surrogate model on data points
    [IEADstar_all] = tsgEvaluate(lGrid_IEADstar, [Log10_Te_Ti Psi]);
    
    % ====================================================================
    %            Compute Moments per Case for Transformation
    % ====================================================================
    
    % Reload 2D moments grid
    gridname = ['LS_2D_Grid_Moments_level_' num2str(sglevel)];
    [lGrid_Moments] = tsgReloadGrid(gridname);
    
    % Moment
    [Moments] = tsgEvaluate(lGrid_Moments, [Log10_Te_Ti Psi]); 
    
    % Assign moments to arrays
    thetabar_mean = Moments(:,1);
    Ebar_mean = 10.^Moments(:,2);
    Theta_mean_11 = 10.^Moments(:,3);
    Theta_mean_22 = 10.^Moments(:,4);
    Theta_mean_12 = Moments(:,5);
    Theta_mean = [Theta_mean_11 Theta_mean_22 Theta_mean_12]; % [Theta(1,1) Theta(2,2) Theta(1,2)];

    % ====================================================================
    %              Read Binning in Transformed Coordinates
    % ====================================================================

    % Create elements for patch function
    [~, ~, ~, ~, xstar, ystar, ~, ~, ~] = gridpolygonsTC;

    % ====================================================================
    %                      Compute Moments 
    % ====================================================================

    % Moments
    Moments_Theta_E = zeros(npoints,4);

    % Ion temperature
    Ti = 10;

    % Create reference non-uniform grid in the transformed coordinates
    [~, ~, ~, ~, xNU_ref, yNU_ref, xNUpolyarray_ref, yNUpolyarray_ref, ~] = gridpolygonsTC;

    % Run over samples
    for i = 1:npoints

        fprintf('Sample %5g of %5g \n',i,npoints)

        % Read IEAD* for the ith sample
        IEADstar = IEADstar_all(i,:);

        % Moments to convert data
        thetabar = thetabar_mean(i);
        Ebar = Ebar_mean(i);
        Theta = [Theta_mean(i,1) Theta_mean(i,3); Theta_mean(i,3) Theta_mean(i,2)];

        % Electron temperature
        Te = Ti*10^Log10_Te_Ti(i);

        % Transform grid in the transformed coordinates to the original coordinates
        [xNU, yNU, xNUpolyarray, yNUpolyarray] = transform_gridpolygonsTC(xNU_ref, yNU_ref, xNUpolyarray_ref, yNUpolyarray_ref, thetabar, Ebar, Theta);

        % Evaluate IEAD from IEAD*
        [IEAD_NU] = evaluate_IEAD_from_IEADstar(xNU(:),yNU(:),thetabar,Ebar,Theta,xstar(:),ystar(:),IEADstar(:));

        % Mean and standard devation of angle & energy
        [~, ~, expected_theta, expected_E, std_theta, std_E] = MP_Expected_Theta_E(xNU,yNU,xNUpolyarray,yNUpolyarray,IEAD_NU,Te);
        Moments_Theta_E(i,:) = [expected_theta expected_E std_theta std_E];

    end

    % ====================================================================
    %                    Plot Contour Moments
    % ====================================================================

    % Number of contour lines
    Nlines = 35;
    
    % Initialize figure counter
    nfig = 0;

    % ------------
    %  theta mean
    % ------------

    nfig = nfig + 1;

    figure(nfig)

    % Data
    data = Moments_Theta_E(:,1);
    moment_plot = reshape(data,Npoints,Npoints);
    xmesh_plot = reshape(Log10_Te_Ti,Npoints,Npoints);
    ymesh_plot = reshape(Psi,Npoints,Npoints);

    % Plot contours
    ztext = '$\mu_{\theta}$';
    title_text = '';
    plot_moment_contour(xmesh_plot,ymesh_plot,moment_plot,ztext,title_text,Nlines)

    set(gcf,'Position',[500 800 560 420])

    % ------------
    %    E mean
    % ------------

    nfig = nfig + 1;

    figure(nfig)

    % Data
    data = Moments_Theta_E(:,2);
    moment_plot = reshape(data,Npoints,Npoints);
    xmesh_plot = reshape(Log10_Te_Ti,Npoints,Npoints);
    ymesh_plot = reshape(Psi,Npoints,Npoints);

    % Plot contours
    ztext = '$\mu_{E}$';
    title_text = '';
    plot_moment_contour(xmesh_plot,ymesh_plot,moment_plot,ztext,title_text,Nlines)

    set(gcf,'Position',[1200 800 560 420])

    % ------------
    %  theta std.
    % ------------

    nfig = nfig + 1;

    figure(nfig)

    % Data
    data = Moments_Theta_E(:,3);
    moment_plot = reshape(data,Npoints,Npoints);
    xmesh_plot = reshape(Log10_Te_Ti,Npoints,Npoints);
    ymesh_plot = reshape(Psi,Npoints,Npoints);

    % Plot contours
    ztext = '$\sigma_{\theta}$';
    title_text = '';
    plot_moment_contour(xmesh_plot,ymesh_plot,moment_plot,ztext,title_text,Nlines)

    set(gcf,'Position',[500 200 560 420])

    % ------------
    %    E std.
    % ------------

    nfig = nfig + 1;

    figure(nfig)

    % Data
    data = Moments_Theta_E(:,4);
    moment_plot = reshape(data,Npoints,Npoints);
    xmesh_plot = reshape(Log10_Te_Ti,Npoints,Npoints);
    ymesh_plot = reshape(Psi,Npoints,Npoints);

    % Plot contours
    ztext = '$\sigma_{E}$';
    title_text = '';
    plot_moment_contour(xmesh_plot,ymesh_plot,moment_plot,ztext,title_text,Nlines)
    
    set(gcf,'Position',[1200 200 560 420])

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The function plot_moment_contour generates a contour plot for given 
% moment data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Input
% -----
% - xmesh_fit: x-coordinates of mesh
% - ymesh_fit: y-coordinates of mesh
% - fit: moment data at the mesh points
% - ztext: z-label
% - title_text: plot title
% - Nlines: number of contour lines

% Author: Pablo Seleson
% ------

% Last Modified: February 10, 2022
% -------------

function plot_moment_contour(xmesh_fit,ymesh_fit,fit,ztext,title_text,Nlines)

    % Plot contours
    contourf(xmesh_fit,ymesh_fit,fit,Nlines)

    % Labels
    set(gca,'fontsize', 20)
    xlabel('$\log_{10}(T_e/T_i)$','Interpreter','latex','Fontsize',24)
    ylabel('$\Psi$','Interpreter','latex','Fontsize',24)
    zlabel(ztext,'Interpreter','latex','Fontsize',24)
    box('on')
    title(title_text,'Interpreter','latex','FontSize',26)

    % Colorbar
    hcb = colorbar;

    % Colorbar title
    colorTitleHandle = get(hcb,'Title');
    titleString = ztext;
    set(colorTitleHandle ,'String',titleString,'interpreter','latex','FontSize',24)
    
end
