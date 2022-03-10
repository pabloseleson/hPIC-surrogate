
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The function MomentsSurrogate computes and fits the moments (angle & 
% energy means and convariance tensor) to transform energies and angles
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Input
% -----
% sglevel: data sparse grid level
% dim: dimension of parameter space
%      2D: dim = 2
%      4D: dim = 4
% Directory: data directory
% flag_plot: if flag_plot == 1 & dim == 2, then plots moments for the 2D case

% NOTE: The Directory should contain:
%       (1) the cases log file: for example, "cases_nD_2_level_7.log"
%           for the level 7 data sparse grid in 2D
%       (2) particle data inside subfolder "data/" given by:
%           [filename '_WallParticleList_sp0.dat']
%           for different file names

% Author: Pablo Seleson
% ------

% Last Modified: February 1, 2022
% -------------

function MomentsSurrogate(sglevel,dim,Directory,flag_plot)

    % Check sglevel input
    if sglevel~= 7 && sglevel~= 10 && sglevel~= 13
        error('sglevel should be 7, 10, or 13.')
    end

    % Check dim input
    if dim~= 2 && dim~= 4
        error('dim should be 2 or 4.')
    end

    % ====================================================================
    %                          Read Cases
    % ====================================================================
     
    % Filename of cases data
    cases_file = [Directory 'cases_nD_' num2str(dim) '_level_' num2str(sglevel) '.log'];

    % Read cases
    if dim == 2
        % Read data for 2D case
        [filename,case_number,rep_number,Log10_Te_Ti,Te,~,Psi,B,~] = read_2D_cases(cases_file);
    else
        % Read data for 4D case
        [filename,case_number,rep_number,Log10_Te_Ti,Te,~,Psi,B,Log10_n_values,~] = read_4D_cases(cases_file);
    end

    % Number to cases
    ncases = max(case_number) + 1;
    
    % Number of data
    ndata = length(filename);
    
    % Initialize points array
    points = zeros(ncases,dim);
    
    % ====================================================================
    %      Read Data to Compute Moments per Case for Transformation
    % ====================================================================
    
    % Initialize moments
    thetabar_mean = zeros(ncases,1);
    Ebar_mean = zeros(ncases,1);
    Theta_mean = zeros(ncases,3); % [Theta(1,1) Theta(2,2) Theta(1,2)];
    
    % Initialize array of number of runs per case
    nruns_array = zeros(ncases,1);
    
    % Run over data
    for n = 1:ndata

        % Data file
        data_file = [Directory 'data/' char(filename(n)) '_WallParticleList_sp0.dat'];

        % -----------------------------
        %    Check if file is empty
        % -----------------------------

        % Open file
        fid = fopen(data_file);

        % Check if file is empty
        if all(fgetl(fid) == -1)

            % -------------
            % file IS empty
            % -------------

            % Read case and repetition numbers
            ncase = case_number(n) + 1;
            rep_value = rep_number(n);

            % Print empty file information
            fprintf('File is empty: ncase = %g ; rep_value = %g \n',ncase,rep_value)
            fprintf('File name: %s \n',char(filename(n)))

        else

            % -----------------
            % file IS NOT empty
            % -----------------

            % Read case and repetition numbers
            ncase = case_number(n) + 1;
            rep_value = rep_number(n);

            if rep_value == 1
                fprintf('case: %4g  of  %4g \n',ncase,ncases)

                if dim == 2
                    % 2D case
                    points(ncase,1) = Log10_Te_Ti(n);
                    points(ncase,2) = Psi(n);
                else
                    % 4D case
                    points(ncase,1) = Log10_Te_Ti(n);
                    points(ncase,2) = Psi(n);
                    points(ncase,3) = B(n);
                    points(ncase,4) = Log10_n_values(n);
                end
            end

            % --------------------------------------------------
            %  Read particle velocities and compute E and theta
            % --------------------------------------------------

            % Read particle energies and angles
            [E_values,theta_values] = read_E_theta(data_file);

            % Rescale energy by the electron temperature
            Te_value = Te(n);
            E_values = E_values/Te_value;

            % -------------------------------------------------
            %               Compute Moments
            % -------------------------------------------------

            [thetabar,Ebar,Theta] = moments_E_theta(theta_values,E_values);

            % Compute average thetabar
            thetabar_mean(ncase) = thetabar_mean(ncase) + thetabar;

            % Compute average Ebar
            Ebar_mean(ncase) = Ebar_mean(ncase) + Ebar;

            % Compute average Theta
            Theta_mean(ncase,:) = Theta_mean(ncase,:) + [Theta(1,1) Theta(2,2) Theta(1,2)];

            % Compute number of runs per case
            nruns_array(ncase) = nruns_array(ncase) + 1;

        end

        % Close file
        fclose(fid);

    end
    
    % Finalize averaging
    thetabar_mean = thetabar_mean./nruns_array;
    Ebar_mean = Ebar_mean./nruns_array;
    Theta_mean = Theta_mean./repmat(nruns_array,1,3);

    % ====================================================================
    %                         Fit Moments
    % ====================================================================

    % --------------------------------------------
    % Create Basis Functions: Lagrange Polynomials
    % --------------------------------------------

    iDepth_basis = 4;
    sType = 'iptotal';
    iOut = 5; % 5 moment elements

    gridname = ['LS_' num2str(dim) 'D_Grid_Moments_level_' num2str(sglevel)];

    if dim == 2
        % 2D case
        [lGrid_Moments, ~] = tsgMakeGlobal(gridname,2,iOut,'clenshaw-curtis',sType,iDepth_basis,[min(points(:,1)) max(points(:,1)); min(points(:,2)) max(points(:,2))]);
    else
        % 4D case
        [lGrid_Moments, ~] = tsgMakeGlobal(gridname,4,iOut,'clenshaw-curtis',sType,iDepth_basis,[min(points(:,1)) max(points(:,1)); min(points(:,2)) max(points(:,2)); min(points(:,3)) max(points(:,3)); min(points(:,4)) max(points(:,4))]);
    end

    % --------------------------------------------
    %          Find Coefficients
    % --------------------------------------------
    
    % Evaluate basis functions at grid points
    % Note: B has dimensions [Npoints x Nbasis]
    % ----
    if dim == 2
        % 2D case
        B = tsgEvaluateHierarchy(lGrid_Moments,[points(:,1) points(:,2)]);
    else
        % 4D case
        B = tsgEvaluateHierarchy(lGrid_Moments,[points(:,1) points(:,2) points(:,3) points(:,4)]);
    end
    
    % Solve coefficients
    % Note: rhs has dimensions [Npoints x iOut]
    % ----
    % Two cases:
    % (1) If A is an N-by-N matrix and B is a column vector with N components, 
    %     or a matrix with several such columns, then X = A\B is the solution 
    %     to the equation A*X = B.
    % (2) If A is an M-by-N matrix with M < or > N and B is a column vector 
    %     with M components, or a matrix with several such columns, then 
    %     X = A\B is the solution in the least-squares sense to the under- 
    %     or over-determined system of equations A*X = B.
    
    rhs = [thetabar_mean log10(Ebar_mean) log10(Theta_mean(:,1)) log10(Theta_mean(:,2))  Theta_mean(:,3)];
    a = B\rhs; % This solves B*a = rhs 
    
    % Load coefficients to the grid object
    tsgLoadHCoefficients(lGrid_Moments, a);
    
    % ====================================================================
    %                    Compute Fitting Error
    % ====================================================================
    
    % Evaluate fitted model on data points
    if dim == 2
        % 2D case
        [result_fit] = tsgEvaluate(lGrid_Moments, [points(:,1) points(:,2)]); 
    else
        % 4D case
        [result_fit] = tsgEvaluate(lGrid_Moments, [points(:,1) points(:,2) points(:,3) points(:,4)]); 
    end

    % Read fit 
    thetabar_mean_fit = result_fit(:,1);
    Ebar_mean_fit = 10.^result_fit(:,2);
    Theta11_mean_fit = 10.^result_fit(:,3);
    Theta22_mean_fit = 10.^result_fit(:,4);
    Theta12_mean_fit = result_fit(:,5);
    
    Rel_Error_thetabar = norm(thetabar_mean - thetabar_mean_fit,2)/norm(thetabar_mean,2);
    Rel_Error_Ebar = norm(Ebar_mean - Ebar_mean_fit,2)/norm(Ebar_mean,2);
    Rel_Error_Theta11 = norm(Theta_mean(:,1) - Theta11_mean_fit,2)/norm(Theta_mean(:,1),2);
    Rel_Error_Theta22 = norm(Theta_mean(:,2) - Theta22_mean_fit,2)/norm(Theta_mean(:,2),2);
    Rel_Error_Theta12 = norm(Theta_mean(:,3) - Theta12_mean_fit,2)/norm(Theta_mean(:,3),2);

    fprintf('\nRelative errors: \n')
    fprintf('Mean theta: %f \n',Rel_Error_thetabar)
    fprintf('Mean E    : %f \n',Rel_Error_Ebar)
    fprintf('Theta_{11}: %f \n',Rel_Error_Theta11)    
    fprintf('Theta_{22}: %f \n',Rel_Error_Theta22)   
    fprintf('Theta_{12}: %f \n\n',Rel_Error_Theta12)   

    % ====================================================================
    %             Compute Covariance Matrix Determinant
    % ====================================================================

    det_Theta = Theta11_mean_fit.*Theta22_mean_fit - Theta12_mean_fit.^2;

    if all(det_Theta > 0)
        fprintf('All determinant values are positive! \n\n')

    else
        find(det_Theta < 1e-14)
    end

    % ====================================================================
    %                      Plot Fitted Moments
    % ====================================================================
    
    if flag_plot == 1 && dim == 2

        % ----------------------------------------------
        %    Create grid to evaluate fitted model
        % ----------------------------------------------
        NLog10_Te_Ti = 100;
        NPsi = 100;

        % Array of Log10_Te_Ti values
        Log10_Te_Ti_EDGES = linspace(min(points(:,1)), max(points(:,1)), NLog10_Te_Ti+1);
        Log10_Te_Ti_NODES = (Log10_Te_Ti_EDGES(1:end-1) + Log10_Te_Ti_EDGES(2:end))/2;
        dimTe_Ti = length(Log10_Te_Ti_NODES);

        % Array of angle values
        Psi_EDGES = linspace(min(points(:,2)), max(points(:,2)), NPsi+1);
        Psi_NODES = (Psi_EDGES(1:end-1) + Psi_EDGES(2:end))/2;
        dimPsi = length(Psi_NODES);

        % Grid arrays
        Log10_Te_Ti_array = zeros(dimTe_Ti*dimPsi,1);
        Psi_array = zeros(dimTe_Ti*dimPsi,1);

        count = 0;
        for i = 1:dimTe_Ti
            for j = 1:dimPsi
                count = count + 1;
                Log10_Te_Ti_array(count) = Log10_Te_Ti_NODES(i);
                Psi_array(count) = Psi_NODES(j);
            end
        end

        % Create mesh to plot
        Log10_Te_Ti_mesh = reshape(Log10_Te_Ti_array,dimTe_Ti,dimPsi);
        Psi_mesh = reshape(Psi_array,dimTe_Ti,dimPsi);

        % ----------------------------------------------
        % Evaluate fitted model on reconstruction points
        % ----------------------------------------------
        [result] = tsgEvaluate(lGrid_Moments, [Log10_Te_Ti_array Psi_array]);

        % ----------------------------------------------
        %              Plot Moments 
        % ----------------------------------------------

        % Limits of parameter space
        limits = [min(points(:,1)) max(points(:,1)) min(points(:,2)) max(points(:,2))];

        % Initialize figure number
        nfig = 0;

        % ---------------
        %  Plot thetabar
        % ---------------

        % Data
        data = thetabar_mean;

        % Read and mesh fit
        thetabar_mean_fit = result(:,1);
        fit = reshape(thetabar_mean_fit,dimTe_Ti,dimPsi);

        % Plot
        nfig = nfig + 1;
        title_text = ['Rel. Err. = ' num2str(Rel_Error_thetabar,2)];
        plot_moment(nfig,points,data,Log10_Te_Ti_mesh,Psi_mesh,fit,'$\bar{\theta}$',title_text,limits)
        
        % Position figure
        set(gcf, 'Position',  [500, 800, 600, 500])

        % ---------------
        %    Plot Ebar
        % ---------------

        % Data
        data = Ebar_mean;

        % Read and mesh fit
        Ebar_mean_fit = 10.^result(:,2);
        fit = reshape(Ebar_mean_fit,dimTe_Ti,dimPsi);

        % Plot
        nfig = nfig + 1;
        title_text = ['Rel. Err. = ' num2str(Rel_Error_Ebar,2)];
        plot_moment(nfig,points,data,Log10_Te_Ti_mesh,Psi_mesh,fit,'$\bar{E}/T_e$',title_text,limits)
        
        % Position figure
        set(gcf, 'Position',  [1300, 800, 600, 500])

        % ---------------
        % Plot Theta_{11}
        % ---------------

        % Data
        data = Theta_mean(:,1);

        % Read and mesh fit
        Theta11_mean_fit = 10.^result(:,3);
        fit = reshape(Theta11_mean_fit,dimTe_Ti,dimPsi);

        % Plot
        nfig = nfig + 1;
        title_text = ['Rel. Err. = ' num2str(Rel_Error_Theta11,2)];
        plot_moment(nfig,points,data,Log10_Te_Ti_mesh,Psi_mesh,fit,'$\Theta_{11}$',title_text,limits)
        
        % Position figure
        set(gcf, 'Position',  [200, 200, 600, 500])

        % ---------------
        % Plot Theta_{22}
        % ---------------

        % Data
        data = Theta_mean(:,2);

        % Read and mesh fit
        Theta22_mean_fit = 10.^result(:,4);
        fit = reshape(Theta22_mean_fit,dimTe_Ti,dimPsi);

        % Plot
        nfig = nfig + 1;
        title_text = ['Rel. Err. = ' num2str(Rel_Error_Theta22,2)];
        plot_moment(nfig,points,data,Log10_Te_Ti_mesh,Psi_mesh,fit,'$\Theta_{22}$',title_text,limits)

        % Position figure
        set(gcf, 'Position',  [900, 200, 600, 500])

        % ---------------
        % Plot Theta_{12}
        % ---------------

        % Data
        data = Theta_mean(:,3);

        % Read and mesh fit
        Theta12_mean_fit = result(:,5);
        fit = reshape(Theta12_mean_fit,dimTe_Ti,dimPsi);

        % Plot
        nfig = nfig + 1;
        title_text = ['Rel. Err. = ' num2str(Rel_Error_Theta12,2)];
        plot_moment(nfig,points,data,Log10_Te_Ti_mesh,Psi_mesh,fit,'$\Theta_{12}$',title_text,limits)
        
        % Position figure
        set(gcf, 'Position',  [1600, 200, 600, 500])

    end
  
end

% ========================================================================
%                    Function to Plot Moments
% ========================================================================

% Input
% -----
% figid: figure id
% points_data: grid points in parameter space
% data: moment data
% xmesh_fit: x-coordinates of fitting mesh in parameter space
% ymesh_fit: y-coordinates of fitting mesh in parameter space
% fit: fitted moment
% ztext: zlabel 
% title_text: plot title
% limits: plot limits of parameters

% Author: Pablo Seleson
% ------

% Last Modified: February 1, 2022
% -------------

function plot_moment(figid,points_data,data,xmesh_fit,ymesh_fit,fit,ztext,title_text,limits)

    % Create figure
    figure_n = figure(figid);
    axes_n = axes('Parent',figure_n);
    hold(axes_n,'all');
    set(axes_n,'FontSize',20);
    
    % Colorbar
    colormap('parula');
    colorbar('FontSize',20);
    colorbar
    
    % Plot data
    scatter3(points_data(:,1),points_data(:,2),data,40,data,'o','filled');
    
    % Plot fit
    surf(xmesh_fit,ymesh_fit,fit,'FaceAlpha',0.5,'EdgeColor','none');
    view(-45,20)

    % Labels
    xlabel('$\log_{10}(T_e/T_i)$','Interpreter','latex','Fontsize',24)
    ylabel('$\Psi$','Interpreter','latex','Fontsize',24)
    zlabel(ztext,'Interpreter','latex','Fontsize',24)

    % Title
    title(title_text,'Interpreter','latex','FontSize',26)

    % Other settings
    axis(limits)
    box('on')
    pbaspect([1 1 1])
    
end
