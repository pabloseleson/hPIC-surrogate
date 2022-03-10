
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The function IEADstart_data_computation computes and saves data for the
% IEAD and IEAD* (IEAD in the transformed coordinates). The 
% transformation uses analytical functions to approximate the data moments. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Input
% -----
% sglevel: data sparse grid level
% dim: dimension of parameter space
%      2D: dim = 2
%      4D: dim = 4
% Directory: data directory
% data_type: 'train' for training data and 'test' for testing data

% NOTE: The Directory should contain:
%       (1) the cases log file: for example, "cases_nD_2_level_7.log"
%           for the level 7 data sparse grid in 2D
%       (2) particle data inside subfolder "data/" given by:
%           [filename '_WallParticleList_sp0.dat']
%           for different file names

% Author: Pablo Seleson
% ------

% Last Modified: March 9, 2022
% -------------

function IEADstart_data_computation(sglevel,dim,Directory,data_type)

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

    if strcmp(data_type,'train')
        % Filename of cases data: training data
        cases_file = [Directory 'cases_nD_' num2str(dim) '_level_' num2str(sglevel) '.log'];

    elseif strcmp(data_type,'test')
        if dim == 2
            % Filename of 2D cases data: testing data 
            cases_file = [Directory 'uniform_grid_2D_20.log'];
        else
            % Filename of 4D cases data: testing data 
            cases_file = [Directory 'uniform_grid_4D_5.log'];
        end
    else
        error('data_type uknown')
    end

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
       
    % ====================================================================
    %               Create Meshes for IEAD and IEAD*
    % ====================================================================
    
    % ------------------
    %   Original grid
    % ------------------
    
    % Number of points in each dimension
    Ntheta = 90;
    NE = 240; 
    
    % Array of theta values
    XEDGES = linspace(0, pi*90/180, Ntheta+1);
    XNODES = (XEDGES(1:end-1) + XEDGES(2:end))/2;
    dimX = length(XNODES);
    
    % Array of energy values
    YEDGES = linspace(0, 24, NE+1); % Max energy is 24*Te
    YNODES = (YEDGES(1:end-1) + YEDGES(2:end))/2;
    dimY = length(YNODES);
    
    % ------------------
    %  Transformed grid
    % ------------------
    
    N = 25; % Number of points in [-1,1];
    xL = 5; % Number of standard deviations
    
    % Array of theta* values
    [XNODESstar] = normaldistgrid(N,xL);
    dimXstar = length(XNODESstar);
    
    % Array of E* values
    YNODESstar = XNODESstar;
    dimYstar = length(YNODESstar);
    
    % ====================================================================
    %            Compute Moments per Case for Transformation
    %               based on the Moments Surrogate Model
    % ====================================================================
    
    % Reload moments grid
    gridname = ['LS_' num2str(dim) 'D_Grid_Moments_level_' num2str(sglevel)];
    [lGrid_Moments] = tsgReloadGrid(gridname);
    
    % Evaluate fitted model on data points
    if dim == 2
        % 2D case
        [result] = tsgEvaluate(lGrid_Moments, [Log10_Te_Ti Psi]);
    else
        % 4D case
        [result] = tsgEvaluate(lGrid_Moments, [Log10_Te_Ti Psi B Log10_n_values]);
    end
    
    % Assign moments to arrays
	thetabar_mean = result(:,1);
	Ebar_mean = 10.^result(:,2);
    Theta_mean_11 = 10.^result(:,3);
    Theta_mean_22 = 10.^result(:,4);
    Theta_mean_12 = result(:,5);
	Theta_mean = [Theta_mean_11 Theta_mean_22 Theta_mean_12]; % [Theta(1,1) Theta(2,2) Theta(1,2)];
    
    % ====================================================================
    %             Read Data and Compute IEAD and IEAD*
    % ====================================================================
        
    % Initialize array of number of runs per case
    nruns_array = zeros(ncases,1);
    
    % Initialize array of sparse grid points
    points = zeros(ncases,dim);
    
    % Initialize IEAD and IEAD* arrays  
    dim_bins = dimXstar * dimYstar; % Number of bins
    IEAD_array = zeros(ncases,dimX*dimY);
    IEADstar_array = zeros(ncases,dim_bins);
    
    % Initialize moments arrays
    thetabar_array = zeros(ncases,1);
    Ebar_array = zeros(ncases,1);
    Theta_array = zeros(ncases,3);
    
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

            % Data
            if dim == 2

                % -------
                % 2D case 
                % -------

                % Find temperature ratio and magnetic field angle
                Log10_Te_Ti_value = Log10_Te_Ti(n);
                Psi_value = Psi(n);

                % Assign data
                if rep_value == 1
                    points(ncase,1) = Log10_Te_Ti_value;
                    points(ncase,2) = Psi_value;

                    thetabar_array(ncase) = thetabar_mean(n);
                    Ebar_array(ncase) = Ebar_mean(n);
                    Theta_array(ncase,:) = [Theta_mean(n,1) Theta_mean(n,2) Theta_mean(n,3)];

                else
                    if abs(points(ncase,1) - Log10_Te_Ti_value) > 1e-14 || abs(points(ncase,2) - Psi_value) > 1e-14
                        error('Case with different physical parameters.')
                    end
                end

            else

                % -------
                % 4D case
                % -------

                % Find temperature ratio, magnetic field angle, magnetic field magnitude, and density
                Log10_Te_Ti_value = Log10_Te_Ti(n);
                Psi_value = Psi(n);
                B_value = B(n);
                Log10_n_value = Log10_n_values(n);

                % Assign data
                if rep_value == 1
                    points(ncase,1) = Log10_Te_Ti_value;
                    points(ncase,2) = Psi_value;
                    points(ncase,3) = B_value;
                    points(ncase,4) = Log10_n_value;

                    thetabar_array(ncase) = thetabar_mean(n);
                    Ebar_array(ncase) = Ebar_mean(n);
                    Theta_array(ncase,:) = [Theta_mean(n,1) Theta_mean(n,2) Theta_mean(n,3)];

                else
                    if abs(points(ncase,1) - Log10_Te_Ti_value) > 1e-14 || abs(points(ncase,2) - Psi_value) > 1e-14 || abs(points(ncase,3) - B_value) > 1e-14 || abs(points(ncase,4) - Log10_n_value) > 1e-14
                        error('Case with different physical parameters.')
                    end
                end

            end

            % Compute number of runs per case
            nruns_array(ncase) = nruns_array(ncase) + 1;

            if rep_value == 1
                fprintf('\n case: %4g  of  %4g - nrun = %g',ncase,ncases,rep_value)
            else
                fprintf(' %g',rep_value)
            end

            % ------------------------------------------------------------
            %                  Compute E and theta
            % ------------------------------------------------------------

            % Read particle energies and angles
            [E_values,theta_values] = read_E_theta(data_file);

            % Rescale energy by electron temperature
            Te_value = Te(n);
            E_values = E_values/Te_value;

            % ------------------------------------------------------------
            %                     Compute IEAD
            % ------------------------------------------------------------

            % Distribution domain limits
            xmin = min(XEDGES);
            xmax = max(XEDGES);
            ymin = min(YEDGES);
            ymax = max(YEDGES);

            [~,~,f,~] = compute_2D_NU_distribution(theta_values,E_values,XNODES,YNODES,xmin,xmax,ymin,ymax);

            % Normalize distribution
            dim_data = length(theta_values);
            f = f/dim_data;

            % Compute IEAD mean
            data = reshape(f,1,dimX*dimY);
            IEAD_array(ncase,:) = IEAD_array(ncase,:) + data;

            % ------------------------------------------------------------
            %                 Transform theta and E
            % ------------------------------------------------------------

            % Assign moments from fit
            thetabar_mean_value = thetabar_array(ncase);
            Ebar_mean_value = Ebar_array(ncase);
            Theta_mean_value = [Theta_array(ncase,1) Theta_array(ncase,3); Theta_array(ncase,3) Theta_array(ncase,2)];

            % Perform transformation
            [thetastar_values,Estar_values] = transform_E_theta(theta_values,E_values,thetabar_mean_value,Ebar_mean_value,Theta_mean_value);

            % ------------------------------------------------------------
            %                     Compute IEAD*
            % ------------------------------------------------------------

            % Distribution domain limits
            xstarmin = -5;
            xstarmax =  5;
            ystarmin = -5;
            ystarmax =  5;

            [~,~,fstar,~] = compute_2D_NU_distribution(thetastar_values,Estar_values,XNODESstar,YNODESstar,xstarmin,xstarmax,ystarmin,ystarmax);

            % Normalize distribution
            dim_datastar = length(thetastar_values);
            fstar = fstar/dim_datastar;

            % Computations for IEAD* mean and standard deviation
            datastar = reshape(fstar,1,dim_bins);
            IEADstar_array(ncase,:) = IEADstar_array(ncase,:) + datastar;

        end

        % Close file
        fclose(fid);

    end

    % Finalize averaging
    for ncase = 1:ncases    
        % IEAD*
        IEADstar_array(ncase,:) = IEADstar_array(ncase,:)/nruns_array(ncase);

        % IEAD      
        IEAD_array(ncase,:) = IEAD_array(ncase,:)/nruns_array(ncase);
    end
  
    % ====================================================================
    %                          Save Data 
    % ====================================================================

    IEAD_data_filename = ['IEAD_data_' num2str(dim) 'D/Level_' num2str(sglevel) '_' data_type '.mat'];
    save(IEAD_data_filename,'points','XNODES','YNODES','IEAD_array','XNODESstar','YNODESstar','IEADstar_array')
    
end
