
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The function main runs the sequence of codes needed to create (train) and 
% test the surrogate ion energy-angle distribution (IEAD) models as well as 
% to compute the sensitivity indices
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Discussion
% ----------

% Tasmanian directory: 
%       - This directory is where the tasmanian sparse grids containing the 
%       moments and IEAD surrogate models are stored.
%       - The main function uses this directory to check whether a 
%       surrogate model already exists.

% Data: 
%       The data needs to be stored in the directories defined under the 
%       section "Data Directories". The data contains:
%       (1) .log files with information about the different cases 
%           (each case represents a different combination of physical
%           parameters). These files should be placed directly in the
%           chosen directories.
%       (2) particle data which should be placed within a subdirectory 
%           "data/" inside the chosen directories.

% Functions: 
%       The functions are as follows:
%       - "MomentsSurrogate": Creates the moments surrogate models.
%       - "IEADstart_data_computation": Computes (training or testing) 
%          data in the form of IEAD and IEAD* (IEAD in the transformed 
%          coordinates).
%       - "IEADstarSurrogate": Creates the IEAD surrogate models.
%       - "IEADstarSurrogateError": Computes IEAD surrogate (training 
%          or testing) errors.
%       - "plot_IEADstar_model_vs_data": Plots comparisons between the 
%         IEAD surrogate models and the (training or testing) IEAD data.
%       - "GSA_Si": Computes sensitivity indices.
%       - "GSA_Moments_2D": Computes angle and energy moments profiles 
%         (mean and standard deviation) using the 2D IEAD surrogate models
%         (only for 2D).

% Note: 
%       The various inputs for the functions are as follows:
%       - sglevel: (training) data sparse grid level: "7", "10", or "13"
%       - dim: "2" (for 2D) or "4" (for 4D)
%       - Directory: "Directory_train" for training data
%                    "Directory_test" for testing data
%       - flag_plot: flag to produce plots
%       - data_type: 'train' for training and 'test' for testing 
%       - N: number of samples 
%       - Npoints: number of samples per dimension

% Input:
%       The code runs by default for sglevel = 7 & dim = 2. 
%       These values are hardcoded but can be modified under "Input".

% Author: Pablo Seleson
% ------

% Last Modified: March 10, 2022
% -------------

function main

    clear all
    close all
    clc

    % ====================================================================
    %                              Input
    % ====================================================================

    % Parameter space dimension
    dim = 2;

    % Data sparse grid level
    sglevel = 7;

    % Flag to plot   
    flag_plot = 1;
    
    % Check dimension
    if dim ~= 2 && dim ~= 4
        error('dim should be 2 or 4.')
    end

    % Sparse grid level
    if sglevel ~= 7 && sglevel ~= 10 && sglevel ~= 13
        error('sglevel should be 7, 10, or 13.')
    end

    fprintf(' ======================================================== \n')
    fprintf('             %gD hPIC IEAD Surrogate Model \n',dim)
    fprintf('                Sparse Grid Level: %g   \n',sglevel)
    fprintf(' ======================================================== \n\n')

    % ====================================================================
    %                       Tasmanian Directory
    % ====================================================================

    Directory_tsg = '';

    % ====================================================================
    %                       Data Directories
    % ====================================================================

    % Training data directory
    if dim == 2
        % Directory of data for 2D case
        Directory_train = '';
    else
        % Directory of data for 4D case
        Directory_train = '';
    end

    % Testing data directory
    if dim == 2
        % Directory of data for 2D case
        Directory_test = '';
    else
        % Directory of data for 4D case
        Directory_test = '';
    end

    % ====================================================================
    %                    Create data folders
    % ====================================================================

    if dim == 2
        if ~exist('IEAD_data_2D', 'dir')
            mkdir('IEAD_data_2D')
        end
    else
        if ~exist('IEAD_data_4D', 'dir')
            mkdir('IEAD_data_4D')
        end
    end

    % ====================================================================
    %                  Create moments surrogate model
    % ====================================================================

    fprintf(' ---------------- Moments surrogate model --------------- \n\n')

    % Check if moments surrogate model already exists
    Moments_surrogate_filename = [Directory_tsg 'LS_' num2str(dim) 'D_Grid_Moments_level_' num2str(sglevel) '_FileG']; 

    % Initialize default answer
    str = 'N';

    if isfile(Moments_surrogate_filename)
        % File exists
        fprintf('Moments surrogate model for %gD level %g exists: \n',dim,sglevel)
        prompt = 'Build the moments surrogate model again? [Y/N]: ';
        str = input(prompt,'s');
        fprintf('\n')
    end

    if strcmp( str, 'Y' ) || ~isfile(Moments_surrogate_filename)
        tic
        fprintf('Building moments surrogate model ... \n\n')

        % Run function
        MomentsSurrogate(sglevel,dim,Directory_train,flag_plot)
        
        fprintf('Moments surrogate model built: time = %f (sec)  \n\n', toc);
    end

    % ====================================================================
    %                   Create IEAD training data
    % ====================================================================

    fprintf(' ----------------   IEAD training data    --------------- \n\n')

    % Data type
    data_type = 'train';

    % Check if IEAD training data already exists
    IEAD_traindata_filename = ['IEAD_data_' num2str(dim) 'D/Level_' num2str(sglevel) '_' data_type '.mat'];

    % Initialize default answer
    str = 'N';

    if isfile(IEAD_traindata_filename)
        % File exists
        fprintf('IEAD training data for %gD level %g exists: \n',dim,sglevel)
        prompt = 'Create IEAD training data again? [Y/N]: ';
        str = input(prompt,'s');
        fprintf('\n')
    end

    if strcmp( str, 'Y' ) || ~isfile(IEAD_traindata_filename)
        tic
        fprintf('Generating IEAD training data ... \n\n')

        % Run function
        IEADstart_data_computation(sglevel,dim,Directory_train,data_type)

        fprintf('\nIEAD training data generated: time = %f (sec)  \n\n', toc);
    end

    % ====================================================================
    %                 Generate IEAD surrogate model
    % ====================================================================

    fprintf(' ----------------  IEAD surrogate model   --------------- \n\n')

    % Check if IEAD* surrogate model already exists
    IEADstar_surrogate_filename = [Directory_tsg 'LS_' num2str(dim) 'D_Grid_IEADstar_level_' num2str(sglevel) '_FileG']; 

   % Initialize default answer
    str = 'N';

    if isfile(IEADstar_surrogate_filename)
        % File exists
        fprintf('IEAD surrogate model for %gD level %g exists: \n',dim,sglevel)
        prompt = 'Build the IEAD surrogate model again? [Y/N]: ';
        str = input(prompt,'s');
        fprintf('\n')
    end

    if strcmp( str, 'Y' ) || ~isfile(IEADstar_surrogate_filename)
        tic
        fprintf('Building IEAD surrogate model ... \n\n')

        % Run function
        IEADstarSurrogate(sglevel,dim)

        fprintf('IEAD surrogate model built: time = %f (sec)  \n\n', toc);
    end

    % ====================================================================
    %              Compute and plot IEAD training error
    % ====================================================================

    fprintf(' ----------------  IEAD training error    --------------- \n\n')

    prompt = 'Would you like to compute the IEAD training error? [Y/N]: ';
    str = input(prompt,'s');
    fprintf('\n')

    if strcmp( str, 'Y' )
        tic
        fprintf('Computing IEAD training error ... \n\n')

        % Run function
        data_type = 'train';
        IEADstarSurrogateError(sglevel,dim,data_type)

        fprintf('IEAD training error computed: time = %f (sec)  \n\n', toc);
    end

    % ===================================================================
    %           Compare IEAD model with the training data 
    % ===================================================================

    fprintf(' ------------  IEAD model vs. (training) data ----------- \n\n')

    prompt = 'Would you like to compare the IEAD model with the (training) data in the transformed coordinates? [Y/N]: ';
    str = input(prompt,'s');
    fprintf('\n')

    if strcmp( str, 'Y' )
        close all
        data_type = 'train';
        % IEAD*
        plot_IEADstar_model_vs_data(sglevel,dim,data_type)
    end

    prompt = 'Would you like to compare the IEAD model with the (training) data in the original coordinates? [Y/N]: ';
    str = input(prompt,'s');
    fprintf('\n')

    if strcmp( str, 'Y' )
        close all
        data_type = 'train';
        % IEAD
        plot_IEAD_model_vs_data(sglevel,dim,data_type)
    end

    % ===================================================================
    %                   Create IEAD testing data
    % ===================================================================
 
    fprintf(' ----------------   IEAD testing data     --------------- \n\n')

    % Data type
    data_type = 'test';

    % Check if IEAD testing data already exists
    IEAD_testdata_filename = ['IEAD_data_' num2str(dim) 'D/Level_' num2str(sglevel) '_' data_type '.mat'];

    % Initialize default answer
    str = 'N';

    if isfile(IEAD_testdata_filename)
        % File exists
        fprintf('IEAD testing data for %gD level %g exists: \n',dim,sglevel)
        prompt = 'Create IEAD testing data again? [Y/N]: ';
        str = input(prompt,'s');
        fprintf('\n')
    end

    if strcmp( str, 'Y' ) || ~isfile(IEAD_testdata_filename)
        tic
        fprintf('Generating IEAD testing data ... \n\n')

        % Run function
        IEADstart_data_computation(sglevel,dim,Directory_test,data_type)

        fprintf('\nIEAD testing data generated: time = %f (sec)  \n\n', toc);
    end

    % ====================================================================
    %              Compute and plot IEAD testing error
    % ====================================================================

    fprintf(' ----------------   IEAD testing error    --------------- \n\n')

    prompt = 'Would you like to compute the IEAD testing error? [Y/N]: ';
    str = input(prompt,'s');
    fprintf('\n')

    if strcmp( str, 'Y' )
        tic
        fprintf('Computing IEAD testing error ... \n\n')

        % Run function
        data_type = 'test';
        IEADstarSurrogateError(sglevel,dim,data_type)

        fprintf('IEAD testing error computed: time = %f (sec)  \n\n', toc);
    end

    % ====================================================================
    %           Compare IEAD model with the testing data
    % ====================================================================

    fprintf(' ------------  IEAD model vs. (testing) data  ----------- \n\n')

    prompt = 'Would you like to compare the IEAD model with the (testing) data in the transformed coordinates? [Y/N]: ';
    str = input(prompt,'s');
    fprintf('\n')

    if strcmp( str, 'Y' )
        close all
        data_type = 'test';
        % IEAD*
        plot_IEADstar_model_vs_data(sglevel,dim,data_type)
    end

    prompt = 'Would you like to compare the IEAD model with the (testing) data in the original coordinates? [Y/N]: ';
    str = input(prompt,'s');
    fprintf('\n')

    if strcmp( str, 'Y' )
        close all
        data_type = 'test';
        % IEAD
        plot_IEAD_model_vs_data(sglevel,dim,data_type)
    end

    % ===================================================================
    %               Compute sensitivity indices
    % =================================================================== 

    fprintf(' ---------------   Sensitivity indices    --------------- \n\n')

    prompt = 'Would you like to compute sensitivity indices? [Y/N]: ';
    str = input(prompt,'s');
    fprintf('\n')

    if strcmp( str, 'Y' )
        close all

        prompt = 'Enter number of samples N (default 50000): ';
        str = input(prompt,'s');
        fprintf('\n')

        % Number of samples
        N = str2num(str);

        tic
        fprintf('Computing IEAD sensitivity indices ... \n\n')

        % Run function
        GSA_Si(sglevel,dim,N)

        fprintf('\nSensitivity indices computed: time = %f (sec)  \n\n', toc);

    end

    % ===================================================================
    %           Compute angle and energy moments profiles 
    % =================================================================== 

    if dim == 2

        fprintf(' --------------- Angle and energy moments --------------- \n\n')

        prompt = 'Would you like to compute angle and energy moments profiles? [Y/N]: ';
        str = input(prompt,'s');
        fprintf('\n')

        if strcmp( str, 'Y' )
            close all

            prompt = 'Enter number of samples per dimension Npoints (default 250): ';
            str = input(prompt,'s');
            fprintf('\n')

            % Number of samples per dimension in 2D
            Npoints = str2num(str);

            tic
            fprintf('Computing angle and energy moments ... \n\n')

            % Run function
            GSA_Moments_2D(sglevel, Npoints)

            fprintf('\nAngle and energy moments computed: time = %f (sec)  \n\n', toc);

        end

    end

end
