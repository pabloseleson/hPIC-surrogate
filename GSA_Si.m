
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The function GSA_Si performs Global Sensitivity Analysis (GSA) based on 
% Variance-Based Methods, computing sensitivity indices following
%
% Saltelli, A., Ratto, M., Andres, T., Campolongo, F., Cariboni, J., 
% Gatelli, D., Saisana, M., Tarantola, S., 2008. Global Sensitivity Analysis. 
% The Primer. Wiley
%
% Chapter 4.6: How to Compute Sensitivity Indices (p. 164)
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Input
% -----
% sglevel: data sparse grid level
% dim: dimension of parameter space
%      2D: dim = 2
%      4D: dim = 4
% N: number of samples

% Author: Pablo Seleson
% ------

% Last Modified: February 9, 2022
% -------------

function GSA_Si(sglevel,dim,N)

    % Check sglevel input
    if sglevel~= 7 && sglevel~= 10 && sglevel~= 13
        error('sglevel should be 7, 10, or 13.')
    end

    % Check dim input
    if dim~= 2 && dim~= 4
        error('dim should be 2 or 4.')
    end

    % ====================================================================
    %                 Generate random matrices
    % ====================================================================

    if dim == 2
        
        % 2D case: grid limits
        min_Log10_Te_Ti = -1.25;
        max_Log10_Te_Ti = 1.25;
        
        min_Psi = 0;
        max_Psi = 87;

        % Quasi-random numbers
        Q = qrandstream('sobol',4);
        P = qrand(Q,N);

        % Arrays of physical parameters
        P(:,[1 3]) = min_Log10_Te_Ti + (max_Log10_Te_Ti-min_Log10_Te_Ti)*P(:,[1 3]);
        P(:,[2 4]) = min_Psi + (max_Psi-min_Psi)*P(:,[2 4]);

        % Matrices to compute sensitivity indices
        A = P(:,1:2);
        B = P(:,3:4);

        C1 = [A(:,1) B(:,2)];
        C2 = [B(:,1) A(:,2)];
        
    else
        
        % 4D case: grid limits
        min_Log10_Te_Ti = -1.25;
        max_Log10_Te_Ti = 1.25;

        min_Psi = 0;
        max_Psi = 87;

        min_B = 1;
        max_B = 15;

        min_Log10_n_values = 16;
        max_Log10_n_values = 20;

        % Quasi-random numbers
        Q = qrandstream('sobol',8);
        P = qrand(Q,N);

        % Arrays of physical parameters
        P(:,[1 5]) = min_Log10_Te_Ti + (max_Log10_Te_Ti-min_Log10_Te_Ti)*P(:,[1 5]);
        P(:,[2 6]) = min_Psi + (max_Psi-min_Psi)*P(:,[2 6]);
        P(:,[3 7]) = min_B + (max_B-min_B)*P(:,[3 7]);
        P(:,[4 8]) = min_Log10_n_values + (max_Log10_n_values-min_Log10_n_values)*P(:,[4 8]);

        % Matrices to compute sensitivity indices
        A = P(:,1:4);
        B = P(:,5:8);

        C1 = [A(:,1) B(:,2) B(:,3) B(:,4)];
        C2 = [B(:,1) A(:,2) B(:,3) B(:,4)];
        C3 = [B(:,1) B(:,2) A(:,3) B(:,4)];
        C4 = [B(:,1) B(:,2) B(:,3) A(:,4)];
        
    end

    fprintf('Random matrices generated \n')
    
    % ====================================================================
    %                    Evaluate IEAD* model
    % ====================================================================
    
    % Reload IEAD* fit grid
    gridname = ['LS_' num2str(dim) 'D_Grid_IEADstar_level_' num2str(sglevel)];
    [lGrid_IEADstar] = tsgReloadGrid(gridname);
    
    % Evaluate IEAD* surrogate model on points
    if dim == 2
        % Matrix A
        [IEADstar_all_A] = tsgEvaluate(lGrid_IEADstar, [A(:,1) A(:,2)]);
        % Matrix B
        [IEADstar_all_B] = tsgEvaluate(lGrid_IEADstar, [B(:,1) B(:,2)]);

        % Matrix C_1
        [IEADstar_all_C1] = tsgEvaluate(lGrid_IEADstar, [C1(:,1) C1(:,2)]);
        % Matrix C_2
        [IEADstar_all_C2] = tsgEvaluate(lGrid_IEADstar, [C2(:,1) C2(:,2)]);
    else
        % Matrix A
        [IEADstar_all_A] = tsgEvaluate(lGrid_IEADstar, [A(:,1) A(:,2) A(:,3) A(:,4)]);
        % Matrix B
        [IEADstar_all_B] = tsgEvaluate(lGrid_IEADstar, [B(:,1) B(:,2) B(:,3) B(:,4)]);

        % Matrix C_1
        [IEADstar_all_C1] = tsgEvaluate(lGrid_IEADstar, [C1(:,1) C1(:,2) C1(:,3) C1(:,4)]);
        % Matrix C_2
        [IEADstar_all_C2] = tsgEvaluate(lGrid_IEADstar, [C2(:,1) C2(:,2) C2(:,3) C2(:,4)]);
        % Matrix C_3
        [IEADstar_all_C3] = tsgEvaluate(lGrid_IEADstar, [C3(:,1) C3(:,2) C3(:,3) C3(:,4)]);
        % Matrix C_4
        [IEADstar_all_C4] = tsgEvaluate(lGrid_IEADstar, [C4(:,1) C4(:,2) C4(:,3) C4(:,4)]);
    end
 
    fprintf('IEAD* surrogate model evaluated \n')
    
    % ====================================================================
    %            Compute Moments per Case for Transformation
    % ====================================================================
    
    % Reload moments grid
    gridname = ['LS_' num2str(dim) 'D_Grid_Moments_level_' num2str(sglevel)];
    [lGrid_Moments] = tsgReloadGrid(gridname);
    
    % Evaluate moments surrogate model on points
    if dim == 2
        % Matrix A
        [Moments_A] = tsgEvaluate(lGrid_Moments, [A(:,1) A(:,2)]); 
        % Matrix B
        [Moments_B] = tsgEvaluate(lGrid_Moments, [B(:,1) B(:,2)]); 
        
        % Matrix C_1
        [Moments_C1] = tsgEvaluate(lGrid_Moments, [C1(:,1) C1(:,2)]); 
        % Matrix C_2
        [Moments_C2] = tsgEvaluate(lGrid_Moments, [C2(:,1) C2(:,2)]); 
    else
        % Matrix A
        [Moments_A] = tsgEvaluate(lGrid_Moments, [A(:,1) A(:,2) A(:,3) A(:,4)]); 
        % Matrix B
        [Moments_B] = tsgEvaluate(lGrid_Moments, [B(:,1) B(:,2) B(:,3) B(:,4)]); 

        % Matrix C_1
        [Moments_C1] = tsgEvaluate(lGrid_Moments, [C1(:,1) C1(:,2) C1(:,3) C1(:,4)]); 
        % Matrix C_2
        [Moments_C2] = tsgEvaluate(lGrid_Moments, [C2(:,1) C2(:,2) C2(:,3) C2(:,4)]); 
        % Matrix C_3
        [Moments_C3] = tsgEvaluate(lGrid_Moments, [C3(:,1) C3(:,2) C3(:,3) C3(:,4)]); 
        % Matrix C_4
        [Moments_C4] = tsgEvaluate(lGrid_Moments, [C4(:,1) C4(:,2) C4(:,3) C4(:,4)]); 
    end
    
    % Assign moments to arrays
    
    % Moments for Matrix A
    thetabar_mean_A = Moments_A(:,1);
    Ebar_mean_A = 10.^Moments_A(:,2);
    Theta_mean_11_A = 10.^Moments_A(:,3);
    Theta_mean_22_A = 10.^Moments_A(:,4);
    Theta_mean_12_A = Moments_A(:,5);
    Theta_mean_A = [Theta_mean_11_A Theta_mean_22_A Theta_mean_12_A]; % [Theta(1,1) Theta(2,2) Theta(1,2)];
    
    % Moments for Matrix B
    thetabar_mean_B = Moments_B(:,1);
    Ebar_mean_B = 10.^Moments_B(:,2);
    Theta_mean_11_B = 10.^Moments_B(:,3);
    Theta_mean_22_B = 10.^Moments_B(:,4);
    Theta_mean_12_B = Moments_B(:,5);
    Theta_mean_B = [Theta_mean_11_B Theta_mean_22_B Theta_mean_12_B]; % [Theta(1,1) Theta(2,2) Theta(1,2)];
    
    % Moments for Matrix C_1
    thetabar_mean_C1 = Moments_C1(:,1);
    Ebar_mean_C1 = 10.^Moments_C1(:,2);
    Theta_mean_11_C1 = 10.^Moments_C1(:,3);
    Theta_mean_22_C1 = 10.^Moments_C1(:,4);
    Theta_mean_12_C1 = Moments_C1(:,5);
    Theta_mean_C1 = [Theta_mean_11_C1 Theta_mean_22_C1 Theta_mean_12_C1]; % [Theta(1,1) Theta(2,2) Theta(1,2)];

    % Moments for Matrix C_2
    thetabar_mean_C2 = Moments_C2(:,1);
    Ebar_mean_C2 = 10.^Moments_C2(:,2);
    Theta_mean_11_C2 = 10.^Moments_C2(:,3);
    Theta_mean_22_C2 = 10.^Moments_C2(:,4);
    Theta_mean_12_C2 = Moments_C2(:,5);
    Theta_mean_C2 = [Theta_mean_11_C2 Theta_mean_22_C2 Theta_mean_12_C2]; % [Theta(1,1) Theta(2,2) Theta(1,2)];
    
    if dim == 4
        % Moments for Matrix C_3
        thetabar_mean_C3 = Moments_C3(:,1);
        Ebar_mean_C3 = 10.^Moments_C3(:,2);
        Theta_mean_11_C3 = 10.^Moments_C3(:,3);
        Theta_mean_22_C3 = 10.^Moments_C3(:,4);
        Theta_mean_12_C3 = Moments_C3(:,5);
        Theta_mean_C3 = [Theta_mean_11_C3 Theta_mean_22_C3 Theta_mean_12_C3]; % [Theta(1,1) Theta(2,2) Theta(1,2)];
        
        % Moments for Matrix C_4
        thetabar_mean_C4 = Moments_C4(:,1);
        Ebar_mean_C4 = 10.^Moments_C4(:,2);
        Theta_mean_11_C4 = 10.^Moments_C4(:,3);
        Theta_mean_22_C4 = 10.^Moments_C4(:,4);
        Theta_mean_12_C4 = Moments_C4(:,5);
        Theta_mean_C4 = [Theta_mean_11_C4 Theta_mean_22_C4 Theta_mean_12_C4]; % [Theta(1,1) Theta(2,2) Theta(1,2)];
    end
    
    fprintf('Moments surrogate model evaluated \n')

    % ====================================================================
    %              Read Binning in Transformed Coordinates
    % ====================================================================
    
    [~, ~, ~, ~, xstar, ystar, ~, ~, ~] = gridpolygonsTC;

    % ====================================================================
    %                    Compute Moments 
    % ====================================================================

    % Moments for Matrix A
    Moments_Theta_E_A = zeros(N,4);
    % Moments for Matrix B
    Moments_Theta_E_B = zeros(N,4);
    
    % Moments for Matrix C_1
    Moments_Theta_E_C1 = zeros(N,4);
    % Moments for Matrix C_2
    Moments_Theta_E_C2 = zeros(N,4);
    if dim == 4
        % Moments for Matrix C_3
        Moments_Theta_E_C3 = zeros(N,4);
        % Moments for Matrix C_4
        Moments_Theta_E_C4 = zeros(N,4);
    end
    
    % Ion temperature
    Ti = 10;
    
    % Create reference non-uniform grid in the transformed coordinates
    [~, ~, ~, ~, xNU_ref, yNU_ref, xNUpolyarray_ref, yNUpolyarray_ref, ~] = gridpolygonsTC;
    
    fprintf('Converting IEAD* to IEAD and computing moments ... \n\n')
    
    % Run over samples
    for i = 1:N
        
        fprintf('Sample %5g of %5g\n',i,N)
        
        % --------
        % Matrix A
        % --------
        
        % Read IEAD* for the ith sample 
        IEADstar_A = IEADstar_all_A(i,:);
        
        % Moments to convert data
        thetabar = thetabar_mean_A(i);
        Ebar = Ebar_mean_A(i);
        Theta = [Theta_mean_A(i,1) Theta_mean_A(i,3); Theta_mean_A(i,3) Theta_mean_A(i,2)];
        
        % Electron temperature
        Log10_Te_Ti = A(i,1);
        Te = Ti*10^Log10_Te_Ti;

        % Transform grid in the transformed coordinates to the original coordinates
        [xNU, yNU, xNUpolyarray, yNUpolyarray] = transform_gridpolygonsTC(xNU_ref, yNU_ref, xNUpolyarray_ref, yNUpolyarray_ref, thetabar, Ebar, Theta);
        
        % Evaluate IEAD from IEAD*
        [IEAD_A_NU] = evaluate_IEAD_from_IEADstar(xNU(:),yNU(:),thetabar,Ebar,Theta,xstar(:),ystar(:),IEADstar_A(:));
        
        % Mean and standard devation of angle & energy
        [~, ~, expected_theta, expected_E, std_theta, std_E] = MP_Expected_Theta_E(xNU,yNU,xNUpolyarray,yNUpolyarray,IEAD_A_NU,Te);
        Moments_Theta_E_A(i,:) = [expected_theta expected_E std_theta std_E];

        % --------
        % Matrix B
        % --------
        
        % Read IEAD* for the ith sample 
        IEADstar_B = IEADstar_all_B(i,:);
        
        % Moments to convert data
        thetabar = thetabar_mean_B(i);
        Ebar = Ebar_mean_B(i);
        Theta = [Theta_mean_B(i,1) Theta_mean_B(i,3); Theta_mean_B(i,3) Theta_mean_B(i,2)];
        
        % Electron temperature
        Log10_Te_Ti = B(i,1);
        Te = Ti*10^Log10_Te_Ti;

        % Transform grid in the transformed coordinates to the original coordinates
        [xNU, yNU, xNUpolyarray, yNUpolyarray] = transform_gridpolygonsTC(xNU_ref, yNU_ref, xNUpolyarray_ref, yNUpolyarray_ref, thetabar, Ebar, Theta);
        
        % Evaluate IEAD from IEAD*
        [IEAD_B_NU] = evaluate_IEAD_from_IEADstar(xNU(:),yNU(:),thetabar,Ebar,Theta,xstar(:),ystar(:),IEADstar_B(:));
        
        % Mean and standard devation of angle & energy
        [~, ~, expected_theta, expected_E, std_theta, std_E] = MP_Expected_Theta_E(xNU,yNU,xNUpolyarray,yNUpolyarray,IEAD_B_NU,Te);
        Moments_Theta_E_B(i,:) = [expected_theta expected_E std_theta std_E];

        % ----------
        % Matrix C_1
        % ----------
        
        % Read IEAD* for the ith sample 
        IEADstar_C1 = IEADstar_all_C1(i,:);
        
        % Moments to convert data
        thetabar = thetabar_mean_C1(i);
        Ebar = Ebar_mean_C1(i);
        Theta = [Theta_mean_C1(i,1) Theta_mean_C1(i,3); Theta_mean_C1(i,3) Theta_mean_C1(i,2)];     
 
        % Electron temperature
        Log10_Te_Ti = C1(i,1);
        Te = Ti*10^Log10_Te_Ti;

        % Transform grid in the transformed coordinates to the original coordinates
        [xNU, yNU, xNUpolyarray, yNUpolyarray] = transform_gridpolygonsTC(xNU_ref, yNU_ref, xNUpolyarray_ref, yNUpolyarray_ref, thetabar, Ebar, Theta);
        
        % Evaluate IEAD from IEAD*
        [IEAD_C1_NU] = evaluate_IEAD_from_IEADstar(xNU(:),yNU(:),thetabar,Ebar,Theta,xstar(:),ystar(:),IEADstar_C1(:));
        
        % Mean and standard devation of angle & energy
        [~, ~, expected_theta, expected_E, std_theta, std_E] = MP_Expected_Theta_E(xNU,yNU,xNUpolyarray,yNUpolyarray,IEAD_C1_NU,Te);
        Moments_Theta_E_C1(i,:) = [expected_theta expected_E std_theta std_E];

        % ----------
        % Matrix C_2
        % ----------
        
        % Read IEAD* for the ith sample 
        IEADstar_C2 = IEADstar_all_C2(i,:);
        
        % Moments to convert data
        thetabar = thetabar_mean_C2(i);
        Ebar = Ebar_mean_C2(i);
        Theta = [Theta_mean_C2(i,1) Theta_mean_C2(i,3); Theta_mean_C2(i,3) Theta_mean_C2(i,2)];
        
        % Electron temperature
        Log10_Te_Ti = C2(i,1);
        Te = Ti*10^Log10_Te_Ti;

        % Transform grid in the transformed coordinates to the original coordinates
        [xNU, yNU, xNUpolyarray, yNUpolyarray] = transform_gridpolygonsTC(xNU_ref, yNU_ref, xNUpolyarray_ref, yNUpolyarray_ref, thetabar, Ebar, Theta);
        
        % Evaluate IEAD from IEAD*
        [IEAD_C2_NU] = evaluate_IEAD_from_IEADstar(xNU(:),yNU(:),thetabar,Ebar,Theta,xstar(:),ystar(:),IEADstar_C2(:));
        
        % Mean and standard devation of angle & energy
        [~, ~, expected_theta, expected_E, std_theta, std_E] = MP_Expected_Theta_E(xNU,yNU,xNUpolyarray,yNUpolyarray,IEAD_C2_NU,Te);
        Moments_Theta_E_C2(i,:) = [expected_theta expected_E std_theta std_E];
        
        if dim == 4
            % ----------
            % Matrix C_3
            % ----------
            
            % Read IEAD* for the ith sample
            IEADstar_C3 = IEADstar_all_C3(i,:);
            
            % Moments to convert data
            thetabar = thetabar_mean_C3(i);
            Ebar = Ebar_mean_C3(i);
            Theta = [Theta_mean_C3(i,1) Theta_mean_C3(i,3); Theta_mean_C3(i,3) Theta_mean_C3(i,2)];
            
            % Electron temperature
            Log10_Te_Ti = C3(i,1);
            Te = Ti*10^Log10_Te_Ti;
            
            % Transform grid in the transformed coordinates to the original coordinates
            [xNU, yNU, xNUpolyarray, yNUpolyarray] = transform_gridpolygonsTC(xNU_ref, yNU_ref, xNUpolyarray_ref, yNUpolyarray_ref, thetabar, Ebar, Theta);
            
            % Evaluate IEAD from IEAD*
            [IEAD_C3_NU] = evaluate_IEAD_from_IEADstar(xNU(:),yNU(:),thetabar,Ebar,Theta,xstar(:),ystar(:),IEADstar_C3(:));
            
            % Mean and standard devation of angle & energy
            [~, ~, expected_theta, expected_E, std_theta, std_E] = MP_Expected_Theta_E(xNU,yNU,xNUpolyarray,yNUpolyarray,IEAD_C3_NU,Te);
            Moments_Theta_E_C3(i,:) = [expected_theta expected_E std_theta std_E];
            
            % ----------
            % Matrix C_4
            % ----------
            
            % Read IEAD* for the ith sample
            IEADstar_C4 = IEADstar_all_C4(i,:);
            
            % Moments to convert data
            thetabar = thetabar_mean_C4(i);
            Ebar = Ebar_mean_C4(i);
            Theta = [Theta_mean_C4(i,1) Theta_mean_C4(i,3); Theta_mean_C4(i,3) Theta_mean_C4(i,2)];
            
            % Electron temperature
            Log10_Te_Ti = C4(i,1);
            Te = Ti*10^Log10_Te_Ti;
            
            % Transform grid in the transformed coordinates to the original coordinates
            [xNU, yNU, xNUpolyarray, yNUpolyarray] = transform_gridpolygonsTC(xNU_ref, yNU_ref, xNUpolyarray_ref, yNUpolyarray_ref, thetabar, Ebar, Theta);
            
            % Evaluate IEAD from IEAD*
            [IEAD_C4_NU] = evaluate_IEAD_from_IEADstar(xNU(:),yNU(:),thetabar,Ebar,Theta,xstar(:),ystar(:),IEADstar_C4(:));
            
            % Mean and standard devation of angle & energy
            [~, ~, expected_theta, expected_E, std_theta, std_E] = MP_Expected_Theta_E(xNU,yNU,xNUpolyarray,yNUpolyarray,IEAD_C4_NU,Te);
            Moments_Theta_E_C4(i,:) = [expected_theta expected_E std_theta std_E];  
        end
        
    end

    % ====================================================================
    %   Compute Sensitivity Indices for Angle and Energy Moments
    % ====================================================================

    % First-order indices
    
    f0sq_Moments = ((1/N)*sum(Moments_Theta_E_A,1)).^2;
    
    Si_denom_Moments = (1/N)*sum(Moments_Theta_E_A.^2,1) - f0sq_Moments;
    
    S1_Moments = ((1/N)*sum(Moments_Theta_E_A.*Moments_Theta_E_C1,1) - f0sq_Moments)./Si_denom_Moments;
    S2_Moments = ((1/N)*sum(Moments_Theta_E_A.*Moments_Theta_E_C2,1) - f0sq_Moments)./Si_denom_Moments;
    if dim == 4
        S3_Moments = ((1/N)*sum(Moments_Theta_E_A.*Moments_Theta_E_C3,1) - f0sq_Moments)./Si_denom_Moments;
        S4_Moments = ((1/N)*sum(Moments_Theta_E_A.*Moments_Theta_E_C4,1) - f0sq_Moments)./Si_denom_Moments;
    end
            
    % Total-effect indices
    if dim == 4    
        ST1_Moments = 1 - ((1/N)*sum(Moments_Theta_E_B.*Moments_Theta_E_C1,1) - f0sq_Moments)./Si_denom_Moments;
        ST2_Moments = 1 - ((1/N)*sum(Moments_Theta_E_B.*Moments_Theta_E_C2,1) - f0sq_Moments)./Si_denom_Moments;
        ST3_Moments = 1 - ((1/N)*sum(Moments_Theta_E_B.*Moments_Theta_E_C3,1) - f0sq_Moments)./Si_denom_Moments;
        ST4_Moments = 1 - ((1/N)*sum(Moments_Theta_E_B.*Moments_Theta_E_C4,1) - f0sq_Moments)./Si_denom_Moments;
    end
    
    % ====================================================================
    %      Plot Sensitivity Indices for Moments of Angle and Energy
    % ====================================================================    

    figure
    set(gcf,'Position',[500 500 1400 500])

    % -----------
    %     S_1
    % -----------
    if dim == 2
        subplot(1,3,1)
    else
        subplot(2,4,1)  
    end
    
    hold on
    
    % Bar width
    wd = 0.5;
    
    x = [1 2 3 4];
    bar(x,[S1_Moments(1) nan nan nan],wd,'facecolor','r');
    bar(x,[nan S1_Moments(2) nan nan],wd,'facecolor','b');
    bar(x,[nan nan S1_Moments(3) nan],wd,'facecolor','g');
    bar(x,[nan nan nan S1_Moments(4)],wd,'facecolor','m');
    
    set(gca,'xtick',x,'xticklabel',{'$\mu_{\theta}$';'$\mu_{E}$';'$\sigma_{\theta}$';'$\sigma_{E}$'})

    title('$S_1$','interpreter','latex','FontSize',30)
    set(gca,'FontSize',20)
    set(gca,'TickLabelInterpreter','latex')
    
    box('on')
    ylim([0 1])
    pbaspect([1 1 1])
    
    % -----------
    %     S_2
    % -----------
    if dim == 2
        subplot(1,3,2)
    else
        subplot(2,4,2)  
    end
    
    hold on
    
    % Bar width
    wd = 0.5;
    
    x = [1 2 3 4];
    bar(x,[S2_Moments(1) nan nan nan],wd,'facecolor','r');
    bar(x,[nan S2_Moments(2) nan nan],wd,'facecolor','b');
    bar(x,[nan nan S2_Moments(3) nan],wd,'facecolor','g');
    bar(x,[nan nan nan S2_Moments(4)],wd,'facecolor','m');
    
    set(gca,'xtick',x,'xticklabel',{'$\mu_{\theta}$';'$\mu_{E}$';'$\sigma_{\theta}$';'$\sigma_{E}$'})

    title('$S_2$','interpreter','latex','FontSize',30)
    set(gca,'FontSize',20)
    set(gca,'TickLabelInterpreter','latex')
    
    box('on')
    ylim([0 1])
    pbaspect([1 1 1])
    
    % -----------
    %     S_12   (only for dim = 2)
    % -----------
    if dim == 2
        
        S12_Moments = 1 - S1_Moments - S2_Moments;
        
        subplot(1,3,3)
        hold on
        
        % Bar width
        wd = 0.5;
        
        x = [1 2 3 4];
        bar(x,[S12_Moments(1) nan nan nan],wd,'facecolor','r');
        bar(x,[nan S12_Moments(2) nan nan],wd,'facecolor','b');
        bar(x,[nan nan S12_Moments(3) nan],wd,'facecolor','g');
        bar(x,[nan nan nan S12_Moments(4)],wd,'facecolor','m');
        
        set(gca,'xtick',x,'xticklabel',{'$\mu_{\theta}$';'$\mu_{E}$';'$\sigma_{\theta}$';'$\sigma_{E}$'})
        
        
        title('$S_{12}$','interpreter','latex','FontSize',30)
        set(gca,'FontSize',20)
        set(gca,'TickLabelInterpreter','latex')
        
        box('on')
        ylim([0 1])
        pbaspect([1 1 1])
        
    else
        
        % -----------
        %     S_3
        % -----------
        subplot(2,4,3)
        hold on

        % Bar width
        wd = 0.5;

        x = [1 2 3 4];
        bar(x,[S3_Moments(1) nan nan nan],wd,'facecolor','r');
        bar(x,[nan S3_Moments(2) nan nan],wd,'facecolor','b');
        bar(x,[nan nan S3_Moments(3) nan],wd,'facecolor','g');
        bar(x,[nan nan nan S3_Moments(4)],wd,'facecolor','m');

        set(gca,'xtick',x,'xticklabel',{'$\mu_{\theta}$';'$\mu_{E}$';'$\sigma_{\theta}$';'$\sigma_{E}$'})

        title('$S_3$','interpreter','latex','FontSize',30)
        set(gca,'FontSize',20)
        set(gca,'TickLabelInterpreter','latex')

        box('on')
        ylim([0 1])
        pbaspect([1 1 1])

        % -----------
        %     S_4
        % -----------
        subplot(2,4,4)
        hold on

        % Bar width
        wd = 0.5;

        x = [1 2 3 4];
        bar(x,[S4_Moments(1) nan nan nan],wd,'facecolor','r');
        bar(x,[nan S4_Moments(2) nan nan],wd,'facecolor','b');
        bar(x,[nan nan S4_Moments(3) nan],wd,'facecolor','g');
        bar(x,[nan nan nan S4_Moments(4)],wd,'facecolor','m');

        set(gca,'xtick',x,'xticklabel',{'$\mu_{\theta}$';'$\mu_{E}$';'$\sigma_{\theta}$';'$\sigma_{E}$'})

        title('$S_4$','interpreter','latex','FontSize',30)
        set(gca,'FontSize',20)
        set(gca,'TickLabelInterpreter','latex')

        box('on')
        ylim([0 1])
        pbaspect([1 1 1])
        
        % -----------
        %  S_{T_1}
        % -----------
        subplot(2,4,5)
        hold on

        % Bar width
        wd = 0.5;

        x = [1 2 3 4];
        bar(x,[ST1_Moments(1) nan nan nan],wd,'facecolor','r');
        bar(x,[nan ST1_Moments(2) nan nan],wd,'facecolor','b');
        bar(x,[nan nan ST1_Moments(3) nan],wd,'facecolor','g');
        bar(x,[nan nan nan ST1_Moments(4)],wd,'facecolor','m');

        set(gca,'xtick',x,'xticklabel',{'$\mu_{\theta}$';'$\mu_{E}$';'$\sigma_{\theta}$';'$\sigma_{E}$'})

        title('$S_{T_1}$','interpreter','latex','FontSize',30)
        set(gca,'FontSize',20)
        set(gca,'TickLabelInterpreter','latex')

        box('on')
        ylim([0 1])
        pbaspect([1 1 1])

        % -----------
        %  S_{T_2}
        % -----------
        subplot(2,4,6)
        hold on

        % Bar width
        wd = 0.5;

        x = [1 2 3 4];
        bar(x,[ST2_Moments(1) nan nan nan],wd,'facecolor','r');
        bar(x,[nan ST2_Moments(2) nan nan],wd,'facecolor','b');
        bar(x,[nan nan ST2_Moments(3) nan],wd,'facecolor','g');
        bar(x,[nan nan nan ST2_Moments(4)],wd,'facecolor','m');

        set(gca,'xtick',x,'xticklabel',{'$\mu_{\theta}$';'$\mu_{E}$';'$\sigma_{\theta}$';'$\sigma_{E}$'})

        title('$S_{T_2}$','interpreter','latex','FontSize',30)
        set(gca,'FontSize',20)
        set(gca,'TickLabelInterpreter','latex')

        box('on')
        ylim([0 1])
        pbaspect([1 1 1])

        % -----------
        %  S_{T_3}
        % -----------
        subplot(2,4,7)
        hold on

        % Bar width
        wd = 0.5;

        x = [1 2 3 4];
        bar(x,[ST3_Moments(1) nan nan nan],wd,'facecolor','r');
        bar(x,[nan ST3_Moments(2) nan nan],wd,'facecolor','b');
        bar(x,[nan nan ST3_Moments(3) nan],wd,'facecolor','g');
        bar(x,[nan nan nan ST3_Moments(4)],wd,'facecolor','m');

        set(gca,'xtick',x,'xticklabel',{'$\mu_{\theta}$';'$\mu_{E}$';'$\sigma_{\theta}$';'$\sigma_{E}$'})

        title('$S_{T_3}$','interpreter','latex','FontSize',30)
        set(gca,'FontSize',20)
        set(gca,'TickLabelInterpreter','latex')

        box('on')
        ylim([0 1])
        pbaspect([1 1 1])

        % -----------
        %  S_{T_4}
        % -----------
        subplot(2,4,8)
        hold on

        % Bar width
        wd = 0.5;

        x = [1 2 3 4];
        bar(x,[ST4_Moments(1) nan nan nan],wd,'facecolor','r');
        bar(x,[nan ST4_Moments(2) nan nan],wd,'facecolor','b');
        bar(x,[nan nan ST4_Moments(3) nan],wd,'facecolor','g');
        bar(x,[nan nan nan ST4_Moments(4)],wd,'facecolor','m');

        set(gca,'xtick',x,'xticklabel',{'$\mu_{\theta}$';'$\mu_{E}$';'$\sigma_{\theta}$';'$\sigma_{E}$'})

        title('$S_{T_4}$','interpreter','latex','FontSize',30)
        set(gca,'FontSize',20)
        set(gca,'TickLabelInterpreter','latex')

        box('on')
        ylim([0 1])
        pbaspect([1 1 1])
        
    end

end
