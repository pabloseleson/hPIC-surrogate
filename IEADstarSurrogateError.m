
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The function IEADstarSurrogateError computes and plots the error 
% of the IEAD* surrogate model against (training or testing) data 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Input
% -----
% sglevel: data sparse grid level
% dim: dimension of parameter space
%      2D: dim = 2
%      4D: dim = 4
% data_type: 'train' for training data and 'test' for testing data

% Author: Pablo Seleson
% ------

% Last Modified: March 9, 2022
% -------------

function IEADstarSurrogateError(sglevel,dim,data_type)

     close all
    
     % ====================================================================
     %                      Load IEAD data
     % ====================================================================
     
     IEAD_data_filename = ['IEAD_data_' num2str(dim) 'D/Level_' num2str(sglevel) '_' data_type '.mat'];
     
     load(IEAD_data_filename)
    
     % Number of bins
     dimXstar = length(XNODESstar);
     dimYstar = length(YNODESstar);

    % ====================================================================
    %       Evaluate IEAD* surrogate model at the grid points
    % ====================================================================

    % Reload IEAD* fit grid
    gridname = ['LS_' num2str(dim) 'D_Grid_IEADstar_level_' num2str(sglevel)];
    [lGrid_IEADstar] = tsgReloadGrid(gridname);

    if dim == 2
        % Evaluate 2D IEAD* surrogate model at the 2D grid points
        [result_all] = tsgEvaluate(lGrid_IEADstar, [points(:,1) points(:,2)]);
    else
        % Evaluate 4D IEAD* surrogate model at the 4D grid points
        [result_all] = tsgEvaluate(lGrid_IEADstar, [points(:,1) points(:,2) points(:,3) points(:,4)]);
    end

    % Number of cases
    ncases = length(points);

    % Initialize error norm array
    l2_err_norm = zeros(ncases,1);

    % Compute bin areas
    xstarmin = -5;
    xstarmax =  5;
    ystarmin = -5;
    ystarmax =  5;
    [~,~,~,pareastar] = compute_2D_NU_distribution(0,0,XNODESstar,YNODESstar,xstarmin,xstarmax,ystarmin,ystarmax);

    % Run over cases
    for ncase = 1:ncases

        % ----------------------------------
        %        Read data and fit
        % ----------------------------------

        % Data
        IEADstar = reshape(IEADstar_array(ncase,:),dimXstar,dimYstar);

        % Fit
        result = result_all(ncase,:);
        IEADstar_fit = reshape(result,dimXstar,dimYstar);

        % ----------------------------------
        %         Compute error
        % ----------------------------------

        Err = IEADstar_fit - IEADstar;

        % Note: The norms apply weights given by the bin areas
        l2_err_norm(ncase) = sqrt(sum(sum((Err.^2).*pareastar)))/sqrt(sum(sum((IEADstar.^2).*pareastar)));

    end

    % ====================================================================
    %                         Plot errors
    % ====================================================================

    % --------------------------------
    %       Error histogram
    % --------------------------------

    % Create figure
    fig1 = figure;
    axes1 = axes('Parent',fig1);
    hold(axes1,'all');
    set(axes1,'FontSize',20);

    % Plot histogram
    histogram(l2_err_norm);

    % xlabel
    xlabel('Rel. $L^2$ norm of IEAD* Error','Interpreter','latex','Fontsize',24)
    box('on')

    % Position figure
    set(gcf, 'Position',  [600, 500, 600, 500])

    % Draw vertical lines at mean and mean +/- 1 standard deviation
    yl = ylim;
    mu = mean(l2_err_norm);
    sigma = std(l2_err_norm);
    plot([mu mu],yl,'--r','LineWidth',3);
    plot([mu+sigma mu+sigma],yl,':r','LineWidth',2);
    plot([mu-sigma mu-sigma],yl,':r','LineWidth',2);

    % Set axis limits
    if strcmp(data_type,'train')
        if dim == 2
            xlim([-.002 0.146])
            xticks([0 0.05 0.1])
        else
            xlim([-.002 0.48])
        end
    elseif strcmp(data_type,'test')
        if dim == 2
            xlim([0.0045 0.251])
        else
            xlim([-.002 0.55])
        end
    else
        error('data_type should be train or test.')
    end

    % --------------------------------
    %  Errors scatter plot on 2D grid
    % --------------------------------

    if dim == 2

        % Create figure
        fig2 = figure;
        axes2 = axes('Parent',fig2);
        hold(axes2,'all');
        set(axes2,'FontSize',20);

        % Set colormap and color limits
        cm = colormap(cool);

        if strcmp(data_type,'train')
            maxz = 0.15;
        elseif strcmp(data_type,'test')
            maxz = 0.22;
        else
            error('data_type should be train or test.')
        end

        caxis([0 maxz]);

        % Create colors based on data
        zref = linspace(0,maxz,length(cm));
        z = l2_err_norm;
        [~, Ivals] = min(abs(z' - zref'));
        pColor = cm(Ivals,:);

        % Points sizes
        pSize = 200*Ivals/length(cm);

        % Scatter plot
        scatter(points(:,1),points(:,2),pSize,pColor,'s','filled');

        % Other settings
        axis([min(points(:,1)) max(points(:,1)) min(points(:,2)) max(points(:,2))])
        xlabel('$\log_{10}(T_e/T_i)$','Interpreter','latex','Fontsize',24)
        ylabel('$\Psi$','Interpreter','latex','Fontsize',24)
        colorbar
        box('on')

        % Colorbar title
        hcb=colorbar;
        hcb.Title.String = "$e$";
        set(hcb.Title,'FontSize',30,'Interpreter','latex')

        % Position figure
        set(gcf, 'Position',  [1300, 500, 600, 500])
        
    end

end


     
