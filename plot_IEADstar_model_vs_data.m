
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The function plot_IEADstar_model_vs_data plots a comparison between the 
% IEAD* from the model and the data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Input
% ----
% sglevel: data sparse grid level
% dim: dimension of parameter space
%      2D: dim = 2
%      4D: dim = 4
% data_type: 'train' for training data and 'test' for testing data

% Author: Pablo Seleson
% ------

% Last Modified: February 1, 2022
% -------------

function plot_IEADstar_model_vs_data(sglevel,dim,data_type)

    % ===================================================================
    %                        Load IEAD Data
    % ===================================================================
    
    IEAD_data_filename = ['IEAD_data_' num2str(dim) 'D/Level_' num2str(sglevel) '_' data_type '.mat'];
    
    load(IEAD_data_filename)
    
    % Number of bins
    dimXstar = length(XNODESstar);
    dimYstar = length(YNODESstar);
    
    % Number of cases
    ncases = length(points);

    % ===================================================================
    %        Evaluate IEAD* surrogate model at the grid points
    % ===================================================================

    % Reload IEAD* fit grid
    gridname = ['LS_' num2str(dim) 'D_Grid_IEADstar_level_' num2str(sglevel)];
    [lGrid_IEADstar] = tsgReloadGrid(gridname);

    if dim == 2
        % Evaluate 2D IEAD* surrogate model at the 2D grid points
        [IEADstar_all] = tsgEvaluate(lGrid_IEADstar, [points(:,1) points(:,2)]);
    else
        % Evaluate 4D IEAD* surrogate model at the 4D grid points
        [IEADstar_all] = tsgEvaluate(lGrid_IEADstar, [points(:,1) points(:,2) points(:,3) points(:,4)]);
    end

    % ------------------------------------------------------
    %                   Create grid
    % ------------------------------------------------------

    % Create grid in the transformed coordinates
    [XNODESstar, YNODESstar, dimXstar, dimYstar, xstar, ystar, xstarpolyarray, ystarpolyarray, elementstar] = gridpolygonsTC;

    % ------------------------------------------------------
    %                   Open figure
    % ------------------------------------------------------

    figure
    set(gcf,'Position',[300 400 1800 700])

    % Run over cases
    for ncase = 1:ncases
        
        % ----------------------------------
        % Find values of physical parameters
        % ----------------------------------
        if dim == 2
            Log10_Te_Ti_value = points(ncase,1);
            Psi_value = points(ncase,2);
        else
            Log10_Te_Ti_value = points(ncase,1);
            Psi_value = points(ncase,2);
            B_value = points(ncase,3);
            Log10_n_value = points(ncase,4);
        end

        % ----------------------------------
        %             Title
        % ----------------------------------

        string_fig = ['IEAD* -- Case = ' num2str(ncase) ' of ' num2str(ncases)];
        if dim == 2
            string_case = ['$\log_{10}(T_e/T_i)$: ' num2str(Log10_Te_Ti_value,'%.2f')  ' $\Psi$: ' num2str(Psi_value,'%.1f')];
        else
            string_case = ['$\log_{10}(T_e/T_i)$: ' num2str(Log10_Te_Ti_value,'%.2f')  ' $\Psi$: ' num2str(Psi_value,'%.1f') ' $B$: ' num2str(B_value,'%.1f') ' $\log_{10}(n)$: ' num2str(Log10_n_value,'%.1f')];
        end
        sgtitle({string_fig,string_case},'interpreter','latex','FontSize',24)

        % ----------------------------------
        %           Plot fit
        % ----------------------------------

        sp1 = subplot(1,2,1);

        % Fit
        IEADstar_fit = IEADstar_all(ncase,:);

        % Find IEAD* limits
        max_value = max(IEADstar_fit);
        min_value = min(IEADstar_fit);
        caxis([min_value max_value]);

        if ncase > 1
            delete(p1);
        end

        % Plot IEAD* fit
        p1 = patch(xstarpolyarray, ystarpolyarray,IEADstar_fit,'LineStyle','none');
        axis([min(xstarpolyarray(:)) max(xstarpolyarray(:)) min(ystarpolyarray(:)) max(ystarpolyarray(:))])
        
        % Colorbar
        colorbar

        % Other settings
        set(sp1,'FontSize',20)
        title('Model','interpreter','latex','FontSize',24)
        xlabel('$\theta^*$','Interpreter','latex','Fontsize',24)
        ylabel('$E^*$','Interpreter','latex','Fontsize',24)
        pbaspect([1 1 1])
        box('on')

        % Colorbar title
        hcb=colorbar;
        hcb.Title.String = "$f^*$";
        set(hcb.Title,'FontSize',30,'Interpreter','latex')

        % ----------------------------------
        %           Plot data
        % ----------------------------------

        sp2 = subplot(1,2,2);

        % Data
        IEADstar = reshape(IEADstar_array(ncase,:),dimXstar,dimYstar);

        % IEAD* limits
        caxis([min_value max_value]);

        if ncase > 1
            delete(p2);
        end

        % Plot IEAD* fit
        p2 = patch(xstarpolyarray, ystarpolyarray,IEADstar(:),'LineStyle','none');
        axis([min(xstarpolyarray(:)) max(xstarpolyarray(:)) min(ystarpolyarray(:)) max(ystarpolyarray(:))])
        
        % Colorbar
        colorbar

        % Other settings
        set(sp2,'FontSize',20)
        title('Data','interpreter','latex','FontSize',24)
        xlabel('$\theta^*$','Interpreter','latex','Fontsize',24)
        ylabel('$E^*$','Interpreter','latex','Fontsize',24)
        pbaspect([1 1 1])
        box('on')

        % Colorbar title
        hcb=colorbar;
        hcb.Title.String = "$f^*$";
        set(hcb.Title,'FontSize',30,'Interpreter','latex')

        drawnow

    end

end