%% STA & SCA Plotting Code - 1
% 26,01,2021 Onwards


%% Plotting Choices

% Choose which cell RF to plot
Cell_to_Plot  = 19; % 7
RF_Ident_Plot = RF_Ident{Cell_to_Plot,1};

% Choose whether to use B&W or coloured heatmaps
% 1 = B&W;
% 2 = coloured.
Heat_Map_Colour_Choice = 2;

% Allow for 2 colour cases: monochromatic and tetrachromatic (RGBUV)
Spectral_Names = {'R','G','B','UV'};

colorMap_arr        = NaN(256,3,4);
colorMap_arr(:,:,1) = [linspace(0,1,256)', zeros(256,2)];                     % Red
colorMap_arr(:,:,2) = [zeros(256,1),linspace(0,1,256)', zeros(256,1)];        % Green
colorMap_arr(:,:,3) = [zeros(256,2), linspace(0,1,256)'];                     % Blue
colorMap_arr(:,:,4) = [linspace(0,1,256)', zeros(256,1), linspace(0,1,256)']; % UV
%colorMapBlack      = [0 0 0];

First_Bin = 1;
Last_Bin  = 10;
Num_Selected_Bins = Last_Bin-First_Bin+1;

%% Plot

% STA Full
if p.Plot_Choice == 1
    if p.RF_Ident_Meth_vec(1) == 1 % STA-SD method;
        figSTA = figure;
        if p.Spectral_Dim == 1
        elseif p.Spectral_Dim == 4
            if Heat_Map_Colour_Choice == 1
                colormap(figSTA,gray(256));
            end
            for i = 1:p.Spectral_Dim
                for j = First_Bin:Last_Bin%1:p.Num_STE_bins
                    loop_plot = subplot(p.Spectral_Dim,Num_Selected_Bins,Num_Selected_Bins*(i-1)+j);%loop_plot = subplot(p.Spectral_Dim,p.Num_STE_bins,p.Num_STE_bins*(i-1)+j);
                    imagesc(RF_Ident_Plot.STA(:,:,j,i)); hold on;
                    if Heat_Map_Colour_Choice == 2
                        colormap(loop_plot,colorMap_arr(:,:,i));
                        if p.RF_Type(2) == 1 % All Significant Pixels
                            if ~isempty(RF_Ident_Plot.STASD_ASP_Num_RF_pixels{i})
                                scatter(RF_Ident_Plot.STASD_ASP_RF_coords{i}(:,2),RF_Ident_Plot.STASD_ASP_RF_coords{i}(:,1),'wo','LineWidth',1.5);
                            end
                        end
                    else % Heat_Map_Colour_Choice == 1
                        if p.RF_Type(2) == 1 % All Significant Pixels
                            if ~isempty(RF_Ident_Plot.STASD_ASP_Num_RF_pixels{i})
                                scatter(RF_Ident_Plot.STASD_ASP_RF_coords{i}(:,2),RF_Ident_Plot.STASD_ASP_RF_coords{i}(:,1),'go','LineWidth',1.5);
                            end
                        end
                    end
                    axis equal; axis tight;
                    if i==p.Spectral_Dim
                        xlabel('x');
                    end
                    if j==1
                        ylabel('y');
                    end
                    title(Spectral_Names{i});
                    if j == Num_Selected_Bins%p.Num_STE_bins
                        colorbar;
                    end
                    set(gca,'FontSize',12);
                end
            end
        end
        annotation('textbox',[.5 .975 0 0],'String','STA with STA-SD ASP RF','FitBoxToText','on','FontSize',12,'FontWeight','bold','LineStyle','none','HorizontalAlignment','center','VerticalAlignment','middle');
        set(gcf,'color','w');
    end
end


% SCA Full
if p.Plot_Choice == 1
    if p.RF_Ident_Meth_vec(4) == 1 % SCA method;
        figSCA = figure;
        if p.Spectral_Dim == 1
        elseif p.Spectral_Dim == 4
            if Heat_Map_Colour_Choice == 1
                colormap(figSCA,gray(256));
            end
            for i = 1:p.Spectral_Dim
                for j = First_Bin:Last_Bin%1:p.Num_STE_bins
                    loop_plot = subplot(p.Spectral_Dim,Num_Selected_Bins,Num_Selected_Bins*(i-1)+j); %subplot(p.Spectral_Dim,p.Num_STE_bins,p.Num_STE_bins*(i-1)+j);
                    imagesc(RF_Ident_Plot.SCA_Stixel_covar(:,:,j,i)); hold on;
                    if Heat_Map_Colour_Choice == 2
                        colormap(loop_plot,colorMap_arr(:,:,i));
                        if p.RF_Type(2) == 1 % All Significant Pixels
                            if ~isempty(RF_Ident_Plot.SC_ASP_Num_RF_pixels{i})
                                scatter(RF_Ident_Plot.SC_ASP_RF_coords{i}(:,2),RF_Ident_Plot.SC_ASP_RF_coords{i}(:,1),'wo','LineWidth',1.5);
                            end
                        end
                    else % Heat_Map_Colour_Choice == 1
                        if p.RF_Type(2) == 1 % All Significant Pixels
                            if ~isempty(RF_Ident_Plot.SC_ASP_Num_RF_pixels{i})
                                scatter(RF_Ident_Plot.SC_ASP_RF_coords{i}(:,2),RF_Ident_Plot.SC_ASP_RF_coords{i}(:,1),'go','LineWidth',1.5);
                            end
                        end
                    end
                    axis equal; axis tight;
                    if i==p.Spectral_Dim
                        xlabel('x');
                    end
                    if j==1
                        ylabel('y');
                    end
                    title(Spectral_Names{i});
                    if j == Num_Selected_Bins%p.Num_STE_bins
                        colorbar;
                    end
                    set(gca,'FontSize',12);
                end
            end
        end
        annotation('textbox',[.5 .975 0 0],'String','Self Covar. with SC ASP RF','FitBoxToText','on','FontSize',12,'FontWeight','bold','LineStyle','none','HorizontalAlignment','center','VerticalAlignment','middle');
        set(gcf,'color','w');
    end
end


