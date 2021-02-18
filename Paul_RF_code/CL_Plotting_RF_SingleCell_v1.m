%% Plot RF and Associated Summary Statistics - Single Cell
function out = CL_Plotting_RF_SingleCell_v1 (savepath, add_info)
% Choose which cell RF to plot
% THis can either be the cell that is highlighted in the table, or a user
% selection by ticking the column "Plot selection"
cells_selected = add_info.tables.RF_single_cell.Data.("Plot selection");
if any(cells_selected)

    cell_indices = add_info.tables.RF_single_cell.Data.("Cell idx");

    cell_to_plot = cell_indices(cells_selected);
else
    cell_to_plot  = add_info.settings.kernel_new.cell_to_plot;
end
Spectral_Dim = 4;
Spectral_Names = {'Red','Green','Blue','VS or UV'};
%RF_Ident_Plot = RF_Ident{Cell_to_Plot,1};

%Delete existing figures
delete(add_info.panels.BoxRF.Children)
delete(add_info.panels.AllPxRF.Children)
delete(add_info.panels.GaussRF.Children)
drawnow

% add_info.panels.BoxRF.Title = [];
% add_info.panels.AllPxRF.Title = [];
% add_info.panels.GaussRF.Title = [];

%Colormap
colorMap_arr        = NaN(256,3,4);
colorMap_arr(:,:,1) = [linspace(0,1,256)', zeros(256,2)];                     % Red
colorMap_arr(:,:,2) = [zeros(256,1),linspace(0,1,256)', zeros(256,1)];        % Green
colorMap_arr(:,:,3) = [zeros(256,2), linspace(0,1,256)'];                     % Blue
colorMap_arr(:,:,4) = [linspace(0,1,256)', zeros(256,1), linspace(0,1,256)']; % UV



for ii = 1:size(add_info.settings.kernel_new.folderplot,2)
for kk = 1:length(add_info.settings.kernel_new.cell_to_plot)
%try
plot_folder = add_info.settings.kernel_new.folderplot{:,ii};
plot_folder_legend = strrep(plot_folder,'_',' ');
plot_folder_legend = ['Cell ',num2str(cell_to_plot(kk,1)),' ',plot_folder_legend];

%% Decide which data shall be plotted

if contains(plot_folder,'SS')
        
    RF_type = 1;
elseif contains(plot_folder,'SC')
     
    RF_type = 2;
end






%% Plot selected folder and selected RF type
%We have 3 different options: Box, all Pixel, Gaussian
%Check which of those the user wants to be plotted and plot. The legend of
%the plots is according to the name of the folder selected by the user.
%First, we load the overview variable, as much smaller and we can get 
%the index in the structure than we index directly into the .matfile.

% Load data
L = load(findfile_app(add_info.stim_idx,savepath,"RF_overview.mat",'subfolder',plot_folder));
index = find([L.RF_overview.cell_idx] == cell_to_plot(kk,1));
RF_file = strcat(L.RF_overview(index).file,'.bin');

fileID = fopen(RF_file,'r');
RF_Ident = fread(fileID);
RF_Ident = hlp_deserialize(RF_Ident);

RF_Ident_Plot = RF_Ident.RF_results;


if RF_type == 1
if add_info.settings.kernel_new.plot_Box  % Box
        % Organize figures
    if size(add_info.settings.kernel_new.folderplot,2)*...
            length(add_info.settings.kernel_new.cell_to_plot) == 1 ...
            && ~add_info.settings.kernel_new.plot_in_new_window
       fig1 = add_info.panels.BoxRF;
    else
       fig1 = figure;
    end
    %fig1 = figure;
    if Spectral_Dim == 1 %has to be not hard coded
        colormap(fig1,gray(256));
%                 imagesc(RF_Ident_Plot.STA_SD); hold on;
%                 scatter(RF_Ident_Plot.RF_coords_centre(2),RF_Ident_Plot.RF_coords_centre(1),'rx','LineWidth',1.5);
%                 scatter(RF_Ident_Plot.RF_coords_noncentre(:,2),RF_Ident_Plot.RF_coords_noncentre(:,1),'go','LineWidth',1.5);
        axis equal; axis tight;
        xlabel('x');
        ylabel('y');
        title('RF');
        if add_info.settings.kernel_new.plot_colorbar
            colorbar(loop_plot,'southoutside');
        end
        set(gca,'FontSize',12);
    elseif Spectral_Dim == 4
%         if add_info.settings.kernel_new.Heat_gray
%             colormap(add_info.panels.BoxRF,gray(256));
%         end
        for i = 1:Spectral_Dim
            
            loop_plot = subplot(2,2,i,'parent',fig1);
            imagesc(loop_plot,RF_Ident_Plot.STA_SD(:,:,i));
            hold(loop_plot,'on');
            if add_info.settings.kernel_new.Heat_colour
                colormap(loop_plot,colorMap_arr(:,:,i));
            end
             if add_info.settings.kernel_new.Heat_gray
                colormap(loop_plot,gray(256));
             end
            if ~isempty(RF_Ident_Plot.STASD_Box_Num_RF_pixels{i})
                scatter(loop_plot,RF_Ident_Plot.STASD_Box_RF_coords_centre{i}(2),...
                    RF_Ident_Plot.STASD_Box_RF_coords_centre{i}(1),'rx','LineWidth',1.5,...
                    'MarkerFaceAlpha',0.1);
                scatter(loop_plot,RF_Ident_Plot.STASD_Box_RF_coords_noncentre{i}(:,2)...
                    ,RF_Ident_Plot.STASD_Box_RF_coords_noncentre{i}(:,1),'go','LineWidth',1.5,...
                    'MarkerFaceAlpha',0.1);
            end
            axis(loop_plot,'equal'); axis(loop_plot, 'tight');
            if i>2
                xlabel(loop_plot,'x');
            end
            if i==1||i==3
                ylabel(loop_plot,'y');
            end
            %title(Spectral_Names{i});
            %if i == 4
                original_ax = get(loop_plot,'Position');
            if add_info.settings.kernel_new.plot_colorbar
                colorbar(loop_plot,'southoutside');
            end
                set(loop_plot,'Position',original_ax);
            %end
            title(loop_plot,Spectral_Names{i});
            set(loop_plot,'FontSize',12);
        end
        hold(loop_plot,'off');
       
    end
    
    
    if strcmp(get(fig1,'type'),'uipanel')
         sgtitle(fig1,plot_folder_legend);
    else
        fig1.Name = plot_folder_legend;
        fig1.NumberTitle = 'off';
        sgtitle(fig1,'Box');

    end
%     annotation(fig1,'textbox',[.5 .975 0 0],'String',plot_folder_legend,...
%         'FitBoxToText','on','FontSize',12,'FontWeight','bold','LineStyle',...
%         'none','HorizontalAlignment','center','VerticalAlignment','middle');
    %set(gcf,'color','w');

end

if add_info.settings.kernel_new.plot_Allpixel  % All Significant Pixels

    if size(add_info.settings.kernel_new.folderplot,2)*...
            length(add_info.settings.kernel_new.cell_to_plot) == 1 ...
            && ~add_info.settings.kernel_new.plot_in_new_window
       fig2 = add_info.panels.AllPxRF;
    else
       fig2 = figure;
    end
    
    
    if Spectral_Dim == 1
        colormap(fig2,gray(256));
%                 imagesc(RF_Ident_Plot.STA_SD); hold on;
%                 scatter(RF_Ident_Plot.RF_coords(:,2),RF_Ident_Plot.RF_coords(:,1),'go','LineWidth',1.5);
        axis equal; axis tight;
        xlabel('x');
        ylabel('y');
        title('RF');
        if add_info.settings.kernel_new.plot_colorbar
           colorbar(loop_plot,'southoutside');
        end
        set(gca,'FontSize',12);
    elseif Spectral_Dim == 4
        
        for i = 1:Spectral_Dim
            loop_plot = subplot(2,2,i,'parent',fig2);
            imagesc(loop_plot,RF_Ident_Plot.STA_SD(:,:,i)); 
            hold(loop_plot,'on');
            if add_info.settings.kernel_new.Heat_colour
                colormap(loop_plot,colorMap_arr(:,:,i));
                if ~isempty(RF_Ident_Plot.STASD_ASP_Num_RF_pixels{i})
                    scatter(loop_plot,RF_Ident_Plot.STASD_ASP_RF_coords{i}(:,2),RF_Ident_Plot.STASD_ASP_RF_coords{i}(:,1),'wo','LineWidth',1.5);
                end
            else % Heat_Map_Colour_Choice == 1
                if ~isempty(RF_Ident_Plot.STASD_ASP_Num_RF_pixels{i})
                    scatter(loop_plot,RF_Ident_Plot.STASD_ASP_RF_coords{i}(:,2),RF_Ident_Plot.STASD_ASP_RF_coords{i}(:,1),'go','LineWidth',1.5);
                end
                
                if add_info.settings.kernel_new.Heat_gray
                    colormap(loop_plot,gray(256));
                end
                if add_info.settings.kernel_new.Heat_colour
                    colormap(loop_plot,colorMap_arr(:,:,i));
                end
             
             
            end
            axis(loop_plot,'equal'); axis(loop_plot,'tight');
            if i>2
                xlabel(loop_plot,'x');
            end
            if i==1||i==3
                ylabel(loop_plot,'y');
            end
            title(loop_plot,Spectral_Names{i});
            %if i == 4
            original_ax = get(loop_plot,'Position');
            if add_info.settings.kernel_new.plot_colorbar
                colorbar(loop_plot,'southoutside');
            end
            set(loop_plot,'Position',original_ax);
            %end
            set(loop_plot,'FontSize',12);
        end
    end
    hold(loop_plot,'off');
    
    if strcmp(get(fig2,'type'),'uipanel')
         sgtitle(fig2,plot_folder_legend);
    else
        fig2.Name = plot_folder_legend;
        fig2.NumberTitle = 'off';
        sgtitle(fig2,'All sig. Pixel');

    end
%     annotation(fig2,'textbox',[.5 .975 0 0],'String',plot_folder_legend,...
%         'FitBoxToText','on','FontSize',12,'FontWeight','bold','LineStyle',...
%         'none','HorizontalAlignment','center','VerticalAlignment','middle');
    %set(gcf,'color','w');

end

if add_info.settings.kernel_new.plot_Gauss % Gaussian

     stim_rows = size(RF_Ident_Plot.STA_SD,1);
     stim_columns = size(RF_Ident_Plot.STA_SD,2);

    if size(add_info.settings.kernel_new.folderplot,2)*...
            length(add_info.settings.kernel_new.cell_to_plot) == 1 ...
            && ~add_info.settings.kernel_new.plot_in_new_window
       fig3 = add_info.panels.GaussRF;
    else
       fig3 = figure;
    end
    
    if Spectral_Dim == 1
        colormap(fig3,gray(256));
%                 imagesc(RF_Ident_Plot.STA_SD,[min(RF_Ident_Plot.STA_SD,[],'all') max(RF_Ident_Plot.STA_SD,[],'all')]); hold on;
%                 scatter(RF_Ident_Plot.STASD_Gaus_RF_coords(:,2),RF_Ident_Plot.STASD_Gaus_RF_coords(:,1),'go','LineWidth',1.5);
%                 plot(RF_Ident_Plot.STASD_Gaus_Gaussian_mean(1),RF_Ident_Plot.STASD_Gaus_Gaussian_mean(2),'rx','LineWidth',1.5); %
%                 fcontour(RF_Ident_Plot.STASD_Gaus_gmPDF,[0.5 (p.stim_columns+0.5) 0.5 (p.stim_rows+0.5)],'r','LineWidth',1.5,'LevelList',RF_Ident_Plot.Thresh_height); % 'LineColor',LightestBlue,'LineWidth',2
        axis equal; axis tight;
        xlabel('x');
        ylabel('y');
        title('RF');
        if add_info.settings.kernel_new.plot_colorbar
            colorbar(loop_plot,'southoutside');
        end
        set(gca,'FontSize',12);
    elseif Spectral_Dim == 4
        if add_info.settings.kernel_new.Heat_gray
            colormap(loop_plot,gray(256));
        end
        %colormap(loop_plot,gray(256));
        for i = 1:Spectral_Dim
            loop_plot = subplot(2,2,i,'parent',fig3);
            imagesc(loop_plot,RF_Ident_Plot.STA_SD(:,:,i),...
                [min(RF_Ident_Plot.STA_SD(:,:,i),[],'all')...
                max(RF_Ident_Plot.STA_SD(:,:,i),[],'all')]);
            hold(loop_plot,'on')
            if add_info.settings.kernel_new.Heat_colour
                colormap(loop_plot,colorMap_arr(:,:,i));
                if ~isempty(RF_Ident_Plot.STASD_Gaus_Num_RF_pixels{i})
                    scatter(loop_plot,RF_Ident_Plot.STASD_Gaus_RF_coords{i}(:,2),...
                        RF_Ident_Plot.STASD_Gaus_RF_coords{i}(:,1),'wo','LineWidth',1.5);
                    
                    plot(loop_plot,RF_Ident_Plot.STASD_Gaus_Gaussian_mean{i}(1),...
                        RF_Ident_Plot.STASD_Gaus_Gaussian_mean{i}(2),'wx','LineWidth',1.5);
                    
                    fcontour(loop_plot,RF_Ident_Plot.STASD_Gaus_gmPDF{i},...
                        [0.5 (stim_columns+0.5) 0.5 (stim_rows+0.5)],'w',...
                        'LineWidth',1.5,'LevelList',RF_Ident_Plot.STASD_Gaus_Thresh_height{i});
                end
            else % Heat_Map_Colour_Choice == 1
                if ~isempty(RF_Ident_Plot.STASD_Gaus_Num_RF_pixels{i})
                    scatter(loop_plot,RF_Ident_Plot.STASD_Gaus_RF_coords{i}(:,2),...
                        RF_Ident_Plot.STASD_Gaus_RF_coords{i}(:,1),'go','LineWidth',1.5);
                    plot(loop_plot,RF_Ident_Plot.STASD_Gaus_Gaussian_mean{i}(1),...
                        RF_Ident_Plot.STASD_Gaus_Gaussian_mean{i}(2),'rx','LineWidth',1.5);
                    fcontour(loop_plot,RF_Ident_Plot.STASD_Gaus_gmPDF{i},...
                        [0.5 (stim_columns+0.5) 0.5 (stim_rows+0.5)],'r',...
                        'LineWidth',1.5,'LevelList',RF_Ident_Plot.STASD_Gaus_Thresh_height{i});
                end
                if add_info.settings.kernel_new.Heat_gray
                   colormap(loop_plot,gray(256));
                end
            end
            axis(loop_plot,'equal'); axis(loop_plot,'tight');
            if i>2
                xlabel(loop_plot,'x');
            end
            if i==1||i==3
                ylabel(loop_plot,'y');
            end
            title(loop_plot,Spectral_Names{i});
            %if i == 4
            original_ax = get(loop_plot,'Position');
            if add_info.settings.kernel_new.plot_colorbar
                colorbar(loop_plot,'southoutside');
            end
            set(loop_plot,'Position',original_ax);
            %end
            set(loop_plot,'FontSize',12);
        end
    end
    hold(loop_plot,'off');
    if strcmp(get(fig3,'type'),'uipanel')
         sgtitle(fig3,plot_folder_legend);
    else
            fig3.Name = plot_folder_legend;
            fig3.NumberTitle = 'off';
            sgtitle(fig3,'Gaussian');

    end
%     annotation(fig3,'textbox',[.5 .975 0 0],'String',plot_folder_legend,...
%         'FitBoxToText','on','FontSize',12,'FontWeight','bold','LineStyle',...
%         'none','HorizontalAlignment','center','VerticalAlignment','middle');
    

end
%% SC Method

%PLotting for the SC method

elseif RF_type == 2
   if add_info.settings.kernel_new.plot_Box  % Box
        % Organize figures
    if size(add_info.settings.kernel_new.folderplot,2)*...
            length(add_info.settings.kernel_new.cell_to_plot) == 1 ...
            && ~add_info.settings.kernel_new.plot_in_new_window
       fig1 = add_info.panels.BoxRF;
    else
       fig1 = figure;
    end
    %fig1 = figure;
    if Spectral_Dim == 1 %has to be not hard coded
        colormap(fig1,gray(256));
%                 imagesc(RF_Ident_Plot.STA_SD); hold on;
%                 scatter(RF_Ident_Plot.RF_coords_centre(2),RF_Ident_Plot.RF_coords_centre(1),'rx','LineWidth',1.5);
%                 scatter(RF_Ident_Plot.RF_coords_noncentre(:,2),RF_Ident_Plot.RF_coords_noncentre(:,1),'go','LineWidth',1.5);
        axis(loop_plot,'equal'); axis(loop_plot,'tight');
        xlabel('x');
        ylabel('y');
        title('RF');
        if add_info.settings.kernel_new.plot_colorbar
            colorbar(loop_plot,'southoutside');
        end
        set(gca,'FontSize',12);
    elseif Spectral_Dim == 4
%         if add_info.settings.kernel_new.Heat_gray
%             colormap(add_info.panels.BoxRF,gray(256));
%         end
        for i = 1:Spectral_Dim
            
            loop_plot = subplot(2,2,i,'parent',fig1);
            imagesc(loop_plot,RF_Ident_Plot.SCA_Pixel_covar(:,:,i));
            hold(loop_plot,'on');
            if add_info.settings.kernel_new.Heat_colour
                colormap(loop_plot,colorMap_arr(:,:,i));
            end
             if add_info.settings.kernel_new.Heat_gray
                colormap(loop_plot,gray(256));
             end
            if ~isempty(RF_Ident_Plot.SC_Box_Num_RF_pixels{i})
                scatter(loop_plot,RF_Ident_Plot.SC_Box_RF_coords_centre{i}(2),...
                    RF_Ident_Plot.SC_Box_RF_coords_centre{i}(1),'rx','LineWidth',1.5,...
                    'MarkerFaceAlpha',0.1);
                scatter(loop_plot,RF_Ident_Plot.SC_Box_RF_coords_noncentre{i}(:,2)...
                    ,RF_Ident_Plot.SC_Box_RF_coords_noncentre{i}(:,1),'go','LineWidth',1.5,...
                    'MarkerFaceAlpha',0.1);
            end
            axis(loop_plot,'equal'); axis(loop_plot, 'tight');
            if i>2
                xlabel(loop_plot,'x');
            end
            if i==1||i==3
                ylabel(loop_plot,'y');
            end
            %title(Spectral_Names{i});
            %if i == 4
            original_ax = get(loop_plot,'Position');
            if add_info.settings.kernel_new.plot_colorbar
                colorbar(loop_plot,'southoutside');
            end
            set(loop_plot,'Position',original_ax);
            %end
            title(loop_plot,Spectral_Names{i});
            loop_plot.TitleHorizontalAlignment = 'left';
            set(loop_plot,'FontSize',12);
        end
       
    end
    
    if strcmp(get(fig1,'type'),'uipanel')
         sgtitle(fig1,plot_folder_legend);
    else
            fig1.Name = plot_folder_legend;
            fig1.NumberTitle = 'off';

    end
%     annotation(fig1,'textbox',[.5 .975 0 0],'String',plot_folder_legend,...
%         'FitBoxToText','on','FontSize',12,'FontWeight','bold','LineStyle',...
%         'none','HorizontalAlignment','center','VerticalAlignment','middle');
    %set(gcf,'color','w');

end

if add_info.settings.kernel_new.plot_Allpixel  % All Significant Pixels

   if size(add_info.settings.kernel_new.folderplot,2)*...
            length(add_info.settings.kernel_new.cell_to_plot) == 1 ...
            && ~add_info.settings.kernel_new.plot_in_new_window
       fig2 = add_info.panels.AllPxRF;
    else
       fig2 = figure;
    end
    
    
    if Spectral_Dim == 1
        colormap(fig2,gray(256));
%                 imagesc(RF_Ident_Plot.STA_SD); hold on;
%                 scatter(RF_Ident_Plot.RF_coords(:,2),RF_Ident_Plot.RF_coords(:,1),'go','LineWidth',1.5);
        axis(loop_plot,'equal'); axis(loop_plot,'tight');
        xlabel('x');
        ylabel('y');
        title('RF');
        if add_info.settings.kernel_new.plot_colorbar
            colorbar(loop_plot,'southoutside');
        end
        set(gca,'FontSize',12);
    elseif Spectral_Dim == 4
        
        for i = 1:Spectral_Dim
            loop_plot = subplot(2,2,i,'parent',fig2);
            imagesc(loop_plot,RF_Ident_Plot.SCA_Pixel_covar(:,:,i)); 
            hold(loop_plot,'on');
            if add_info.settings.kernel_new.Heat_colour
                colormap(loop_plot,colorMap_arr(:,:,i));
                if ~isempty(RF_Ident_Plot.SC_ASP_Num_RF_pixels{i})
                    scatter(loop_plot,RF_Ident_Plot.SC_ASP_RF_coords{i}(:,2),RF_Ident_Plot.SC_ASP_RF_coords{i}(:,1),'wo','LineWidth',1.5);
                end
            else % Heat_Map_Colour_Choice == 1
                if ~isempty(RF_Ident_Plot.SC_Gaus_Num_RF_pixels{i})
                    scatter(loop_plot,RF_Ident_Plot.SC_ASP_RF_coords{i}(:,2),RF_Ident_Plot.SC_ASP_RF_coords{i}(:,1),'go','LineWidth',1.5);
                end
                
                if add_info.settings.kernel_new.Heat_gray
                    colormap(loop_plot,gray(256));
                end
                
                if add_info.settings.kernel_new.Heat_colour
                    colormap(loop_plot,colorMap_arr(:,:,i));
                end
             
            end
            axis(loop_plot,'equal'); axis(loop_plot,'tight');
            if i>2
                xlabel(loop_plot,'x');
            end
            if i==1||i==3
                ylabel(loop_plot,'y');
            end
            title(loop_plot,Spectral_Names{i});
            loop_plot.TitleHorizontalAlignment = 'left';
            %if i == 4
            original_ax = get(loop_plot,'Position');
            if add_info.settings.kernel_new.plot_colorbar
               colorbar(loop_plot,'southoutside');
            end
            set(loop_plot,'Position',original_ax);
            %end
            set(loop_plot,'FontSize',12);
        end
    end
    
    
    if strcmp(get(fig2,'type'),'uipanel')
         sgtitle(fig2,plot_folder_legend);
    else
            fig2.Name = plot_folder_legend;
            fig2.NumberTitle = 'off';

    end
%     annotation(fig2,'textbox',[.5 .975 0 0],'String',plot_folder_legend,...
%         'FitBoxToText','on','FontSize',12,'FontWeight','bold','LineStyle',...
%         'none','HorizontalAlignment','center','VerticalAlignment','middle');
    %set(gcf,'color','w');

end

if add_info.settings.kernel_new.plot_Gauss % Gaussian

     stim_rows = size(RF_Ident_Plot.SCA_Pixel_covar,1);
     stim_columns = size(RF_Ident_Plot.SCA_Pixel_covar,2);

    if size(add_info.settings.kernel_new.folderplot,2)*...
            length(add_info.settings.kernel_new.cell_to_plot) == 1 ...
            && ~add_info.settings.kernel_new.plot_in_new_window
       fig3 = add_info.panels.GaussRF;
    else
       fig3 = figure;
    end
    
    if Spectral_Dim == 1
        colormap(fig3,gray(256));
%                 imagesc(RF_Ident_Plot.STA_SD,[min(RF_Ident_Plot.STA_SD,[],'all') max(RF_Ident_Plot.STA_SD,[],'all')]); hold on;
%                 scatter(RF_Ident_Plot.STASD_Gaus_RF_coords(:,2),RF_Ident_Plot.STASD_Gaus_RF_coords(:,1),'go','LineWidth',1.5);
%                 plot(RF_Ident_Plot.STASD_Gaus_Gaussian_mean(1),RF_Ident_Plot.STASD_Gaus_Gaussian_mean(2),'rx','LineWidth',1.5); %
%                 fcontour(RF_Ident_Plot.STASD_Gaus_gmPDF,[0.5 (p.stim_columns+0.5) 0.5 (p.stim_rows+0.5)],'r','LineWidth',1.5,'LevelList',RF_Ident_Plot.Thresh_height); % 'LineColor',LightestBlue,'LineWidth',2
        axis(loop_plot,'equal'); axis(loop_plot,'tight');
        xlabel('x');
        ylabel('y');
        title('RF');
        if add_info.settings.kernel_new.plot_colorbar
            colorbar(loop_plot,'southoutside');
        end
        set(gca,'FontSize',12);
    elseif Spectral_Dim == 4
        if add_info.settings.kernel_new.Heat_gray
            colormap(loop_plot,gray(256));
        end
       
        for i = 1:Spectral_Dim
            loop_plot = subplot(2,2,i,'parent',fig3);
            imagesc(loop_plot,RF_Ident_Plot.SCA_Pixel_covar(:,:,i),...
                [min(RF_Ident_Plot.SCA_Pixel_covar(:,:,i),[],'all')...
                max(RF_Ident_Plot.SCA_Pixel_covar(:,:,i),[],'all')]);
            hold(loop_plot,'on')
            if add_info.settings.kernel_new.Heat_colour
                colormap(loop_plot,colorMap_arr(:,:,i));
                if ~isempty(RF_Ident_Plot.SC_Gaus_Num_RF_pixels{i})
                    scatter(loop_plot,RF_Ident_Plot.SC_Gaus_RF_coords{i}(:,2),...
                        RF_Ident_Plot.SC_Gaus_RF_coords{i}(:,1),'wo','LineWidth',1.5);
                    
                    plot(loop_plot,RF_Ident_Plot.SC_Gaus_Gaussian_mean{i}(1),...
                        RF_Ident_Plot.SC_Gaus_Gaussian_mean{i}(2),'wx','LineWidth',1.5);
                    
                    fcontour(loop_plot,RF_Ident_Plot.SC_Gaus_gmPDF{i},...
                        [0.5 (stim_columns+0.5) 0.5 (stim_rows+0.5)],'w',...
                        'LineWidth',1.5,'LevelList',RF_Ident_Plot.SC_Gaus_Thresh_height{i});
                end
            else % Heat_Map_Colour_Choice == 1
                if ~isempty(RF_Ident_Plot.SC_Gaus_Num_RF_pixels{i})
                    scatter(loop_plot,RF_Ident_Plot.SC_Gaus_RF_coords{i}(:,2),...
                        RF_Ident_Plot.SC_Gaus_RF_coords{i}(:,1),'go','LineWidth',1.5);
                    plot(loop_plot,RF_Ident_Plot.SC_Gaus_Gaussian_mean{i}(1),...
                        RF_Ident_Plot.SC_Gaus_Gaussian_mean{i}(2),'rx','LineWidth',1.5);
                    fcontour(loop_plot,RF_Ident_Plot.SC_Gaus_gmPDF{i},...
                        [0.5 (stim_columns+0.5) 0.5 (stim_rows+0.5)],'r',...
                        'LineWidth',1.5,'LevelList',RF_Ident_Plot.SC_Gaus_Thresh_height{i});
                end
                if add_info.settings.kernel_new.Heat_gray
                   colormap(loop_plot,gray(256));
                end
            end
            axis(loop_plot,'equal'); axis(loop_plot,'tight');
            if i>2
                xlabel(loop_plot,'x');
            end
            if i==1||i==3
                ylabel(loop_plot,'y');
            end
            title(loop_plot,Spectral_Names{i});
            loop_plot.TitleHorizontalAlignment = 'left';
            %if i == 4
            original_ax = get(loop_plot,'Position');
            if add_info.settings.kernel_new.plot_colorbar
                colorbar(loop_plot,'southoutside');
            end
            set(loop_plot,'Position',original_ax);
            %end
            set(loop_plot,'FontSize',12);
        end
    end
    
    if strcmp(get(fig3,'type'),'uipanel')
         sgtitle(fig3,plot_folder_legend);
    else
            fig3.Name = plot_folder_legend;
            fig3.NumberTitle = 'off';

    end
    

end
    
    
    
    
    
    
end
%catch
   % warning('No receptive field data found for given cell and threshold');
%end
if add_info.settings.kernel_new.plot_all_STA || add_info.settings.kernel_new.plot_STA_sig_pix  
    Plot_STA_app(RF_Ident_Plot,add_info.settings,RF_type);
end
end
end

%% Enhance axes interactions and performance
% Every second child object in the panel is an axis







%% Plot STA if wanted





if ~strcmp(get(fig3,'type'),'uipanel')
    autoArrangeFigures
end
out = 1;

end



