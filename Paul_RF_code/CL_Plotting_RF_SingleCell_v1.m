%% Plot RF and Associated Summary Statistics - Single Cell
function out = CL_Plotting_RF_SingleCell_v1 (savepath, add_info)
% Choose which cell RF to plot
cell_to_plot  = add_info.settings.kernel_new.cell_to_plot;
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
       
    end
    
    
    if strcmp(get(fig1,'type'),'uipanel')
         sgtitle(fig1,plot_folder_legend);
    else
            fig.Name = plot_folder_legend;
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
        colormap(loop_plot,gray(256));
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
    
    if strcmp(get(fig3,'type'),'uipanel')
         sgtitle(fig3,plot_folder_legend);
    else
            fig3.Name = plot_folder_legend;
            fig3.NumberTitle = 'off';

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


%     %% RF Plots
%     
%     if add_info.settings.kernel_new.singlecplot.SS_STA % STA-SD method;
%         %Load the data
%         %First, we load the overview variable, as much smaller and we can get 
%         %the index in the structure than we index directly into the .matfile.
%         L = load(findfile_app(add_info.stim_idx,savepath,"RF_Ident_SS_Std.mat"),"RF_overview");
%         index = find([L.RF_overview.cell_idx] == cell_to_plot);
%         M = matfile(findfile_app(add_info.stim_idx,savepath,"RF_Ident_SS_Std.mat"),...
%             'Writable',true);
%         cell_file = M.RF_Ident_SS_Std(1,index);
%         RF_Ident_Plot = cell_file.RF_results;
%        
%         
%         if add_info.settings.kernel_new.singlecplot.Box  % Box
%             
%             fig1 = figure;
%             if Spectral_Dim == 1 %has to be not hard coded
%                 colormap(fig1,gray(256));
% %                 imagesc(RF_Ident_Plot.STA_SD); hold on;
% %                 scatter(RF_Ident_Plot.RF_coords_centre(2),RF_Ident_Plot.RF_coords_centre(1),'rx','LineWidth',1.5);
% %                 scatter(RF_Ident_Plot.RF_coords_noncentre(:,2),RF_Ident_Plot.RF_coords_noncentre(:,1),'go','LineWidth',1.5);
%                 axis equal; axis tight;
%                 xlabel('x');
%                 ylabel('y');
%                 title('RF');
%                 colorbar;
%                 set(gca,'FontSize',12);
%             elseif Spectral_Dim == 4
%                 if add_info.settings.kernel_new.singlecplot.Heat_gray
%                     colormap(fig1,gray(256));
%                 end
%                 for i = 1:Spectral_Dim
%                     loop_plot = subplot(2,2,i);
%                     imagesc(RF_Ident_Plot.STA_SD(:,:,i)); hold on;
%                     if add_info.settings.kernel_new.singlecplot.Heat_colour
%                         colormap(loop_plot,colorMap_arr(:,:,i));
%                     end
%                     if ~isempty(RF_Ident_Plot.STASD_Box_Num_RF_pixels{i})
%                         scatter(RF_Ident_Plot.STASD_Box_RF_coords_centre{i}(2),RF_Ident_Plot.STASD_Box_RF_coords_centre{i}(1),'rx','LineWidth',1.5);
%                         scatter(RF_Ident_Plot.STASD_Box_RF_coords_noncentre{i}(:,2),RF_Ident_Plot.STASD_Box_RF_coords_noncentre{i}(:,1),'go','LineWidth',1.5);
%                     end
%                     axis equal; axis tight;
%                     if i>2
%                         xlabel('x');
%                     end
%                     if i==1||i==3
%                         ylabel('y');
%                     end
%                     %title(Spectral_Names{i});
%                     colorbar;
%                     set(gca,'FontSize',12);
%                 end
%             end
%             annotation('textbox',[.5 .975 0 0],'String','STA-SD -- Box RF','FitBoxToText','on','FontSize',12,'FontWeight','bold','LineStyle','none','HorizontalAlignment','center','VerticalAlignment','middle');
%             set(gcf,'color','w');
%             
%         end
%         
%         if add_info.settings.kernel_new.singlecplot.Allpix  % All Significant Pixels
%             
%             fig2 = figure;
%             if Spectral_Dim == 1
%                 colormap(fig2,gray(256));
% %                 imagesc(RF_Ident_Plot.STA_SD); hold on;
% %                 scatter(RF_Ident_Plot.RF_coords(:,2),RF_Ident_Plot.RF_coords(:,1),'go','LineWidth',1.5);
%                 axis equal; axis tight;
%                 xlabel('x');
%                 ylabel('y');
%                 title('RF');
%                 colorbar;
%                 set(gca,'FontSize',12);
%             elseif Spectral_Dim == 4
%                 if add_info.settings.kernel_new.singlecplot.Heat_gray
%                     colormap(fig2,gray(256));
%                 end
%                 for i = 1:Spectral_Dim
%                     loop_plot = subplot(2,2,i);
%                     imagesc(RF_Ident_Plot.STA_SD(:,:,i)); hold on;
%                     if add_info.settings.kernel_new.singlecplot.Heat_colour
%                         colormap(loop_plot,colorMap_arr(:,:,i));
%                         if ~isempty(RF_Ident_Plot.STASD_ASP_Num_RF_pixels{i})
%                             scatter(RF_Ident_Plot.STASD_ASP_RF_coords{i}(:,2),RF_Ident_Plot.STASD_ASP_RF_coords{i}(:,1),'wo','LineWidth',1.5);
%                         end
%                     else % Heat_Map_Colour_Choice == 1
%                         if ~isempty(RF_Ident_Plot.STASD_ASP_Num_RF_pixels{i})
%                             scatter(RF_Ident_Plot.STASD_ASP_RF_coords{i}(:,2),RF_Ident_Plot.STASD_ASP_RF_coords{i}(:,1),'go','LineWidth',1.5);
%                         end
%                     end
%                     axis equal; axis tight;
%                     if i>2
%                         xlabel('x');
%                     end
%                     if i==1||i==3
%                         ylabel('y');
%                     end
%                     title(Spectral_Names{i});
%                     colorbar;
%                     set(gca,'FontSize',12);
%                 end
%             end
%             annotation('textbox',[.5 .975 0 0],'String','STA-SD -- All Sig. Pix. RF','FitBoxToText','on','FontSize',12,'FontWeight','bold','LineStyle','none','HorizontalAlignment','center','VerticalAlignment','middle');
%             set(gcf,'color','w');
%             
%         end
%         
%         if add_info.settings.kernel_new.singlecplot.Gaussian % Gaussian
%             
%              stim_rows = size(RF_Ident_Plot.STA_SD,1);
%              stim_columns = size(RF_Ident_Plot.STA_SD,2);
%             
%             fig3 = figure;
%             if Spectral_Dim == 1
%                 colormap(fig3,gray(256));
% %                 imagesc(RF_Ident_Plot.STA_SD,[min(RF_Ident_Plot.STA_SD,[],'all') max(RF_Ident_Plot.STA_SD,[],'all')]); hold on;
% %                 scatter(RF_Ident_Plot.STASD_Gaus_RF_coords(:,2),RF_Ident_Plot.STASD_Gaus_RF_coords(:,1),'go','LineWidth',1.5);
% %                 plot(RF_Ident_Plot.STASD_Gaus_Gaussian_mean(1),RF_Ident_Plot.STASD_Gaus_Gaussian_mean(2),'rx','LineWidth',1.5); %
% %                 fcontour(RF_Ident_Plot.STASD_Gaus_gmPDF,[0.5 (p.stim_columns+0.5) 0.5 (p.stim_rows+0.5)],'r','LineWidth',1.5,'LevelList',RF_Ident_Plot.Thresh_height); % 'LineColor',LightestBlue,'LineWidth',2
%                 axis equal; axis tight;
%                 xlabel('x');
%                 ylabel('y');
%                 title('RF');
%                 colorbar;
%                 set(gca,'FontSize',12);
%             elseif Spectral_Dim == 4
%                 if add_info.settings.kernel_new.singlecplot.Heat_gray
%                     colormap(fig3,gray(256));
%                 end
%                 for i = 1:Spectral_Dim
%                     loop_plot = subplot(2,2,i);
%                     imagesc(RF_Ident_Plot.STA_SD(:,:,i),[min(RF_Ident_Plot.STA_SD(:,:,i),[],'all') max(RF_Ident_Plot.STA_SD(:,:,i),[],'all')]); hold on;
%                     if add_info.settings.kernel_new.singlecplot.Heat_colour
%                         colormap(loop_plot,colorMap_arr(:,:,i));
%                         if ~isempty(RF_Ident_Plot.STASD_Gaus_Num_RF_pixels{i})
%                             scatter(RF_Ident_Plot.STASD_Gaus_RF_coords{i}(:,2),RF_Ident_Plot.STASD_Gaus_RF_coords{i}(:,1),'wo','LineWidth',1.5);
%                             plot(RF_Ident_Plot.STASD_Gaus_Gaussian_mean{i}(1),RF_Ident_Plot.STASD_Gaus_Gaussian_mean{i}(2),'wx','LineWidth',1.5);
%                             fcontour(RF_Ident_Plot.STASD_Gaus_gmPDF{i},[0.5 (stim_columns+0.5) 0.5 (stim_rows+0.5)],'w','LineWidth',1.5,'LevelList',RF_Ident_Plot.STASD_Gaus_Thresh_height{i});
%                         end
%                     else % Heat_Map_Colour_Choice == 1
%                         if ~isempty(RF_Ident_Plot.STASD_Gaus_Num_RF_pixels{i})
%                             scatter(RF_Ident_Plot.STASD_Gaus_RF_coords{i}(:,2),RF_Ident_Plot.STASD_Gaus_RF_coords{i}(:,1),'go','LineWidth',1.5);
%                             plot(RF_Ident_Plot.STASD_Gaus_Gaussian_mean{i}(1),RF_Ident_Plot.STASD_Gaus_Gaussian_mean{i}(2),'rx','LineWidth',1.5);
%                             fcontour(RF_Ident_Plot.STASD_Gaus_gmPDF{i},[0.5 (stim_columns+0.5) 0.5 (stim_rows+0.5)],'r','LineWidth',1.5,'LevelList',RF_Ident_Plot.STASD_Gaus_Thresh_height{i});
%                         end
%                     end
%                     axis equal; axis tight;
%                     if i>2
%                         xlabel('x');
%                     end
%                     if i==1||i==3
%                         ylabel('y');
%                     end
%                     title(Spectral_Names{i});
%                     colorbar;
%                     set(gca,'FontSize',12);
%                 end
%             end
%             annotation('textbox',[.5 .975 0 0],'String','STA-SD -- Gaussian RF','FitBoxToText','on','FontSize',12,'FontWeight','bold','LineStyle','none','HorizontalAlignment','center','VerticalAlignment','middle');
%             set(gcf,'color','w');
%             
%         end
%         
%     end
%     
%     
%     if add_info.settings.kernel_new.singlecplot.SS_CI % STA-SD method;
%         %Load the data
%         %First, we load the overview variable, as much smaller and we can get 
%         %the index in the structure than we index directly into the .matfile.
%         L = load(findfile_app(add_info.stim_idx,savepath,"RF_Ident_SS_CI.mat"),"RF_overview");
%         index = find([L.RF_overview.cell_idx] == cell_to_plot);
%         M = matfile(findfile_app(add_info.stim_idx,savepath,"RF_Ident_SS_CI.mat"),...
%             'Writable',true);
%         cell_file = M.RF_Ident_SS_CI(1,index);
%         RF_Ident_Plot = cell_file.RF_results;
%        
%         
%         if add_info.settings.kernel_new.singlecplot.Box  % Box
%             
%             fig1 = figure;
%             if Spectral_Dim == 1 %has to be not hard coded
%                 colormap(fig1,gray(256));
% %                 imagesc(RF_Ident_Plot.STA_SD); hold on;
% %                 scatter(RF_Ident_Plot.RF_coords_centre(2),RF_Ident_Plot.RF_coords_centre(1),'rx','LineWidth',1.5);
% %                 scatter(RF_Ident_Plot.RF_coords_noncentre(:,2),RF_Ident_Plot.RF_coords_noncentre(:,1),'go','LineWidth',1.5);
%                 axis equal; axis tight;
%                 xlabel('x');
%                 ylabel('y');
%                 title('RF');
%                 colorbar;
%                 set(gca,'FontSize',12);
%             elseif Spectral_Dim == 4
%                 if add_info.settings.kernel_new.singlecplot.Heat_gray
%                     colormap(fig1,gray(256));
%                 end
%                 for i = 1:Spectral_Dim
%                     loop_plot = subplot(2,2,i);
%                     imagesc(RF_Ident_Plot.STA_SD(:,:,i)); hold on;
%                     if add_info.settings.kernel_new.singlecplot.Heat_colour
%                         colormap(loop_plot,colorMap_arr(:,:,i));
%                     end
%                     if ~isempty(RF_Ident_Plot.STASD_Box_Num_RF_pixels{i})
%                         scatter(RF_Ident_Plot.STASD_Box_RF_coords_centre{i}(2),RF_Ident_Plot.STASD_Box_RF_coords_centre{i}(1),'rx','LineWidth',1.5);
%                         scatter(RF_Ident_Plot.STASD_Box_RF_coords_noncentre{i}(:,2),RF_Ident_Plot.STASD_Box_RF_coords_noncentre{i}(:,1),'go','LineWidth',1.5);
%                     end
%                     axis equal; axis tight;
%                     if i>2
%                         xlabel('x');
%                     end
%                     if i==1||i==3
%                         ylabel('y');
%                     end
%                     %title(Spectral_Names{i});
%                     colorbar;
%                     set(gca,'FontSize',12);
%                 end
%             end
%             annotation('textbox',[.5 .975 0 0],'String','STA-SD CI -- Box RF','FitBoxToText','on','FontSize',12,'FontWeight','bold','LineStyle','none','HorizontalAlignment','center','VerticalAlignment','middle');
%             set(gcf,'color','w');
%             
%         end
%         
%         if add_info.settings.kernel_new.singlecplot.Allpix  % All Significant Pixels
%             
%             fig2 = figure;
%             if Spectral_Dim == 1
%                 colormap(fig2,gray(256));
% %                 imagesc(RF_Ident_Plot.STA_SD); hold on;
% %                 scatter(RF_Ident_Plot.RF_coords(:,2),RF_Ident_Plot.RF_coords(:,1),'go','LineWidth',1.5);
%                 axis equal; axis tight;
%                 xlabel('x');
%                 ylabel('y');
%                 title('RF');
%                 colorbar;
%                 set(gca,'FontSize',12);
%             elseif Spectral_Dim == 4
%                 if add_info.settings.kernel_new.singlecplot.Heat_gray
%                     colormap(fig2,gray(256));
%                 end
%                 for i = 1:Spectral_Dim
%                     loop_plot = subplot(2,2,i);
%                     imagesc(RF_Ident_Plot.STA_SD(:,:,i)); hold on;
%                     if add_info.settings.kernel_new.singlecplot.Heat_colour
%                         colormap(loop_plot,colorMap_arr(:,:,i));
%                         if ~isempty(RF_Ident_Plot.STASD_ASP_Num_RF_pixels{i})
%                             scatter(RF_Ident_Plot.STASD_ASP_RF_coords{i}(:,2),RF_Ident_Plot.STASD_ASP_RF_coords{i}(:,1),'wo','LineWidth',1.5);
%                         end
%                     else % Heat_Map_Colour_Choice == 1
%                         if ~isempty(RF_Ident_Plot.STASD_ASP_Num_RF_pixels{i})
%                             scatter(RF_Ident_Plot.STASD_ASP_RF_coords{i}(:,2),RF_Ident_Plot.STASD_ASP_RF_coords{i}(:,1),'go','LineWidth',1.5);
%                         end
%                     end
%                     axis equal; axis tight;
%                     if i>2
%                         xlabel('x');
%                     end
%                     if i==1||i==3
%                         ylabel('y');
%                     end
%                     title(Spectral_Names{i});
%                     colorbar;
%                     set(gca,'FontSize',12);
%                 end
%             end
%             annotation('textbox',[.5 .975 0 0],'String','STA-SD CI-- All Sig. Pix. RF','FitBoxToText','on','FontSize',12,'FontWeight','bold','LineStyle','none','HorizontalAlignment','center','VerticalAlignment','middle');
%             set(gcf,'color','w');
%             
%         end
%         
%         if add_info.settings.kernel_new.singlecplot.Gaussian % Gaussian
%             
%              stim_rows = size(RF_Ident_Plot.STA_SD,1);
%              stim_columns = size(RF_Ident_Plot.STA_SD,2);
%             
%             fig3 = figure;
%             if Spectral_Dim == 1
%                 colormap(fig3,gray(256));
% %                 imagesc(RF_Ident_Plot.STA_SD,[min(RF_Ident_Plot.STA_SD,[],'all') max(RF_Ident_Plot.STA_SD,[],'all')]); hold on;
% %                 scatter(RF_Ident_Plot.STASD_Gaus_RF_coords(:,2),RF_Ident_Plot.STASD_Gaus_RF_coords(:,1),'go','LineWidth',1.5);
% %                 plot(RF_Ident_Plot.STASD_Gaus_Gaussian_mean(1),RF_Ident_Plot.STASD_Gaus_Gaussian_mean(2),'rx','LineWidth',1.5); %
% %                 fcontour(RF_Ident_Plot.STASD_Gaus_gmPDF,[0.5 (p.stim_columns+0.5) 0.5 (p.stim_rows+0.5)],'r','LineWidth',1.5,'LevelList',RF_Ident_Plot.Thresh_height); % 'LineColor',LightestBlue,'LineWidth',2
%                 axis equal; axis tight;
%                 xlabel('x');
%                 ylabel('y');
%                 title('RF');
%                 colorbar;
%                 set(gca,'FontSize',12);
%             elseif Spectral_Dim == 4
%                 if add_info.settings.kernel_new.singlecplot.Heat_gray
%                     colormap(fig3,gray(256));
%                 end
%                 for i = 1:Spectral_Dim
%                     loop_plot = subplot(2,2,i);
%                     imagesc(RF_Ident_Plot.STA_SD(:,:,i),[min(RF_Ident_Plot.STA_SD(:,:,i),[],'all') max(RF_Ident_Plot.STA_SD(:,:,i),[],'all')]); hold on;
%                     if add_info.settings.kernel_new.singlecplot.Heat_colour
%                         colormap(loop_plot,colorMap_arr(:,:,i));
%                         if ~isempty(RF_Ident_Plot.STASD_Gaus_Num_RF_pixels{i})
%                             scatter(RF_Ident_Plot.STASD_Gaus_RF_coords{i}(:,2),RF_Ident_Plot.STASD_Gaus_RF_coords{i}(:,1),'wo','LineWidth',1.5);
%                             plot(RF_Ident_Plot.STASD_Gaus_Gaussian_mean{i}(1),RF_Ident_Plot.STASD_Gaus_Gaussian_mean{i}(2),'wx','LineWidth',1.5);
%                             fcontour(RF_Ident_Plot.STASD_Gaus_gmPDF{i},[0.5 (stim_columns+0.5) 0.5 (stim_rows+0.5)],'w','LineWidth',1.5,'LevelList',RF_Ident_Plot.STASD_Gaus_Thresh_height{i});
%                         end
%                     else % Heat_Map_Colour_Choice == 1
%                         if ~isempty(RF_Ident_Plot.STASD_Gaus_Num_RF_pixels{i})
%                             scatter(RF_Ident_Plot.STASD_Gaus_RF_coords{i}(:,2),RF_Ident_Plot.STASD_Gaus_RF_coords{i}(:,1),'go','LineWidth',1.5);
%                             plot(RF_Ident_Plot.STASD_Gaus_Gaussian_mean{i}(1),RF_Ident_Plot.STASD_Gaus_Gaussian_mean{i}(2),'rx','LineWidth',1.5);
%                             fcontour(RF_Ident_Plot.STASD_Gaus_gmPDF{i},[0.5 (stim_columns+0.5) 0.5 (stim_rows+0.5)],'r','LineWidth',1.5,'LevelList',RF_Ident_Plot.STASD_Gaus_Thresh_height{i});
%                         end
%                     end
%                     axis equal; axis tight;
%                     if i>2
%                         xlabel('x');
%                     end
%                     if i==1||i==3
%                         ylabel('y');
%                     end
%                     title(Spectral_Names{i});
%                     colorbar;
%                     set(gca,'FontSize',12);
%                 end
%             end
%             annotation('textbox',[.5 .975 0 0],'String','STA-SD CI -- Gaussian RF','FitBoxToText','on','FontSize',12,'FontWeight','bold','LineStyle','none','HorizontalAlignment','center','VerticalAlignment','middle');
%             set(gcf,'color','w');
%             
%         end
%         
%     end
%     
%        
%     
%     
%     if p.RF_Ident_Meth_vec(2) == 1 % LC method (Local Covariance)
%         
%         if p.RF_Type(1) == 1     % Box
%             
%             fig4 = figure;
%             colormap(fig4,gray(256));
%             imagesc(RF_Ident_Plot.Pixel_covar); hold on;
%             scatter(RF_Ident_Plot.RF_coords_centre(2),RF_Ident_Plot.RF_coords_centre(1),'rx','LineWidth',1.5);
%             scatter(RF_Ident_Plot.RF_coords_noncentre(:,2),RF_Ident_Plot.RF_coords_noncentre(:,1),'go','LineWidth',1.5);
%             axis equal; axis tight;
%             xlabel('x');
%             ylabel('y');
%             title('RF');
%             colorbar;
%             set(gca,'FontSize',12);
%             set(gcf,'color','w');
%             
%         end
%         
%         if p.RF_Type(2) == 1 % All Significant Pixels
%             
%             fig5 = figure;
%             colormap(fig5,gray(256));
%             imagesc(RF_Ident_Plot.Pixel_covar); hold on;
%             scatter(RF_Ident_Plot.RF_coords(:,2),RF_Ident_Plot.RF_coords(:,1),'go','LineWidth',1.5);
%             axis equal; axis tight;
%             xlabel('x');
%             ylabel('y');
%             title('RF');
%             colorbar;
%             set(gca,'FontSize',12);
%             set(gcf,'color','w');
%             
%         end
%         
%         if p.RF_Type(3) == 1 % Gaussian
%             
%             fig6 = figure;
%             colormap(fig6,gray(256));
%             imagesc(RF_Ident_Plot.Pixel_covar); hold on;
%             scatter(RF_Ident_Plot.RF_coords(:,2),RF_Ident_Plot.RF_coords(:,1),'go','LineWidth',1.5);
%             plot(RF_Ident_Plot.Gaussian_mean(1),RF_Ident_Plot.Gaussian_mean(2),'rx','LineWidth',1.5); %
%             fcontour(RF_Ident_Plot.gmPDF,[0.5 (p.stim_columns+0.5) 0.5 (p.stim_rows+0.5)],'r','LineWidth',1.5,'LevelList',RF_Ident_Plot.Thresh_height); % 'LineColor',LightestBlue,'LineWidth',2
%             axis equal; axis tight;
%             xlabel('x');
%             ylabel('y');
%             title('RF');
%             colorbar;
%             set(gca,'FontSize',12);
%             set(gcf,'color','w');
%             
%         end
%         
%     end
%     
%     if p.RF_Ident_Meth_vec(3) == 1 % MI method (Mutual Information)
%         
%         if p.RF_Type(1) == 1     % Box
%             
%             fig7 = figure;
%             colormap(fig7,gray(256));
%             %imagesc(RF_Ident_Plot.Pixel_covar); hold on; --> Update this line!!
%             scatter(RF_Ident_Plot.RF_coords_centre(2),RF_Ident_Plot.RF_coords_centre(1),'rx','LineWidth',1.5);
%             scatter(RF_Ident_Plot.RF_coords_noncentre(:,2),RF_Ident_Plot.RF_coords_noncentre(:,1),'go','LineWidth',1.5);
%             axis equal; axis tight;
%             xlabel('x');
%             ylabel('y');
%             title('RF');
%             colorbar;
%             set(gca,'FontSize',12);
%             set(gcf,'color','w');
%             
%         end
%         
%         if p.RF_Type(2) == 1 % All Significant Pixels
%             
%             fig8 = figure;
%             colormap(fig8,gray(256));
%             %imagesc(RF_Ident_Plot.Pixel_covar); hold on; --> Update this line!!
%             scatter(RF_Ident_Plot.RF_coords(:,2),RF_Ident_Plot.RF_coords(:,1),'go','LineWidth',1.5);
%             axis equal; axis tight;
%             xlabel('x');
%             ylabel('y');
%             title('RF');
%             colorbar;
%             set(gca,'FontSize',12);
%             set(gcf,'color','w');
%             
%         end
%         
%         if p.RF_Type(3) == 1 % Gaussian
%             
%             fig9 = figure;
%             colormap(fig9,gray(256));
%             %imagesc(RF_Ident_Plot.Pixel_covar); hold on; --> Update this line!!
%             scatter(RF_Ident_Plot.RF_coords(:,2),RF_Ident_Plot.RF_coords(:,1),'go','LineWidth',1.5);
%             plot(RF_Ident_Plot.Gaussian_mean(1),RF_Ident_Plot.Gaussian_mean(2),'rx','LineWidth',1.5); %
%             fcontour(RF_Ident_Plot.gmPDF,[0.5 (p.stim_columns+0.5) 0.5 (p.stim_rows+0.5)],'r','LineWidth',1.5,'LevelList',RF_Ident_Plot.Thresh_height); % 'LineColor',LightestBlue,'LineWidth',2
%             axis equal; axis tight;
%             xlabel('x');
%             ylabel('y');
%             title('RF');
%             colorbar;
%             set(gca,'FontSize',12);
%             set(gcf,'color','w');
%             
%         end
%         
%     end
%     
%     
%     if p.RF_Ident_Meth_vec(4) == 1 % SC method;
%         
%         if p.RF_Type(1) == 1     % Box
%             
%             fig10 = figure;
%             if p.Spectral_Dim == 1
%                 colormap(fig10,gray(256));
%                 %                 imagesc(RF_Ident_Plot.SCA_Pixel_covar); hold on;
%                 %                 scatter(RF_Ident_Plot.RF_coords_centre(2),RF_Ident_Plot.RF_coords_centre(1),'rx','LineWidth',1.5);
%                 %                 scatter(RF_Ident_Plot.RF_coords_noncentre(:,2),RF_Ident_Plot.RF_coords_noncentre(:,1),'go','LineWidth',1.5);
%                 axis equal; axis tight;
%                 xlabel('x');
%                 ylabel('y');
%                 title('RF');
%                 colorbar;
%                 set(gca,'FontSize',12);
%             elseif p.Spectral_Dim == 4
%                 if Heat_Map_Colour_Choice == 1
%                     colormap(fig10,gray(256));
%                 end
%                 for i = 1:p.Spectral_Dim
%                     loop_plot = subplot(2,2,i);
%                     imagesc(RF_Ident_Plot.SCA_Pixel_covar(:,:,i)); hold on;
%                     if Heat_Map_Colour_Choice == 2
%                         colormap(loop_plot,colorMap_arr(:,:,i));
%                         if ~isempty(RF_Ident_Plot.SC_Box_Num_RF_pixels{i})
%                             scatter(RF_Ident_Plot.SC_Box_RF_coords_centre{i}(2),RF_Ident_Plot.SC_Box_RF_coords_centre{i}(1),'wx','LineWidth',1.5);
%                             scatter(RF_Ident_Plot.SC_Box_RF_coords_noncentre{i}(:,2),RF_Ident_Plot.SC_Box_RF_coords_noncentre{i}(:,1),'wo','LineWidth',1.5);
%                         end
%                     else % Heat_Map_Colour_Choice == 1
%                         if ~isempty(RF_Ident_Plot.SC_Box_Num_RF_pixels{i})
%                             scatter(RF_Ident_Plot.SC_Box_RF_coords_centre{i}(2),RF_Ident_Plot.SC_Box_RF_coords_centre{i}(1),'rx','LineWidth',1.5);
%                             scatter(RF_Ident_Plot.SC_Box_RF_coords_noncentre{i}(:,2),RF_Ident_Plot.SC_Box_RF_coords_noncentre{i}(:,1),'go','LineWidth',1.5);
%                         end
%                     end
%                     axis equal; axis tight;
%                     if i>2
%                         xlabel('x');
%                     end
%                     if i==1||i==3
%                         ylabel('y');
%                     end
%                     title(Spectral_Names{i});
%                     colorbar;
%                     set(gca,'FontSize',12);
%                 end
%             end
%             annotation('textbox',[.5 .975 0 0],'String','Self Covar. -- Box RF','FitBoxToText','on','FontSize',12,'FontWeight','bold','LineStyle','none','HorizontalAlignment','center','VerticalAlignment','middle');
%             set(gcf,'color','w');
%             
%         end
%         
%         if p.RF_Type(2) == 1 % All Significant Pixels
%             
%             fig11 = figure;
%             if p.Spectral_Dim == 1
%                 colormap(fig11,gray(256));
%                 %                 imagesc(RF_Ident_Plot.SCA_Pixel_covar); hold on;
%                 %                 scatter(RF_Ident_Plot.RF_coords(:,2),RF_Ident_Plot.RF_coords(:,1),'go','LineWidth',1.5);
%                 axis equal; axis tight;
%                 xlabel('x');
%                 ylabel('y');
%                 title('RF');
%                 colorbar;
%                 set(gca,'FontSize',12);
%             elseif p.Spectral_Dim == 4
%                 if Heat_Map_Colour_Choice == 1
%                     colormap(fig11,gray(256));
%                 end
%                 for i = 1:p.Spectral_Dim
%                     loop_plot = subplot(2,2,i);
%                     imagesc(RF_Ident_Plot.SCA_Pixel_covar(:,:,i)); hold on;
%                     if Heat_Map_Colour_Choice == 2
%                         colormap(loop_plot,colorMap_arr(:,:,i));
%                         if ~isempty(RF_Ident_Plot.SC_ASP_Num_RF_pixels{i})
%                             scatter(RF_Ident_Plot.SC_ASP_RF_coords{i}(:,2),RF_Ident_Plot.SC_ASP_RF_coords{i}(:,1),'wo','LineWidth',1.5);
%                         end
%                     else % Heat_Map_Colour_Choice == 1
%                         if ~isempty(RF_Ident_Plot.SC_ASP_Num_RF_pixels{i})
%                             scatter(RF_Ident_Plot.SC_ASP_RF_coords{i}(:,2),RF_Ident_Plot.SC_ASP_RF_coords{i}(:,1),'go','LineWidth',1.5);
%                         end
%                     end
%                     axis equal; axis tight;
%                     if i>2
%                         xlabel('x');
%                     end
%                     if i==1||i==3
%                         ylabel('y');
%                     end
%                     title(Spectral_Names{i});
%                     colorbar;
%                     set(gca,'FontSize',12);
%                 end
%             end
%             annotation('textbox',[.5 .975 0 0],'String','Self Covar. -- All Sig. Pix. RF','FitBoxToText','on','FontSize',12,'FontWeight','bold','LineStyle','none','HorizontalAlignment','center','VerticalAlignment','middle');
%             set(gcf,'color','w');
%             
%         end
%         
%         if p.RF_Type(3) == 1 % Gaussian
%             
%             fig12 = figure;
%             if p.Spectral_Dim == 1
%                 colormap(fig12,gray(256));
% %                 imagesc(RF_Ident_Plot.SCA_Pixel_covar,[min(RF_Ident_Plot.SCA_Pixel_covar,[],'all') max(RF_Ident_Plot.SCA_Pixel_covar,[],'all')]); hold on;
% %                 scatter(RF_Ident_Plot.SC_Gaus_RF_coords(:,2),RF_Ident_Plot.SC_Gaus_RF_coords(:,1),'go','LineWidth',1.5);
% %                 plot(RF_Ident_Plot.SC_Gaus_Gaussian_mean(1),RF_Ident_Plot.SC_Gaus_Gaussian_mean(2),'rx','LineWidth',1.5); %
% %                 fcontour(RF_Ident_Plot.SC_Gaus_gmPDF,[0.5 (p.stim_columns+0.5) 0.5 (p.stim_rows+0.5)],'r','LineWidth',1.5,'LevelList',RF_Ident_Plot.Thresh_height); % 'LineColor',LightestBlue,'LineWidth',2
%                 axis equal; axis tight;
%                 xlabel('x');
%                 ylabel('y');
%                 title('RF');
%                 colorbar;
%                 set(gca,'FontSize',12);
%             elseif p.Spectral_Dim == 4
%                 if Heat_Map_Colour_Choice == 1
%                     colormap(fig12,gray(256));
%                 end
%                 for i = 1:p.Spectral_Dim
%                     loop_plot = subplot(2,2,i);
%                     imagesc(RF_Ident_Plot.SCA_Pixel_covar(:,:,i),[min(RF_Ident_Plot.SCA_Pixel_covar(:,:,i),[],'all') max(RF_Ident_Plot.SCA_Pixel_covar(:,:,i),[],'all')]); hold on;
%                     if Heat_Map_Colour_Choice == 2
%                         colormap(loop_plot,colorMap_arr(:,:,i));
%                         if ~isempty(RF_Ident_Plot.SC_Gaus_Num_RF_pixels{i})
%                             scatter(RF_Ident_Plot.SC_Gaus_RF_coords{i}(:,2),RF_Ident_Plot.SC_Gaus_RF_coords{i}(:,1),'wo','LineWidth',1.5);
%                             plot(RF_Ident_Plot.SC_Gaus_Gaussian_mean{i}(1),RF_Ident_Plot.SC_Gaus_Gaussian_mean{i}(2),'wx','LineWidth',1.5);
%                             fcontour(RF_Ident_Plot.SC_Gaus_gmPDF{i},[0.5 (p.stim_columns+0.5) 0.5 (p.stim_rows+0.5)],'w','LineWidth',1.5,'LevelList',RF_Ident_Plot.SC_Gaus_Thresh_height{i});
%                         end
%                     else % Heat_Map_Colour_Choice == 1
%                         if ~isempty(RF_Ident_Plot.SC_Gaus_Num_RF_pixels{i})
%                             scatter(RF_Ident_Plot.SC_Gaus_RF_coords{i}(:,2),RF_Ident_Plot.SC_Gaus_RF_coords{i}(:,1),'go','LineWidth',1.5);
%                             plot(RF_Ident_Plot.SC_Gaus_Gaussian_mean{i}(1),RF_Ident_Plot.SC_Gaus_Gaussian_mean{i}(2),'rx','LineWidth',1.5);
%                             fcontour(RF_Ident_Plot.SC_Gaus_gmPDF{i},[0.5 (p.stim_columns+0.5) 0.5 (p.stim_rows+0.5)],'r','LineWidth',1.5,'LevelList',RF_Ident_Plot.SC_Gaus_Thresh_height{i});
%                         end
%                     end
%                     axis equal; axis tight;
%                     if i>2
%                         xlabel('x');
%                     end
%                     if i==1||i==3
%                         ylabel('y');
%                     end
%                     title(Spectral_Names{i});
%                     colorbar;
%                     set(gca,'FontSize',12);
%                 end
%             end
%             annotation('textbox',[.5 .975 0 0],'String','Self Covar. -- Gaussian RF','FitBoxToText','on','FontSize',12,'FontWeight','bold','LineStyle','none','HorizontalAlignment','center','VerticalAlignment','middle');
%             set(gcf,'color','w');
%             
%         end
%         
%     end
%     
%     
%     
%     
%     
%     
%     
%     
%     %%% RF Significance Plots
%     
%     if p.RF_Ident_Meth_vec(1) == 1 % STA-SD method;
%         
%         figure;
%         for i = 1:p.Spectral_Dim
%             subplot(2,2,i);
%             hist_loop = histogram(RF_Ident_Plot.STA_SD(:,:,i)); hold on;
%             if Hist_Colour_Choice == 2
%                 set(hist_loop,'FaceColor',colorMap_arr(end,:,i));
%             end
%             vline_loop_1 = vline(mean(RF_Ident_Plot.STA_SD(:,:,i),'all'),'r'); % mean
%             vline_loop_2 = vline(RF_Ident_Plot.STA_SD_Sig_Thresh(i),'g'); % significance thrshold
%             vline_loop_3 = vline(hist_loop.BinEdges(end),'b:');   % RH edge of rightmost bin
%             if i>2
%                 xlabel('STA-SD.');
%             end
%             if i==1||i==3
%                 ylabel('num. pixels');
%             end
%             title(Spectral_Names{i});
%             set(gca,'FontSize',12);
%             set(vline_loop_1,'LineWidth',1.5);
%             set(vline_loop_2,'LineWidth',1.5);
%             set(vline_loop_3,'LineWidth',1.5);
%             xlim([hist_loop.BinEdges(1) max(hist_loop.BinEdges(end),RF_Ident_Plot.STA_SD_Sig_Thresh(i))]);
%             ylim([0 max(hist_loop.Values)]);
%         end
%         annotation('textbox',[.5 .975 0 0],'String','STA-SD','FitBoxToText','on','FontSize',12,'FontWeight','bold','LineStyle','none','HorizontalAlignment','center','VerticalAlignment','middle'); % dim: [x y w h] [.45 .7 .3 .3]
%         set(gcf,'color','w');
%         
%     end
%     
%     if p.RF_Ident_Meth_vec(2) == 1 % LC method (Local Covariance)
%         
%         figure;
%         for i = 1:p.Spectral_Dim
%             subplot(2,2,i);
%             hist_loop = histogram(RF_Ident_Plot.LCA_Stixel_covar(:,:,i)); hold on;
%             if Hist_Colour_Choice == 2
%                 set(hist_loop,'FaceColor',colorMap_arr(end,:,i));
%             end
%             vline_loop_1 = vline(mean(RF_Ident_Plot.LCA_Stixel_covar(:,:,i),'all'),'r'); % mean
%             vline_loop_2 = vline(RF_Ident_Plot.LC_Sig_Thresh(i),'g');   % significance thrshold
%             vline_loop_3 = vline(hist_loop.BinEdges(end),'b:'); % RH edge of rightmost bin
%             if i>2
%                 xlabel('STA-SD.');
%             end
%             if i==1||i==3
%                 ylabel('num. pixels');
%             end
%             title(Spectral_Names{i});
%             set(gca,'FontSize',12);
%             set(vline_loop_1,'LineWidth',1.5);
%             set(vline_loop_2,'LineWidth',1.5);
%             set(vline_loop_3,'LineWidth',1.5);
%             xlim([hist_loop.BinEdges(1) max(hist_loop.BinEdges(end),RF_Ident_Plot.LC_Sig_Thresh(i))]);
%             ylim([0 max(hist_loop.Values)]);
%         end
%         annotation('textbox',[.5 .975 0 0],'String','Local Covariance','FitBoxToText','on','FontSize',12,'FontWeight','bold','LineStyle','none','HorizontalAlignment','center','VerticalAlignment','middle');
%         set(gcf,'color','w');
%         
%     end
%     
%     if p.RF_Ident_Meth_vec(3) == 1 % MI method (Mutual Information)
%         
%         figure;
%         for i = 1:p.Spectral_Dim
%             subplot(2,2,i);
%             %hist_loop = histogram(); hold on;
%             if Hist_Colour_Choice == 2
%                 set(hist_loop,'FaceColor',colorMap_arr(end,:,i));
%             end
%             %vline_loop_1 = vline(mean(,'all'),'r'); % mean
%             vline_loop_2 = vline(RF_Ident_Plot.MI_Sig_Thresh(i),'g');   % significance thrshold
%             vline_loop_3 = vline(hist_loop.BinEdges(end),'b:'); % RH edge of rightmost bin
%             if i>2
%                 xlabel('STA-SD.');
%             end
%             if i==1||i==3
%                 ylabel('num. pixels');
%             end
%             title(Spectral_Names{i});
%             set(gca,'FontSize',12);
%             set(vline_loop_1,'LineWidth',1.5);
%             set(vline_loop_2,'LineWidth',1.5);
%             set(vline_loop_3,'LineWidth',1.5);
%             xlim([hist_loop.BinEdges(1) max(hist_loop.BinEdges(end),RF_Ident_Plot.MI_Sig_Thresh(i))]);
%             ylim([0 max(hist_loop.Values)]);
%         end
%         annotation('textbox',[.5 .975 0 0],'String','Mutual Information','FitBoxToText','on','FontSize',12,'FontWeight','bold','LineStyle','none','HorizontalAlignment','center','VerticalAlignment','middle');
%         set(gcf,'color','w');
%         
%     end
%     
%     if p.RF_Ident_Meth_vec(4) == 1 % SC method;
%         
%         figure;
%         for i = 1:p.Spectral_Dim
%             subplot(2,2,i);
%             hist_loop = histogram(RF_Ident_Plot.SCA_Pixel_covar(:,:,i)); hold on;
%             if Hist_Colour_Choice == 2
%                 set(hist_loop,'FaceColor',colorMap_arr(end,:,i));
%             end
%             vline_loop_1 = vline(mean(RF_Ident_Plot.SCA_Pixel_covar(:,:,i),'all'),'r'); % mean
%             vline_loop_2 = vline(RF_Ident_Plot.SC_Sig_Thresh(i),'g'); % significance thrshold
%             vline_loop_3 = vline(hist_loop.BinEdges(1),'b:'); % LH edge of leftmost bin
%             if i>2
%                 xlabel('self-covar.');
%             end
%             if i==1||i==3
%                 ylabel('num. pixels');
%             end
%             title(Spectral_Names{i});
%             set(gca,'FontSize',12);
%             set(vline_loop_1,'LineWidth',1.5);
%             set(vline_loop_2,'LineWidth',1.5);
%             set(vline_loop_3,'LineWidth',1.5);
%             xlim([min(hist_loop.BinEdges(1),RF_Ident_Plot.SC_Sig_Thresh(i)) hist_loop.BinEdges(end)]);
%             ylim([0 max(hist_loop.Values)]);
%         end
%         annotation('textbox',[.5 .975 0 0],'String','Self Covariance','FitBoxToText','on','FontSize',12,'FontWeight','bold','LineStyle','none','HorizontalAlignment','center','VerticalAlignment','middle');
%         set(gcf,'color','w');
%         
%     end
%     
%     
%     
%     
%      %%% Pure Heatmap Plots
%      
%      if p.RF_Ident_Meth_vec(1) == 1 % STA-SD method;
%          
%          fig13 = figure;
%          if p.Spectral_Dim == 1
%              colormap(fig13,gray(256));
%              %                 imagesc(RF_Ident_Plot.STA_SD); hold on;
%              axis equal; axis tight;
%              xlabel('x');
%              ylabel('y');
%              title('RF');
%              colorbar;
%              set(gca,'FontSize',12);
%          elseif p.Spectral_Dim == 4
%              if Heat_Map_Colour_Choice == 1
%                  colormap(fig13,gray(256));
%              end
%              for i = 1:p.Spectral_Dim
%                  loop_plot = subplot(2,2,i);
%                  imagesc(RF_Ident_Plot.STA_SD(:,:,i)); hold on;
%                  if Heat_Map_Colour_Choice == 2
%                      colormap(loop_plot,colorMap_arr(:,:,i));
%                  end
%                  axis equal; axis tight;
%                  if i>2
%                      xlabel('x');
%                  end
%                  if i==1||i==3
%                      ylabel('y');
%                  end
%                  title(Spectral_Names{i});
%                  colorbar;
%                  set(gca,'FontSize',12);
%              end
%          end
%          annotation('textbox',[.5 .975 0 0],'String','STA-SD','FitBoxToText','on','FontSize',12,'FontWeight','bold','LineStyle','none','HorizontalAlignment','center','VerticalAlignment','middle');
%          set(gcf,'color','w');
%          
%      end
%      
%      if p.RF_Ident_Meth_vec(2) == 1 % LC method;
%          
%          fig14 = figure;
%          if p.Spectral_Dim == 1
%              colormap(fig14,gray(256));
%              %                 imagesc(RF_Ident_Plot.Pixel_covar); hold on;
%              axis equal; axis tight;
%              xlabel('x');
%              ylabel('y');
%              title('RF');
%              colorbar;
%              set(gca,'FontSize',12);
%          elseif p.Spectral_Dim == 4
%              if Heat_Map_Colour_Choice == 1
%                  colormap(fig14,gray(256));
%              end
%              for i = 1:p.Spectral_Dim
%                  loop_plot = subplot(2,2,i);
%                  imagesc(RF_Ident_Plot.Pixel_covar(:,:,i)); hold on;
%                  if Heat_Map_Colour_Choice == 2
%                      colormap(loop_plot,colorMap_arr(:,:,i));
%                  end
%                  axis equal; axis tight;
%                  if i>2
%                      xlabel('x');
%                  end
%                  if i==1||i==3
%                      ylabel('y');
%                  end
%                  title(Spectral_Names{i});
%                  colorbar;
%                  set(gca,'FontSize',12);
%              end
%          end
%          annotation('textbox',[.5 .975 0 0],'String','Local Covariance','FitBoxToText','on','FontSize',12,'FontWeight','bold','LineStyle','none','HorizontalAlignment','center','VerticalAlignment','middle');
%          set(gcf,'color','w');
%          
%      end
%      
%      if p.RF_Ident_Meth_vec(3) == 1 % MI method;
%          
%          fig15 = figure;
%          if p.Spectral_Dim == 1
%              colormap(fig15,gray(256));
%              %                 imagesc(RF_Ident_Plot.Pixel_covar); hold on;
%              axis equal; axis tight;
%              xlabel('x');
%              ylabel('y');
%              title('RF');
%              colorbar;
%              set(gca,'FontSize',12);
%          elseif p.Spectral_Dim == 4
%              if Heat_Map_Colour_Choice == 1
%                  colormap(fig15,gray(256));
%              end
%              for i = 1:p.Spectral_Dim
%                  loop_plot = subplot(2,2,i);
%                  %imagesc(RF_Ident_Plot.Pixel_covar(:,:,i)); hold on; Update!!!
%                  if Heat_Map_Colour_Choice == 2
%                      colormap(loop_plot,colorMap_arr(:,:,i));
%                  end
%                  axis equal; axis tight;
%                  if i>2
%                      xlabel('x');
%                  end
%                  if i==1||i==3
%                      ylabel('y');
%                  end
%                  title(Spectral_Names{i});
%                  colorbar;
%                  set(gca,'FontSize',12);
%              end
%          end
%          annotation('textbox',[.5 .975 0 0],'String','Mutual Information','FitBoxToText','on','FontSize',12,'FontWeight','bold','LineStyle','none','HorizontalAlignment','center','VerticalAlignment','middle');
%          set(gcf,'color','w');
%          
%      end
%      
%      if p.RF_Ident_Meth_vec(4) == 1 % SC method;
%          
%          fig16 = figure;
%          if p.Spectral_Dim == 1
%              colormap(fig16,gray(256));
%              %                 imagesc(RF_Ident_Plot.SCA_Pixel_covar); hold on;
%              axis equal; axis tight;
%              xlabel('x');
%              ylabel('y');
%              title('RF');
%              colorbar;
%              set(gca,'FontSize',12);
%          elseif p.Spectral_Dim == 4
%              if Heat_Map_Colour_Choice == 1
%                  colormap(fig16,gray(256));
%              end
%              for i = 1:p.Spectral_Dim
%                  loop_plot = subplot(2,2,i);
%                  imagesc(RF_Ident_Plot.SCA_Pixel_covar(:,:,i)); hold on;
%                  if Heat_Map_Colour_Choice == 2
%                      colormap(loop_plot,colorMap_arr(:,:,i));
%                  end
%                  axis equal; axis tight;
%                  if i>2
%                      xlabel('x');
%                  end
%                  if i==1||i==3
%                      ylabel('y');
%                  end
%                  title(Spectral_Names{i});
%                  colorbar;
%                  set(gca,'FontSize',12);
%              end
%          end
%          annotation('textbox',[.5 .975 0 0],'String','Self Covariance','FitBoxToText','on','FontSize',12,'FontWeight','bold','LineStyle','none','HorizontalAlignment','center','VerticalAlignment','middle');
%          set(gcf,'color','w');
%          
%      end
%     
%     
% end
% 
% %% Save data
% 
% %save data_MChick_1_Cell_7_v_1;
% 
% 
% 
% 
% 
% 
% % figure;
% % colormap(gray(256));
% % subplot(1,5,1);
% % imagesc(RF_Ident_Plot.Stixel_covar(:,:,1)); axis equal; axis tight;
% % subplot(1,5,2);
% % imagesc(RF_Ident_Plot.Stixel_covar(:,:,2)); axis equal; axis tight;
% % subplot(1,5,3);
% % imagesc(RF_Ident_Plot.Stixel_covar(:,:,3)); axis equal; axis tight;
% % subplot(1,5,4);
% % imagesc(RF_Ident_Plot.Stixel_covar(:,:,4)); axis equal; axis tight;
% % subplot(1,5,5);
% % imagesc(RF_Ident_Plot.Stixel_covar(:,:,5)); axis equal; axis tight;
% % set(gcf,'color','w');     
% %      
% %      
% %      
% % fig2 = figure;
% % colormap(fig2,gray(256));
% % xslice = [];   
% % yslice = [];
% % zslice = [1,2,3,4,5];
% % slice(RF_Ident_Plot.Stixel_covar,xslice,yslice,zslice);
