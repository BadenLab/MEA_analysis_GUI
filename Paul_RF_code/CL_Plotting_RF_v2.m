%% Plotting RF Command Line -- ver. 2
% Use this for RF Ident CL 4 onwards

% clear all;
% clc;
% load data_RF_Ident_CL3_MChick_1_Cell_1to5000_v_1.mat;
function out = CL_Plotting_RF_v2 (savepath, add_info)
%% Plotting options


%% Global variables


% Choose whether to use B&W or coloured heatmaps and histograms
% 1 = B&W;
% 2 = coloured.
if add_info.settings.kernel_new.plot_all.heat_colour
    Heat_Map_Colour_Choice = 2;
elseif add_info.settings.kernel_new.plot_all.heat_gray
    Heat_Map_Colour_Choice = 1;
end
if add_info.settings.kernel_new.plot_all.hist_heat_colour
    Hist_Colour_Choice = 2;
elseif add_info.settings.kernel_new.plot_all.hist_heat_gray
    Hist_Colour_Choice = 1;
end

% Choose whether to align histogram bins
% 1: Yes
% 2: No
Hist_Bin_Align_Choice = add_info.settings.kernel_new.plot_all.alignhisto;

% Allow for 2 colour cases: monochromatic and tetrachromatic (RGBUV)
Spectral_Names = {'R','G','B','UV'}; %mars: Has to be not hard coded

colorMap_arr        = NaN(256,3,4);
colorMap_arr(:,:,1) = [linspace(0,1,256)', zeros(256,2)];                     % Red
colorMap_arr(:,:,2) = [zeros(256,1),linspace(0,1,256)', zeros(256,1)];        % Green
colorMap_arr(:,:,3) = [zeros(256,2), linspace(0,1,256)'];                     % Blue
colorMap_arr(:,:,4) = [linspace(0,1,256)', zeros(256,1), linspace(0,1,256)']; % UV
%colorMapBlack      = [0 0 0];

%% Load the dataset information

%Check how many folder shall be plotted
for ff = 1:length(add_info.settings.kernel_new.folderplot)


if ~isfield(add_info.settings.kernel_new,'folderplot')
    error('No analysis folder selected, do that first an run this script afterwards')
end
plot_folder = add_info.settings.kernel_new.folderplot{1,ff};
L = load(findfile_app(add_info.stim_idx,savepath,"RF_overview.mat",'subfolder',plot_folder));
RF_overview = L.RF_overview;
clear L

%Check for how many cells a RF has been found
try
    true_num_cells = nnz([RF_overview.STASD_ASP_RF]);
    true_cells_idx = [RF_overview.STASD_ASP_RF] == 1;
catch
    true_num_cells = nnz([RF_overview.SC_ASP_RF]);
    true_cells_idx = [RF_overview.SC_ASP_RF] == 1;
end

path_to_cells = [RF_overview.file];
path_to_cells = path_to_cells(1,true_cells_idx);
true_cells_cell_idx = [RF_overview.cell_idx];
true_cells_cell_idx = true_cells_cell_idx(1,true_cells_idx);


%Recreate RF.Ident cell
for ii = 1:true_num_cells
    RF_file = strcat(path_to_cells(1,ii),'.bin');
    %Now we load the file and transfer it back from binary to matlab
    %formate
    fileID = fopen(RF_file,'r');
    RF_Ident_temp = fread(fileID);
    fclose(fileID);
    RF_Ident_temp = hlp_deserialize(RF_Ident_temp);
    if ii == 1
        RF_Ident = RF_Ident_temp;
    else
        RF_Ident = [RF_Ident,RF_Ident_temp]; %Build a new structure from the ones loaded
    end
    
end
%Recreate 'p' variable
p = RF_overview.p;
p.Spectral_Dim = 4; %mars: This should be stored in the variable and not be hard coded
%% Extract Data to Plot


if contains(add_info.settings.kernel_new.folderplot{1},'SS')
    
%     if p.RF_Type(1) == 1 % Box %mars: not sure why this is here?
%         
%     end
    
    if add_info.settings.kernel_new.plot_Allpixel   % All Significant Pixels
        
        % Num Spectral and Full RF Pixels
        Num_STASD_ASP_RF_Pixels_mat     = NaN(true_num_cells,p.Spectral_Dim);
        Num_STASD_ASP_FullRF_Pixels_vec = NaN(true_num_cells,1);
        for i = 1:true_num_cells
            for j = 1:p.Spectral_Dim
                if ~isempty(RF_Ident(i).RF_results.STASD_ASP_Num_RF_pixels{j})
                    Num_STASD_ASP_RF_Pixels_mat(i,j) = RF_Ident(i).RF_results.STASD_ASP_Num_RF_pixels{j};
                end
            end
            if ~isempty(RF_Ident(i).RF_results.STASD_ASP_FullRF_Num_pixels)
                Num_STASD_ASP_FullRF_Pixels_vec(i) = RF_Ident(i).RF_results.STASD_ASP_FullRF_Num_pixels;
            end
        end
        Num_STASD_ASP_RF_Pixels_mat(sum(isnan(Num_STASD_ASP_RF_Pixels_mat),2)==p.Spectral_Dim,:)             = []; % Vector with NaN rows removed
        Num_STASD_ASP_FullRF_Pixels_vec(isnan(Num_STASD_ASP_FullRF_Pixels_vec))                              = []; % Vector with NaN rows removed
        Num_STASD_ASP_FullRF_Pixels_vec(Num_STASD_ASP_FullRF_Pixels_vec==0)                                  = []; % Vector with zero rows removed
        
    end
    
    if add_info.settings.kernel_new.plot_Gauss   % Gaussian
        
        % Num Pixels Spectral and Full RFs
        % Abs Ellipse Area Spectral and Full RFs
        % Major Axis Angle Spectral and Full RFs
        % Axis Eval Ratio Spectral and Full RFs
        % Mean Coords Spectral and Full RFs
        STASD_Gaus_Num_RF_pixels_mat           = NaN(true_num_cells,p.Spectral_Dim);
        STASD_Gaus_FullRF_Num_pixels_vec       = NaN(true_num_cells,1);
        STASD_Gaus_Abs_Ellipse_Area_mat        = NaN(true_num_cells,p.Spectral_Dim);
        STASD_Gaus_FullRF_Abs_Ellipse_Area_vec = NaN(true_num_cells,1);
        STASD_Gaus_RF_Angle_Major_Axis_mat     = NaN(true_num_cells,p.Spectral_Dim);
        STASD_Gaus_FullRF_Angle_Major_Axis_vec = NaN(true_num_cells,1);
        STASD_Gaus_Axis_eval_Ratio_mat         = NaN(true_num_cells,p.Spectral_Dim);
        STASD_Gaus_FullRF_Axis_eval_Ratio_vec  = NaN(true_num_cells,1);
        STASD_Gaus_Gaussian_mean_arr           = NaN(true_num_cells,p.Spectral_Dim,2);
        STASD_Gaus_FullRF_Gaussian_mean_mat    = NaN(true_num_cells,2);

        for i = 1:true_num_cells
            for j = 1:p.Spectral_Dim
                if ~isempty(RF_Ident(i).RF_results.STASD_Gaus_Num_RF_pixels{j})
                    STASD_Gaus_Num_RF_pixels_mat(i,j)       = RF_Ident(i).RF_results.STASD_Gaus_Num_RF_pixels{j};
                    STASD_Gaus_Abs_Ellipse_Area_mat(i,j)    = RF_Ident(i).RF_results.STASD_Gaus_Abs_Ellipse_Area(j);
                    STASD_Gaus_RF_Angle_Major_Axis_mat(i,j) = RF_Ident(i).RF_results.STASD_Gaus_Angle_Major_Axis(j);
                    STASD_Gaus_Axis_eval_Ratio_mat(i,j)     = RF_Ident(i).RF_results.STASD_Gaus_Axis_eval_Ratio(j);
                    STASD_Gaus_Gaussian_mean_arr(i,j,:)     = [RF_Ident(i).RF_results.STASD_Gaus_Gaussian_mean{j}(1),...
                        RF_Ident(i).RF_results.STASD_Gaus_Gaussian_mean{j}(2)];
                end
            end
            if RF_Ident(i).RF_results.STASD_Gaus_FullRF_Num_pixels ~=0
                STASD_Gaus_FullRF_Num_pixels_vec(i)       = RF_Ident(i).RF_results.STASD_Gaus_FullRF_Num_pixels;
                STASD_Gaus_FullRF_Abs_Ellipse_Area_vec(i) = RF_Ident(i).RF_results.STASD_Gaus_FullRF_Abs_Ellipse_Area;
                STASD_Gaus_FullRF_Angle_Major_Axis_vec(i) = RF_Ident(i).RF_results.STASD_Gaus_FullRF_Angle_Major_Axis;
                STASD_Gaus_FullRF_Axis_eval_Ratio_vec(i)  = RF_Ident(i).RF_results.STASD_Gaus_FullRF_Axis_eval_Ratio;
                STASD_Gaus_FullRF_Gaussian_mean_mat(i,:)  = [RF_Ident(i).RF_results.STASD_Gaus_FullRF_Gaussian_mean(1),...
                    RF_Ident(i).RF_results.STASD_Gaus_FullRF_Gaussian_mean(2)];
            end
        end
        STASD_Gaus_Num_RF_pixels_mat_2     = STASD_Gaus_Num_RF_pixels_mat;                                           % Full vector
        STASD_Gaus_FullRF_Num_pixels_vec_2 = STASD_Gaus_FullRF_Num_pixels_vec;                                       % Full vector
        STASD_Gaus_Num_RF_pixels_mat(sum(isnan(STASD_Gaus_Num_RF_pixels_mat),2)==p.Spectral_Dim,:)             = []; % Vector with NaN rows removed
        STASD_Gaus_FullRF_Num_pixels_vec(isnan(STASD_Gaus_FullRF_Num_pixels_vec))                              = []; % Vector with NaN rows removed
        STASD_Gaus_Abs_Ellipse_Area_mat(sum(isnan(STASD_Gaus_Abs_Ellipse_Area_mat),2)==p.Spectral_Dim,:)       = [];
        STASD_Gaus_FullRF_Abs_Ellipse_Area_vec(isnan(STASD_Gaus_FullRF_Abs_Ellipse_Area_vec))                  = [];
        STASD_Gaus_RF_Angle_Major_Axis_mat(sum(isnan(STASD_Gaus_RF_Angle_Major_Axis_mat),2)==p.Spectral_Dim,:) = [];
        STASD_Gaus_FullRF_Angle_Major_Axis_vec(isnan(STASD_Gaus_FullRF_Angle_Major_Axis_vec))                  = [];   
        STASD_Gaus_Axis_eval_Ratio_mat(sum(isnan(STASD_Gaus_Axis_eval_Ratio_mat),2)==p.Spectral_Dim,:)         = [];
        STASD_Gaus_FullRF_Axis_eval_Ratio_vec(isnan(STASD_Gaus_FullRF_Axis_eval_Ratio_vec))                    = [];   
        STASD_Gaus_Gaussian_mean_arr(sum(isnan(STASD_Gaus_Gaussian_mean_arr(:,:,1)),2)==p.Spectral_Dim,:,:)    = [];
        STASD_Gaus_FullRF_Gaussian_mean_mat(isnan(STASD_Gaus_FullRF_Gaussian_mean_mat(:,1)),:)                 = [];
  
    end
    
end

% if     p.RF_Ident_Meth_vec(2) == 1 % LC
%     
%     if p.RF_Type(1) == 1 % Box
%         
%     end
%     
%     if p.RF_Type(2) == 1     % All Significant Pixels
%         
%     end
%     
%     if p.RF_Type(3) == 1     % Gaussian
%         
%     end
%     
% end
% 
% if     p.RF_Ident_Meth_vec(3) == 1 % MI
%     
%     if p.RF_Type(1) == 1 % Box
%         
%     end
%     
%     if p.RF_Type(2) == 1     % All Significant Pixels
%         
%     end
%     
%     if p.RF_Type(3) == 1     % Gaussian
%         
%     end
%     
% end

if contains(add_info.settings.kernel_new.folderplot{1},'SC')
    
%     if p.RF_Type(1) == 1 % Box
%         
%     end
    
    if add_info.settings.kernel_new.plot_Allpixel % All Significant Pixels
        
        % Num Spectral and Full RF Pixels
        Num_SC_ASP_RF_Pixels_mat     = NaN(true_num_cells,p.Spectral_Dim);
        Num_SC_ASP_FullRF_Pixels_vec = NaN(true_num_cells,1);
        for i = 1:true_num_cells
            for j = 1:p.Spectral_Dim
                if ~isempty(RF_Ident(i).RF_results.SC_ASP_Num_RF_pixels{j})
                    Num_SC_ASP_RF_Pixels_mat(i,j) = RF_Ident(i).RF_results.SC_ASP_Num_RF_pixels{j};
                end
            end
            if ~isempty(RF_Ident(i).RF_results.SC_ASP_FullRF_Num_pixels)
                Num_SC_ASP_FullRF_Pixels_vec(i) = RF_Ident(i).RF_results.SC_ASP_FullRF_Num_pixels;
            end
        end
        Num_SC_ASP_RF_Pixels_mat(sum(isnan(Num_SC_ASP_RF_Pixels_mat),2)==p.Spectral_Dim,:)             = []; % Vector with NaN rows removed
        Num_SC_ASP_FullRF_Pixels_vec(isnan(Num_SC_ASP_FullRF_Pixels_vec))                              = []; % Vector with NaN rows removed
        Num_SC_ASP_FullRF_Pixels_vec(Num_SC_ASP_FullRF_Pixels_vec==0)                                  = []; % Vector with zero rows removed
        
    end
    
    if p.RF_Type(3) == 1     % Gaussian
        
        % Num Pixels Spectral and Full RFs
        % Abs Ellipse Area Spectral and Full RFs
        % Major Axis Angle Spectral and Full RFs
        % Axis Eval Ratio Spectral and Full RFs
        % Mean Coords Spectral and Full RFs
        SC_Gaus_Num_RF_pixels_mat           = NaN(true_num_cells,p.Spectral_Dim);
        SC_Gaus_FullRF_Num_pixels_vec       = NaN(true_num_cells,1);
        SC_Gaus_Abs_Ellipse_Area_mat        = NaN(true_num_cells,p.Spectral_Dim);
        SC_Gaus_FullRF_Abs_Ellipse_Area_vec = NaN(true_num_cells,1);
        SC_Gaus_RF_Angle_Major_Axis_mat     = NaN(true_num_cells,p.Spectral_Dim);
        SC_Gaus_FullRF_Angle_Major_Axis_vec = NaN(true_num_cells,1);
        SC_Gaus_Axis_eval_Ratio_mat         = NaN(true_num_cells,p.Spectral_Dim);
        SC_Gaus_FullRF_Axis_eval_Ratio_vec  = NaN(true_num_cells,1);
        SC_Gaus_Gaussian_mean_arr           = NaN(true_num_cells,p.Spectral_Dim,2);
        SC_Gaus_FullRF_Gaussian_mean_mat    = NaN(true_num_cells,2);

        for i = 1:true_num_cells
            for j = 1:p.Spectral_Dim
                if ~isempty(RF_Ident(i).RF_results.SC_Gaus_Num_RF_pixels{j})
                    SC_Gaus_Num_RF_pixels_mat(i,j)       = RF_Ident(i).RF_results.SC_Gaus_Num_RF_pixels{j};
                    SC_Gaus_Abs_Ellipse_Area_mat(i,j)    = RF_Ident(i).RF_results.SC_Gaus_Abs_Ellipse_Area(j);
                    SC_Gaus_RF_Angle_Major_Axis_mat(i,j) = RF_Ident(i).RF_results.SC_Gaus_Angle_Major_Axis(j);
                    SC_Gaus_Axis_eval_Ratio_mat(i,j)     = RF_Ident(i).RF_results.SC_Gaus_Axis_eval_Ratio(j);
                    SC_Gaus_Gaussian_mean_arr(i,j,:)     = [RF_Ident(i).RF_results.SC_Gaus_Gaussian_mean{j}(1),RF_Ident(i).RF_results.SC_Gaus_Gaussian_mean{j}(2)];
                end
            end
            if RF_Ident(i).RF_results.SC_Gaus_FullRF_Num_pixels ~=0
                SC_Gaus_FullRF_Num_pixels_vec(i)       = RF_Ident(i).RF_results.SC_Gaus_FullRF_Num_pixels;
                SC_Gaus_FullRF_Abs_Ellipse_Area_vec(i) = RF_Ident(i).RF_results.SC_Gaus_FullRF_Abs_Ellipse_Area;
                SC_Gaus_FullRF_Angle_Major_Axis_vec(i) = RF_Ident(i).RF_results.SC_Gaus_FullRF_Angle_Major_Axis;
                SC_Gaus_FullRF_Axis_eval_Ratio_vec(i)  = RF_Ident(i).RF_results.SC_Gaus_FullRF_Axis_eval_Ratio;
                SC_Gaus_FullRF_Gaussian_mean_mat(i,:)  = [RF_Ident(i).RF_results.SC_Gaus_FullRF_Gaussian_mean(1),RF_Ident(i).RF_results.SC_Gaus_FullRF_Gaussian_mean(2)];
            end
        end
        SC_Gaus_Num_RF_pixels_mat_2     = SC_Gaus_Num_RF_pixels_mat;                                           % Full vector
        SC_Gaus_FullRF_Num_pixels_vec_2 = SC_Gaus_FullRF_Num_pixels_vec;                                       % Full vector
        SC_Gaus_Num_RF_pixels_mat(sum(isnan(SC_Gaus_Num_RF_pixels_mat),2)==p.Spectral_Dim,:)             = []; % Vector with NaN rows removed
        SC_Gaus_FullRF_Num_pixels_vec(isnan(SC_Gaus_FullRF_Num_pixels_vec))                              = []; % Vector with NaN rows removed
        SC_Gaus_Abs_Ellipse_Area_mat(sum(isnan(SC_Gaus_Abs_Ellipse_Area_mat),2)==p.Spectral_Dim,:)       = [];
        SC_Gaus_FullRF_Abs_Ellipse_Area_vec(isnan(SC_Gaus_FullRF_Abs_Ellipse_Area_vec))                  = [];
        SC_Gaus_RF_Angle_Major_Axis_mat(sum(isnan(SC_Gaus_RF_Angle_Major_Axis_mat),2)==p.Spectral_Dim,:) = [];
        SC_Gaus_FullRF_Angle_Major_Axis_vec(isnan(SC_Gaus_FullRF_Angle_Major_Axis_vec))                  = [];   
        SC_Gaus_Axis_eval_Ratio_mat(sum(isnan(SC_Gaus_Axis_eval_Ratio_mat),2)==p.Spectral_Dim,:)         = [];
        SC_Gaus_FullRF_Axis_eval_Ratio_vec(isnan(SC_Gaus_FullRF_Axis_eval_Ratio_vec))                    = [];   
        SC_Gaus_Gaussian_mean_arr(sum(isnan(SC_Gaus_Gaussian_mean_arr(:,:,1)),2)==p.Spectral_Dim,:,:)    = [];
        SC_Gaus_FullRF_Gaussian_mean_mat(isnan(SC_Gaus_FullRF_Gaussian_mean_mat(:,1)),:)                 = [];
        
    end
    
end

    

%% Plot data

%% STA-SD
if contains(add_info.settings.kernel_new.folderplot{1},'SS') % STA-SD
    
    %% Box
%     if p.RF_Type(1) == 1 % Box
%         
%     end
    
    %% ASP
    if add_info.settings.kernel_new.plot_Allpixel   % All Significant Pixels
        
        %% Num Pixels
        
        % Specrtal RF
        figure;
        for i = 1:double(p.Spectral_Dim)
            subplot(2,2,i);
            hist_loop = histogram(Num_STASD_ASP_RF_Pixels_mat(:,i));
            if Hist_Colour_Choice == 2
                set(hist_loop,'FaceColor',colorMap_arr(end,:,i));
            end
            if i>2
                xlabel('num. pixels');
            end
            if i==1||i==3
                ylabel('num. cells');
            end
            title(Spectral_Names{i});
            set(gca,'FontSize',12);
            xlim([hist_loop.BinEdges(1) hist_loop.BinEdges(end)]);
            if exist('hist_loop.Values','var')
                ylim([0 max(hist_loop.Values)]);
            end
        end
        annotation('textbox',[.5 .975 0 0],'String','STA-SD -- ASP RF: Num. Pixels','FitBoxToText','on','FontSize',12,'FontWeight','bold','LineStyle','none','HorizontalAlignment','center','VerticalAlignment','middle');
        set(gcf,'color','w'); 
        
        % Full RF
        figure;
        subplot(1,2,1);
        hist = histogram(Num_STASD_ASP_FullRF_Pixels_vec);
        xlabel('num. pixels');
        ylabel('num. cells');
        title('Full');
        set(gca,'FontSize',12);
        xlim([hist.BinEdges(1) hist.BinEdges(end)]);
        ylim([0 max(hist.Values)]);
        subplot(1,2,2);
        hist_loop_arr = cell(p.Spectral_Dim,1);
        bin_width_vec = NaN(p.Spectral_Dim,1);
        for i = 1:p.Spectral_Dim
            hist_loop_arr{i} = histogram(Num_STASD_ASP_RF_Pixels_mat(:,i),'FaceAlpha',0.5,'DisplayName',Spectral_Names{i},'Normalization','probability'); hold on;
            bin_width_vec(i) = hist_loop_arr{i}.BinWidth;
            if Hist_Colour_Choice == 2
                set(hist_loop_arr{i},'FaceColor',colorMap_arr(end,:,i));
            end
        end
        if Hist_Bin_Align_Choice == 1
            min_bin_width = min(bin_width_vec(bin_width_vec>0));
            for i = 1:p.Spectral_Dim
                hist_loop_arr{i}.BinWidth = min_bin_width;
            end
        end
        subplot(1,2,2);
        xlabel('num. pixels');
        ylabel('pptn. spectral RFs'); % Proportion of RFs of that colour
        legend;
        title('Spectral');
        set(gca,'FontSize',12);
        %xlim([hist_loop.BinEdges(1) hist_loop.BinEdges(end)]);
%         if exist('hist_loop.Values','var')
%             ylim([0 max(hist_loop.Values)]);
%         end
        annotation('textbox',[.5 .975 0 0],'String','STA-SD -- ASP RF: Num. Pixels','FitBoxToText','on','FontSize',12,'FontWeight','bold','LineStyle','none','HorizontalAlignment','center','VerticalAlignment','middle');
        set(gcf,'color','w');
        
    end
    
    %% Gaussian
    if add_info.settings.kernel_new.plot_Gauss % Gaussian
        
        
        %% Number of Pixels
        
        % Spectral RF
        figure;
        for i = 1:double(p.Spectral_Dim)
            subplot(2,2,i);
            hist_loop = histogram(STASD_Gaus_Num_RF_pixels_mat(:,i));
            if Hist_Colour_Choice == 2
                set(hist_loop,'FaceColor',colorMap_arr(end,:,i));
            end
            if i>2
                xlabel('num. pixels');
            end
            if i==1||i==3
                ylabel('num. cells');
            end
            title(Spectral_Names{i});
            set(gca,'FontSize',12);
            xlim([hist_loop.BinEdges(1) hist_loop.BinEdges(end)]);
            if exist('hist_loop.Values','var')
                ylim([0 max(hist_loop.Values)]);
            end
        end
        annotation('textbox',[.5 .975 0 0],'String','STA-SD -- Gaussian RF: Num. Pixels','FitBoxToText','on','FontSize',12,'FontWeight','bold','LineStyle','none','HorizontalAlignment','center','VerticalAlignment','middle');
        set(gcf,'color','w');
        
        % Full RF (use normalisation prob for spectral together)
        figure;
        subplot(1,2,1);
        hist = histogram(STASD_Gaus_FullRF_Num_pixels_vec);
        xlabel('num. pixels');
        ylabel('num. cells');
        title('Full');
        set(gca,'FontSize',12);
        xlim([hist.BinEdges(1) hist.BinEdges(end)]);
        ylim([0 max(hist.Values)]);
        subplot(1,2,2);
        hist_loop_arr = cell(p.Spectral_Dim,1);
        bin_width_vec = NaN(p.Spectral_Dim,1);
        for i = 1:p.Spectral_Dim
            hist_loop_arr{i} = histogram(STASD_Gaus_Num_RF_pixels_mat(STASD_Gaus_Num_RF_pixels_mat(:,i)<=10,i),'FaceAlpha',0.7,'DisplayName',Spectral_Names{i},'Normalization','probability'); hold on; % I have removed unrealistically big RF.
            bin_width_vec(i) = hist_loop_arr{i}.BinWidth;
            if Hist_Colour_Choice == 2
                set(hist_loop_arr{i},'FaceColor',colorMap_arr(end,:,i));
            end
        end
        if Hist_Bin_Align_Choice == 1
            min_bin_width = min(bin_width_vec(bin_width_vec>0));
            for i = 1:p.Spectral_Dim
                hist_loop_arr{i}.BinWidth = min_bin_width;
            end
        end
        subplot(1,2,2);
        xlabel('num. pixels');
        ylabel('pptn. spectral RFs'); % Proportion of RFs of that colour
        legend;
        title('Spectral');
        set(gca,'FontSize',12);
        %xlim([hist_loop.BinEdges(1) hist_loop.BinEdges(end)]);
%         if exist('hist_loop.Values','var')
%             ylim([0 max(hist_loop.Values)]);
%         end
        annotation('textbox',[.5 .975 0 0],'String','STA-SD -- Gaussian RF: Num. Pixels','FitBoxToText','on','FontSize',12,'FontWeight','bold','LineStyle','none','HorizontalAlignment','center','VerticalAlignment','middle');
        set(gcf,'color','w');
        
        
        %% Ellipse Absolute Area
        
        % Spectral RF
        figure;
        for i = 1:double(p.Spectral_Dim)
            subplot(2,2,i);
            hist_loop = histogram(STASD_Gaus_Abs_Ellipse_Area_mat(:,i));
            if Hist_Colour_Choice == 2
                set(hist_loop,'FaceColor',colorMap_arr(end,:,i));
            end
            if i>2
                xlabel('ellipse area (a.u.)');
            end
            if i==1||i==3
                ylabel('num. cells');
            end
            title(Spectral_Names{i});
            set(gca,'FontSize',12);
            xlim([hist_loop.BinEdges(1) hist_loop.BinEdges(end)]);
            if exist('hist_loop.Values','var')
                ylim([0 max(hist_loop.Values)]);
            end
        end
        annotation('textbox',[.5 .975 0 0],'String','STA-SD -- Gaussian RF: Ellipse Area','FitBoxToText','on','FontSize',12,'FontWeight','bold','LineStyle','none','HorizontalAlignment','center','VerticalAlignment','middle');
        set(gcf,'color','w');
        
        % Full RF (use normalisation prob for spectral together)
        figure;
        subplot(1,2,1);
        hist = histogram(STASD_Gaus_FullRF_Abs_Ellipse_Area_vec);
        xlabel('ellipse area (a.u.)');
        ylabel('num. cells');
        title('Full');
        set(gca,'FontSize',12);
        xlim([hist.BinEdges(1) hist.BinEdges(end)]);
        ylim([0 max(hist.Values)]);
        subplot(1,2,2);
        hist_loop_arr = cell(p.Spectral_Dim,1);
        bin_width_vec = NaN(p.Spectral_Dim,1);
        for i = 1:p.Spectral_Dim
            hist_loop_arr{i} = histogram(STASD_Gaus_Abs_Ellipse_Area_mat(STASD_Gaus_Num_RF_pixels_mat(:,i)<=10,i),'FaceAlpha',0.7,'DisplayName',Spectral_Names{i},'Normalization','probability'); hold on; % I have removed unrealistically big RF.
            bin_width_vec(i) = hist_loop_arr{i}.BinWidth;
            if Hist_Colour_Choice == 2
                set(hist_loop_arr{i},'FaceColor',colorMap_arr(end,:,i));
            end
        end
        if Hist_Bin_Align_Choice == 1
            min_bin_width = min(bin_width_vec(bin_width_vec>0));
            for i = 1:p.Spectral_Dim
                hist_loop_arr{i}.BinWidth = min_bin_width;
            end
        end
        subplot(1,2,2);
        xlabel('ellipse area (a.u.)');
        ylabel('pptn. spectral RFs'); % Proportion of RFs of that colour
        legend;
        title('Spectral');
        set(gca,'FontSize',12);
        %xlim([hist_loop.BinEdges(1) hist_loop.BinEdges(end)]);
%         if exist('hist_loop.Values','var')
%             ylim([0 max(hist_loop.Values)]);
%         end
        annotation('textbox',[.5 .975 0 0],'String','STA-SD -- Gaussian RF: Ellipse Area','FitBoxToText','on','FontSize',12,'FontWeight','bold','LineStyle','none','HorizontalAlignment','center','VerticalAlignment','middle');
        set(gcf,'color','w');
        
        
        %% Major Axis Angle
        % plot angles only for elliptical RFs (not circular ones)
        
        % Spectral RF
        
        figure;
        for i = 1:p.Spectral_Dim
            subplot(2,2,i);
            hist_loop = polarhistogram(STASD_Gaus_RF_Angle_Major_Axis_mat(STASD_Gaus_Axis_eval_Ratio_mat(:,i)~=1,i),(pi/180)*[-90 -67.5 -45 -22.5 0 22.5 45 67.5 90]);%,'DisplayName',Spectral_Names{i},'Normalization','probability'
            if Hist_Colour_Choice == 2
                set(hist_loop,'FaceColor',colorMap_arr(end,:,i)); % ,'ThetaLim',[-180 180]
            end
            title(Spectral_Names{i});
            set(gca,'ThetaTick',[0 45 90 270 315]);
            set(gca,'ThetaTickLabel',{'0';'45';'90';'-90';'-45'});
            set(gca,'FontSize',12);
            %hist_loop.ThetaLim = [-180 180];
        end
        %legend;
        annotation('textbox',[.5 .975 0 0],'String','STA-SD -- Gaussian RF: Major Axis Angle','FitBoxToText','on','FontSize',12,'FontWeight','bold','LineStyle','none','HorizontalAlignment','center','VerticalAlignment','middle');
        set(gcf,'color','w');
        %set(gca,'ThetaZeroLocation','right');
        
        
        % Full RF
        
        figure;
        subplot(1,2,1);
        polarhistogram(STASD_Gaus_FullRF_Angle_Major_Axis_vec(STASD_Gaus_FullRF_Axis_eval_Ratio_vec~=1),(pi/180)*[-90 -67.5 -45 -22.5 0 22.5 45 67.5 90]); % ,'FaceColor','m', 'Normalization','probability'
        title('Full'); %,'Interpreter','Latex','FontName','Helvetica','FontWeight','bold'
        legend('Full RF');
        set(gca,'ThetaTick',[0 45 90 270 315]);
        set(gca,'ThetaTickLabel',{'0';'45';'90';'-90';'-45'});
        set(gca,'FontSize',12);
        subplot(1,2,2);
        hist_loop_arr = cell(p.Spectral_Dim,1);
        bin_width_vec = NaN(p.Spectral_Dim,1);
        for i = 1:p.Spectral_Dim%1:p.Spectral_Dim
            hist_loop_arr{i} = polarhistogram(STASD_Gaus_RF_Angle_Major_Axis_mat(STASD_Gaus_Axis_eval_Ratio_mat(:,i)~=1,i),(pi/180)*[-90 -67.5 -45 -22.5 0 22.5 45 67.5 90],'Normalization','probability','FaceAlpha',0.7,'DisplayName',Spectral_Names{i}); hold on; %,'DisplayName',Spectral_Names{i}
            bin_width_vec(i) = hist_loop_arr{i}.BinWidth;
            if Hist_Colour_Choice == 2
                set(hist_loop_arr{i},'FaceColor',colorMap_arr(end,:,i));
            end
            set(gca,'ThetaTick',[0 45 90 270 315]);
            set(gca,'ThetaTickLabel',{'0';'45';'90';'-90';'-45'});
            set(gca,'FontSize',12);
        end
%         if Hist_Bin_Align_Choice == 1
%             min_bin_width = min(bin_width_vec(bin_width_vec>0));
%             for i = 1:p.Spectral_Dim
%                 hist_loop_arr{i}.BinWidth = min_bin_width;
%             end
%         end
        title('Spectral');
        legend;
        annotation('textbox',[.5 .975 0 0],'String','STA-SD -- Gaussian RF: Major Axis Angle','FitBoxToText','on','FontSize',12,'FontWeight','bold','LineStyle','none','HorizontalAlignment','center','VerticalAlignment','middle');
        set(gcf,'color','w');
        
        
        %% RF Ellipticity (Axis Eval Ratio)
        
        % Spectral RF
        figure;
        for i = 1:p.Spectral_Dim
            subplot(2,2,i);
            hist_loop = histogram(STASD_Gaus_Axis_eval_Ratio_mat(STASD_Gaus_Axis_eval_Ratio_mat(:,i)<=20,i)); % I have removed RF ellipticities > 20.
            if Hist_Colour_Choice == 2
                set(hist_loop,'FaceColor',colorMap_arr(end,:,i));
            end
            if i>2
                xlabel('RF ellipticity');
            end
            if i==1||i==3
                ylabel('num. cells');
            end
            title(Spectral_Names{i});
            set(gca,'FontSize',12);
            xlim([hist_loop.BinEdges(1) hist_loop.BinEdges(end)]);
%             if exist('hist_loop.Values','var')
%                 ylim([0 max(hist_loop.Values)]);
%             end
        end
        annotation('textbox',[.5 .975 0 0],'String','STA-SD -- Gaussian RF: Ellipticity','FitBoxToText','on','FontSize',12,'FontWeight','bold','LineStyle','none','HorizontalAlignment','center','VerticalAlignment','middle');
        set(gcf,'color','w');
        
        % Full RF (use normalisation prob for spectral together)
        figure;
        subplot(1,2,1);
        hist = histogram(STASD_Gaus_FullRF_Axis_eval_Ratio_vec(STASD_Gaus_FullRF_Axis_eval_Ratio_vec<=20)); % I have removed RF ellipticities > 20.
        xlabel('RF ellipticity');
        ylabel('num. cells');
        title('Full');
        set(gca,'FontSize',12);
        xlim([hist.BinEdges(1) hist.BinEdges(end)]);
        ylim([0 max(hist.Values)]);
        subplot(1,2,2);
        hist_loop_arr = cell(p.Spectral_Dim,1);
        bin_width_vec = NaN(p.Spectral_Dim,1);
        for i = 1:p.Spectral_Dim
            hist_loop_arr{i} = histogram(STASD_Gaus_Axis_eval_Ratio_mat(STASD_Gaus_Axis_eval_Ratio_mat(:,i)<=20,i),'FaceAlpha',0.7,'DisplayName',Spectral_Names{i},'Normalization','probability'); hold on; % I have removed RF ellipticities > 20.
            bin_width_vec(i) = hist_loop_arr{i}.BinWidth;
            if Hist_Colour_Choice == 2
                set(hist_loop_arr{i},'FaceColor',colorMap_arr(end,:,i));
            end
        end
        if Hist_Bin_Align_Choice == 1
            min_bin_width = min(bin_width_vec(bin_width_vec>0));
            for i = 1:p.Spectral_Dim
                hist_loop_arr{i}.BinWidth = min_bin_width;
            end
        end
        subplot(1,2,2);
        xlabel('RF ellipticity');
        ylabel('pptn. spectral RFs'); % Proportion of RFs of that colour
        legend;
        title('Spectral');
        set(gca,'FontSize',12);
        %xlim([hist_loop.BinEdges(1) hist_loop.BinEdges(end)]);
%         if exist('hist_loop.Values','var')
%             ylim([0 max(hist_loop.Values)]);
%         end
        annotation('textbox',[.5 .975 0 0],'String','STA-SD -- Gaussian RF: Ellipticity','FitBoxToText','on','FontSize',12,'FontWeight','bold','LineStyle','none','HorizontalAlignment','center','VerticalAlignment','middle');
        set(gcf,'color','w');
        
        
        %% RF Ellipses
        
        % Spectral RF
        ellip_var_vec = linspace(0,2*pi,101); % 101
        figure;
        %tic;
        for j = 1:p.Spectral_Dim
            for i = 1:true_num_cells
                subplot(2,2,j);
                if ~isnan(STASD_Gaus_Num_RF_pixels_mat_2(i,j))
                    mu_1_loop             = RF_Ident(i).RF_results.STASD_Gaus_Gaussian_mean{j}(1);
                    mu_2_loop             = RF_Ident(i).RF_results.STASD_Gaus_Gaussian_mean{j}(2);
                    covar_loop            = RF_Ident(i).RF_results.STASD_Gaus_Gaussian_covar{j};
                    eval_descend_vec_loop = RF_Ident(i).RF_results.STASD_Gaus_eval_descend_vec{j};
                    evec_descend_mat_loop = RF_Ident(i).RF_results.STASD_Gaus_evec_descend_mat{j};
                    sigma_x_loop = sqrt(eval_descend_vec_loop(1));
                    sigma_y_loop = sqrt(eval_descend_vec_loop(2));
                    x_vec_loop = p.RF_SDs*(sigma_x_loop*evec_descend_mat_loop(1,1)*cos(ellip_var_vec) + sigma_y_loop*evec_descend_mat_loop(1,2)*sin(ellip_var_vec)) + mu_1_loop;
                    y_vec_loop = p.RF_SDs*(sigma_x_loop*evec_descend_mat_loop(2,1)*cos(ellip_var_vec) + sigma_y_loop*evec_descend_mat_loop(2,2)*sin(ellip_var_vec)) + mu_2_loop;
                    plot(x_vec_loop,y_vec_loop,'Color',colorMap_arr(end,:,j),'LineWidth',1.5); hold on;
                    %fcontour(RF_Ident{i}.STASD_Gaus_gmPDF{j},[0.5 (p.stim_columns+0.5) 0.5 (p.stim_rows+0.5)],'LineColor',colorMap_arr(end,:,j),'LineWidth',1.5,'LevelList',RF_Ident{i}.STASD_Gaus_Thresh_height{j}); hold on;
                    % 'MeshDensity' default 71, 'MeshDensity',500
                end
            end
            axis equal;
            xlim([0.5 (p.stim_columns+0.5)]);
            ylim([0.5 (p.stim_columns+0.5)]);
            if j>2
                xlabel('x');
            end
            if j==1||j==3
                ylabel('y');
            end
            title(Spectral_Names{j});
            set(gca,'FontSize',12);
            set(gca,'Ydir','reverse');
        end
        %toc; % Elapsed time is 22.532888 seconds with plotting over all cells in Marvin chicken data 1.
        annotation('textbox',[.5 .975 0 0],'String','STA-SD -- Gaussian RF: SD Contours','FitBoxToText','on','FontSize',12,'FontWeight','bold','LineStyle','none','HorizontalAlignment','center','VerticalAlignment','middle');
        set(gcf,'color','w');
        
        
        % Full RF (and Spectral Together)
        ellip_var_vec = linspace(0,2*pi,101); % 101
        figure;
        %tic;
        for i = 1:true_num_cells
            subplot(1,2,1);
            if ~isnan(STASD_Gaus_FullRF_Num_pixels_vec_2(i))
                mu_1_loop  = RF_Ident(i).RF_results.STASD_Gaus_FullRF_Gaussian_mean(1);
                mu_2_loop  = RF_Ident(i).RF_results.STASD_Gaus_FullRF_Gaussian_mean(2);
                covar_loop = RF_Ident(i).RF_results.STASD_Gaus_FullRF_Gaussian_covar;
                eval_descend_vec_loop = RF_Ident(i).RF_results.STASD_Gaus_FullRF_eval_descend_vec;
                evec_descend_mat_loop = RF_Ident(i).RF_results.STASD_Gaus_FullRF_evec_descend_mat;
                sigma_x_loop = sqrt(eval_descend_vec_loop(1));
                sigma_y_loop = sqrt(eval_descend_vec_loop(2));
                x_vec_loop = p.RF_SDs*(sigma_x_loop*evec_descend_mat_loop(1,1)*cos(ellip_var_vec) + sigma_y_loop*evec_descend_mat_loop(1,2)*sin(ellip_var_vec)) + mu_1_loop;
                y_vec_loop = p.RF_SDs*(sigma_x_loop*evec_descend_mat_loop(2,1)*cos(ellip_var_vec) + sigma_y_loop*evec_descend_mat_loop(2,2)*sin(ellip_var_vec)) + mu_2_loop;
                plot(x_vec_loop,y_vec_loop,'Color','k','LineWidth',1.5); hold on;
            end
            subplot(1,2,2);
            for j = 1:p.Spectral_Dim
                if ~isnan(STASD_Gaus_Num_RF_pixels_mat_2(i,j))
                    mu_1_loop  = RF_Ident(i).RF_results.STASD_Gaus_Gaussian_mean{j}(1);
                    mu_2_loop  = RF_Ident(i).RF_results.STASD_Gaus_Gaussian_mean{j}(2);
                    covar_loop = RF_Ident(i).RF_results.STASD_Gaus_Gaussian_covar{j};
                    eval_descend_vec_loop = RF_Ident(i).RF_results.STASD_Gaus_eval_descend_vec{j};
                    evec_descend_mat_loop = RF_Ident(i).RF_results.STASD_Gaus_evec_descend_mat{j};
                    sigma_x_loop = sqrt(eval_descend_vec_loop(1));
                    sigma_y_loop = sqrt(eval_descend_vec_loop(2));
                    x_vec_loop = p.RF_SDs*(sigma_x_loop*evec_descend_mat_loop(1,1)*cos(ellip_var_vec) + sigma_y_loop*evec_descend_mat_loop(1,2)*sin(ellip_var_vec)) + mu_1_loop;
                    y_vec_loop = p.RF_SDs*(sigma_x_loop*evec_descend_mat_loop(2,1)*cos(ellip_var_vec) + sigma_y_loop*evec_descend_mat_loop(2,2)*sin(ellip_var_vec)) + mu_2_loop;
                    plot(x_vec_loop,y_vec_loop,'Color',colorMap_arr(end,:,j),'LineWidth',1.5); hold on;
                    %fcontour(RF_Ident{i}.STASD_Gaus_gmPDF{j},[0.5 (p.stim_columns+0.5) 0.5 (p.stim_rows+0.5)],'LineColor',colorMap_arr(end,:,j),'LineWidth',1.5,'LevelList',RF_Ident{i}.STASD_Gaus_Thresh_height{j}); hold on;
                    % 'MeshDensity' default 71, 'MeshDensity',500
                end
                
            end
        end
        %toc; % Elapsed time is 12.102378 seconds with plotting over all cells in Marvin chicken data 1.
        subplot(1,2,1);
        axis equal;
        xlim([0.5 (p.stim_columns+0.5)]);
        ylim([0.5 (p.stim_columns+0.5)]);
        xlabel('x');
        ylabel('y');
        title('Full');
        set(gca,'Ydir','reverse');
        set(gca,'FontSize',12);
        subplot(1,2,2);
        axis equal;
        xlim([0.5 (p.stim_columns+0.5)]);
        ylim([0.5 (p.stim_columns+0.5)]);
        xlabel('x');
        title('Spectral');
        set(gca,'Ydir','reverse');
        set(gca,'FontSize',12);
        annotation('textbox',[.5 .975 0 0],'String','STA-SD -- Gaussian RF: SD Contours','FitBoxToText','on','FontSize',12,'FontWeight','bold','LineStyle','none','HorizontalAlignment','center','VerticalAlignment','middle');
        set(gcf,'color','w');
        % axis equal;  --> equal units (lengths may be different)
        % axis tight;  --> fit window to data range
        % axis square; --> equal axis lengths (units may be different)
        
        
        %% Mean RF Locations
        
        % Spectral RF
        figure;
        %tic; % It's slightly faster to have the for loops this way around
        for j = 1:p.Spectral_Dim
            for i = 1:true_num_cells
                subplot(2,2,j);
                if ~isnan(STASD_Gaus_Num_RF_pixels_mat_2(i,j))
                    plot(RF_Ident(i).RF_results.STASD_Gaus_Gaussian_mean{j}(1),RF_Ident(i).RF_results.STASD_Gaus_Gaussian_mean{j}(2),'Color',colorMap_arr(end,:,j),'Marker','x','LineWidth',1.5); hold on;
                end
            end
            axis equal;
            xlim([0.5 (p.stim_columns+0.5)]);
            ylim([0.5 (p.stim_columns+0.5)]);
            if j>2
                xlabel('x');
            end
            if j==1||j==3
                ylabel('y');
            end
            title(Spectral_Names{j});
            set(gca,'Ydir','reverse');
            set(gca,'FontSize',12);
        end
        %toc; % Elapsed time is 22.773918 seconds with plotting over all cells in Marvin chicken data 1.
        annotation('textbox',[.5 .975 0 0],'String','STA-SD -- Gaussian RF: Means','FitBoxToText','on','FontSize',12,'FontWeight','bold','LineStyle','none','HorizontalAlignment','center','VerticalAlignment','middle');
        set(gcf,'color','w');
        
        
        % Full RF (and Spectral Together)
        figure;
        %tic; % It's slightly faster to have the for loops this way around
        for i = 1:true_num_cells
            subplot(1,2,1);
            if ~isnan(STASD_Gaus_FullRF_Num_pixels_vec_2(i))
                plot(RF_Ident(i).RF_results.STASD_Gaus_FullRF_Gaussian_mean(1),RF_Ident(i).RF_results.STASD_Gaus_FullRF_Gaussian_mean(2),'Color','k','Marker','x','LineWidth',1.5); hold on;
            end
            for j = 1:p.Spectral_Dim
                subplot(1,2,2);
                if ~isnan(STASD_Gaus_Num_RF_pixels_mat_2(i,j))
                    plot(RF_Ident(i).RF_results.STASD_Gaus_Gaussian_mean{j}(1),RF_Ident(i).RF_results.STASD_Gaus_Gaussian_mean{j}(2),'Color',colorMap_arr(end,:,j),'Marker','x','LineWidth',1.5); hold on;
                end
            end
        end
        %toc; %
        subplot(1,2,1);
        axis equal;
        xlim([0.5 (p.stim_columns+0.5)]);
        ylim([0.5 (p.stim_columns+0.5)]);
        xlabel('x');
        ylabel('y');
        title('Full');
        set(gca,'Ydir','reverse');
        set(gca,'FontSize',12);
        subplot(1,2,2);
        axis equal;
        xlim([0.5 (p.stim_columns+0.5)]);
        ylim([0.5 (p.stim_columns+0.5)]);
        xlabel('x');
        title('Spectral');
        set(gca,'Ydir','reverse');
        set(gca,'FontSize',12);
        annotation('textbox',[.5 .975 0 0],'String','STA-SD -- Gaussian RF: Means','FitBoxToText','on','FontSize',12,'FontWeight','bold','LineStyle','none','HorizontalAlignment','center','VerticalAlignment','middle');
        set(gcf,'color','w');
        
        
        %% Mean Location Density (Histograms and Kernel Density Smoothed)
        
        % Spectral RF - Histogram
        figure;
        if Hist_Colour_Choice == 1
            colormap gray;
        end
        for j = 1:p.Spectral_Dim
            ax = subplot(2,2,j);
            histogram2(STASD_Gaus_Gaussian_mean_arr(:,j,1),STASD_Gaus_Gaussian_mean_arr(:,j,2),...
                [10 10],...
                'XBinLimits',[0.5 (p.stim_columns+0.5)],...
                'YBinLimits',[0.5 (p.stim_columns+0.5)],...
                'Normalization','probability',...
                'DisplayStyle','tile',...
                'ShowEmptyBins','on');
            if Hist_Colour_Choice == 2
                colormap(ax,colorMap_arr(:,:,j));
            end
            if j>2
                xlabel('x');
            end
            if j==1||j==3
                ylabel('y');
            end
            axis equal;
            colorbar;
            title(Spectral_Names{j});
            set(gca,'Ydir','reverse');
            set(gca,'FontSize',12);
        end
        annotation('textbox',[.5 .975 0 0],'String','STA-SD -- Gaussian RF: Mean Dist.','FitBoxToText','on','FontSize',12,'FontWeight','bold','LineStyle','none','HorizontalAlignment','center','VerticalAlignment','middle');
        set(gcf,'color','w');
        
        
        % Full RF - Histogram
        figure;
        colormap gray;
        histogram2(STASD_Gaus_FullRF_Gaussian_mean_mat(:,1),STASD_Gaus_FullRF_Gaussian_mean_mat(:,2),...
            [10 10],...
            'XBinLimits',[0.5 (p.stim_columns+0.5)],...
            'YBinLimits',[0.5 (p.stim_columns+0.5)],...
            'Normalization','probability',...
            'DisplayStyle','tile',...
            'ShowEmptyBins','on'); hold on;
        xlabel('x');
        ylabel('y');
        axis equal;
        colorbar;
        title('STA-SD -- Gaussian RF: Full Mean Dist.');
        set(gca,'Ydir','reverse');
        set(gca,'FontSize',12);
        set(gcf,'color','w');
        
        
        % Spectral RF - Kernel Density Smoothed
        KernSmth_xGrid = linspace(0.5,(p.stim_columns+0.5),41);
        KernSmth_yGrid = linspace(0.5,(p.stim_columns+0.5),41);
        KernSmth_xGrid_length = length(KernSmth_xGrid);
        KernSmth_yGrid_length = length(KernSmth_yGrid);
        Bwidth_Choice  = 2; % Option 1: default, Option 2: specified
        Bwidth_x       = 1; % 1, 2
        Bwidth_y       = 1; % 1, 2
        [x1,x2]        = meshgrid(KernSmth_xGrid, KernSmth_yGrid);
        x1             = x1(:);
        x2             = x2(:);
        xi             = [x1 x2];
        KernSmth_2D = cell(p.Spectral_Dim,1);
        for j = 1:p.Spectral_Dim
            if Bwidth_Choice == 1
                [KernSmth,~,bw] = ksdensity([STASD_Gaus_Gaussian_mean_arr(:,j,1),STASD_Gaus_Gaussian_mean_arr(:,j,2)],xi);
            else
                [KernSmth,~,bw] = ksdensity([STASD_Gaus_Gaussian_mean_arr(:,j,1),STASD_Gaus_Gaussian_mean_arr(:,j,2)],xi,'Bandwidth',[Bwidth_x Bwidth_y]);
            end
            KernSmth_2D{j} = reshape(KernSmth,[KernSmth_xGrid_length,KernSmth_yGrid_length]);
        end
        figure;
        if Hist_Colour_Choice == 1
            colormap gray;
        end
        for j = 1:p.Spectral_Dim
            ax = subplot(2,2,j);
            imagesc(KernSmth_xGrid,KernSmth_yGrid,KernSmth_2D{j});
            %set(gca,'YDir','normal'); % don't reverse, correct as is
            if Hist_Colour_Choice == 2
                colormap(ax,colorMap_arr(:,:,j));
            end
            if j>2
                xlabel('x');
            end
            if j==1||j==3
                ylabel('y');
            end
            axis equal;
            xlim([0.5 (p.stim_columns+0.5)]);
            ylim([0.5 (p.stim_columns+0.5)]);
            colorbar;
            title(Spectral_Names{j});
            set(gca,'FontSize',12);
        end
        annotation('textbox',[.5 .975 0 0],'String','STA-SD -- Gaussian RF: Mean Dist.','FitBoxToText','on','FontSize',12,'FontWeight','bold','LineStyle','none','HorizontalAlignment','center','VerticalAlignment','middle');
        set(gcf,'color','w');
        
        
        % Full RF - Kernel Density Smoothed
        KernSmth_xGrid = linspace(0.5,(p.stim_columns+0.5),41);
        KernSmth_yGrid = linspace(0.5,(p.stim_columns+0.5),41);
        KernSmth_xGrid_length = length(KernSmth_xGrid);
        KernSmth_yGrid_length = length(KernSmth_yGrid);
        Bwidth_Choice  = 2; % Option 1: default, Option 2: specified
        Bwidth_x       = 1; % 1, 2
        Bwidth_y       = 1; % 1, 2
        [x1,x2]        = meshgrid(KernSmth_xGrid, KernSmth_yGrid);
        x1             = x1(:);
        x2             = x2(:);
        xi             = [x1 x2];
        if Bwidth_Choice == 1
            [KernSmth,~,bw] = ksdensity([STASD_Gaus_FullRF_Gaussian_mean_mat(:,1),STASD_Gaus_FullRF_Gaussian_mean_mat(:,2)],xi);
        else
            [KernSmth,~,bw] = ksdensity([STASD_Gaus_FullRF_Gaussian_mean_mat(:,1),STASD_Gaus_FullRF_Gaussian_mean_mat(:,2)],xi,'Bandwidth',[Bwidth_x Bwidth_y]);
        end
        KernSmth_2D = reshape(KernSmth,[KernSmth_xGrid_length,KernSmth_yGrid_length]);
        figure;
        colormap gray;
        imagesc(KernSmth_xGrid,KernSmth_yGrid,KernSmth_2D);
        %set(gca,'YDir','normal'); % don't reverse, correct as is
        xlabel('x');
        ylabel('y');
        axis equal;
        xlim([0.5 (p.stim_columns+0.5)]);
        ylim([0.5 (p.stim_columns+0.5)]);
        colorbar;
        title('STA-SD -- Gaussian RF: Full Mean Dist.');
        set(gca,'FontSize',12);
        set(gcf,'color','w');
        
    end
    
end


% %% LC
% if     p.RF_Ident_Meth_vec(2) == 1 % LC
%     
%     if p.RF_Type(1) == 1 % Box
%         
%     end
%     
%     if p.RF_Type(2) == 1     % All Significant Pixels
%         
%     end
%     
%     if p.RF_Type(3) == 1     % Gaussian
%         
%     end
%     
% end
% 
% %% MI
% if     p.RF_Ident_Meth_vec(3) == 1 % MI
%     
%     if p.RF_Type(1) == 1 % Box
%         
%     end
%     
%     if p.RF_Type(2) == 1     % All Significant Pixels
%         
%     end
%     
%     if p.RF_Type(3) == 1     % Gaussian
%         
%     end
%     
% end

%% SC
if contains(add_info.settings.kernel_new.folderplot{1},'SC')
    
%     if p.RF_Type(1) == 1 % Box
%         
%     end
    
    if add_info.settings.kernel_new.plot_Allpixel% All Significant Pixels
        
        %% Num Pixels
        
        % Spectal RF
        figure;
        for i = 1:p.Spectral_Dim
            subplot(2,2,i);
            hist_loop = histogram(Num_SC_ASP_RF_Pixels_mat(:,i));
            if Hist_Colour_Choice == 2
                set(hist_loop,'FaceColor',colorMap_arr(end,:,i));
            end
            if i>2
                xlabel('num. pixels');
            end
            if i==1||i==3
                ylabel('num. cells');
            end
            title(Spectral_Names{i});
            set(gca,'FontSize',12);
            xlim([hist_loop.BinEdges(1) hist_loop.BinEdges(end)]);
            if exist('hist_loop.Values','var')
                ylim([0 max(hist_loop.Values)]);
            end
        end
        annotation('textbox',[.5 .975 0 0],'String','SC -- ASP RF: Num. Pixels','FitBoxToText','on','FontSize',12,'FontWeight','bold','LineStyle','none','HorizontalAlignment','center','VerticalAlignment','middle');
        set(gcf,'color','w'); 
        
        % Full RF
        figure;
        subplot(1,2,1);
        hist = histogram(Num_SC_ASP_FullRF_Pixels_vec);
        xlabel('num. pixels');
        ylabel('num. cells');
        title('Full');
        set(gca,'FontSize',12);
        xlim([hist.BinEdges(1) hist.BinEdges(end)]);
        ylim([0 max(hist.Values)]);
        subplot(1,2,2);
        hist_loop_arr = cell(p.Spectral_Dim,1);
        bin_width_vec = NaN(p.Spectral_Dim,1);
        for i = 1:p.Spectral_Dim
            hist_loop_arr{i} = histogram(Num_SC_ASP_RF_Pixels_mat(:,i),'FaceAlpha',0.5,'DisplayName',Spectral_Names{i},'Normalization','probability'); hold on;
            bin_width_vec(i) = hist_loop_arr{i}.BinWidth;
            if Hist_Colour_Choice == 2
                set(hist_loop_arr{i},'FaceColor',colorMap_arr(end,:,i));
            end
        end
        if Hist_Bin_Align_Choice == 1
            min_bin_width = min(bin_width_vec(bin_width_vec>0));
            for i = 1:p.Spectral_Dim
                hist_loop_arr{i}.BinWidth = min_bin_width;
            end
        end
        subplot(1,2,2);
        xlabel('num. pixels');
        ylabel('pptn. spectral RFs'); % Proportion of RFs of that colour
        legend;
        title('Spectral');
        set(gca,'FontSize',12);
        %xlim([hist_loop.BinEdges(1) hist_loop.BinEdges(end)]);
%         if exist('hist_loop.Values','var')
%             ylim([0 max(hist_loop.Values)]);
%         end
        annotation('textbox',[.5 .975 0 0],'String','SC -- ASP RF: Num. Pixels','FitBoxToText','on','FontSize',12,'FontWeight','bold','LineStyle','none','HorizontalAlignment','center','VerticalAlignment','middle');
        set(gcf,'color','w');
        
    end
    
    if add_info.settings.kernel_new.plot_Gauss    % Gaussian
        
        %% Number of Pixels
        
        % Spectral RF
        figure;
        for i = 1:p.Spectral_Dim
            subplot(2,2,i);
            hist_loop = histogram(SC_Gaus_Num_RF_pixels_mat(SC_Gaus_Num_RF_pixels_mat(:,i)<=10,i)); % I have removed unrealistically big RF.
            if Hist_Colour_Choice == 2
                set(hist_loop,'FaceColor',colorMap_arr(end,:,i));
            end
            if i>2
                xlabel('num. pixels');
            end
            if i==1||i==3
                ylabel('num. cells');
            end
            title(Spectral_Names{i});
            set(gca,'FontSize',12);
            xlim([hist_loop.BinEdges(1) hist_loop.BinEdges(end)]);
            if exist('hist_loop.Values','var')
                ylim([0 max(hist_loop.Values)]);
            end
        end
        annotation('textbox',[.5 .975 0 0],'String','SC -- Gaussian RF: Num. Pixels','FitBoxToText','on','FontSize',12,'FontWeight','bold','LineStyle','none','HorizontalAlignment','center','VerticalAlignment','middle');
        set(gcf,'color','w');
        
        % Full RF (use normalisation prob for spectral together)
        figure;
        subplot(1,2,1);
        hist = histogram(SC_Gaus_FullRF_Num_pixels_vec(logical((SC_Gaus_FullRF_Num_pixels_vec<=11).*(SC_Gaus_FullRF_Axis_eval_Ratio_vec<=10)))); % I have removed unrealistically big RF sizes and ellipticities.
        xlabel('num. pixels');
        ylabel('num. cells');
        title('Full');
        set(gca,'FontSize',12);
        xlim([hist.BinEdges(1) hist.BinEdges(end)]);
        ylim([0 max(hist.Values)]);
        subplot(1,2,2);
        hist_loop_arr = cell(p.Spectral_Dim,1);
        bin_width_vec = NaN(p.Spectral_Dim,1);
        for i = 1:p.Spectral_Dim
            hist_loop_arr{i} = histogram(SC_Gaus_Num_RF_pixels_mat(SC_Gaus_Num_RF_pixels_mat(:,i)<=10,i),'FaceAlpha',0.7,'DisplayName',Spectral_Names{i},'Normalization','probability'); hold on; % I have removed unrealistically big RF.
            bin_width_vec(i) = hist_loop_arr{i}.BinWidth;
            if Hist_Colour_Choice == 2
                set(hist_loop_arr{i},'FaceColor',colorMap_arr(end,:,i));
            end
        end
        if Hist_Bin_Align_Choice == 1
            min_bin_width = min(bin_width_vec(bin_width_vec>0));
            for i = 1:p.Spectral_Dim
                hist_loop_arr{i}.BinWidth = min_bin_width;
            end
        end
        subplot(1,2,2);
        xlabel('num. pixels');
        ylabel('pptn. spectral RFs'); % Proportion of RFs of that colour
        legend;
        title('Spectral');
        set(gca,'FontSize',12);
        %xlim([hist_loop.BinEdges(1) hist_loop.BinEdges(end)]);
%         if exist('hist_loop.Values','var')
%             ylim([0 max(hist_loop.Values)]);
%         end
        annotation('textbox',[.5 .975 0 0],'String','SC -- Gaussian RF: Num. Pixels','FitBoxToText','on','FontSize',12,'FontWeight','bold','LineStyle','none','HorizontalAlignment','center','VerticalAlignment','middle');
        set(gcf,'color','w');
        
        
        %% Ellipse Absolute Area
        
        % Spectral RF
        figure;
        for i = 1:p.Spectral_Dim
            subplot(2,2,i);
            hist_loop = histogram(SC_Gaus_Abs_Ellipse_Area_mat(SC_Gaus_Num_RF_pixels_mat(:,i)<=10,i)); % I have removed unrealistically big RF.
            if Hist_Colour_Choice == 2
                set(hist_loop,'FaceColor',colorMap_arr(end,:,i));
            end
            if i>2
                xlabel('ellipse area (a.u.)');
            end
            if i==1||i==3
                ylabel('num. cells');
            end
            title(Spectral_Names{i});
            set(gca,'FontSize',12);
            xlim([hist_loop.BinEdges(1) hist_loop.BinEdges(end)]);
            if exist('hist_loop.Values','var')
                ylim([0 max(hist_loop.Values)]);
            end
        end
        annotation('textbox',[.5 .975 0 0],'String','SC -- Gaussian RF: Ellipse Area','FitBoxToText','on','FontSize',12,'FontWeight','bold','LineStyle','none','HorizontalAlignment','center','VerticalAlignment','middle');
        set(gcf,'color','w');
        
        % Full RF (use normalisation prob for spectral together)
        figure;
        subplot(1,2,1);
        hist = histogram(SC_Gaus_FullRF_Abs_Ellipse_Area_vec(logical((SC_Gaus_FullRF_Num_pixels_vec<=11).*(SC_Gaus_FullRF_Axis_eval_Ratio_vec<=10)))); % I have removed unrealistically big RF sizes and ellipticities.
        xlabel('ellipse area (a.u.)');
        ylabel('num. cells');
        title('Full');
        set(gca,'FontSize',12);
        xlim([hist.BinEdges(1) hist.BinEdges(end)]);
        ylim([0 max(hist.Values)]);
        subplot(1,2,2);
        hist_loop_arr = cell(p.Spectral_Dim,1);
        bin_width_vec = NaN(p.Spectral_Dim,1);
        for i = 1:p.Spectral_Dim
            hist_loop_arr{i} = histogram(SC_Gaus_Abs_Ellipse_Area_mat(SC_Gaus_Num_RF_pixels_mat(:,i)<=10,i),'FaceAlpha',0.7,'DisplayName',Spectral_Names{i},'Normalization','probability'); hold on; % I have removed unrealistically big RF.
            bin_width_vec(i) = hist_loop_arr{i}.BinWidth;
            if Hist_Colour_Choice == 2
                set(hist_loop_arr{i},'FaceColor',colorMap_arr(end,:,i));
            end
        end
        if Hist_Bin_Align_Choice == 1
            min_bin_width = min(bin_width_vec(bin_width_vec>0));
            for i = 1:p.Spectral_Dim
                hist_loop_arr{i}.BinWidth = min_bin_width;
            end
        end
        subplot(1,2,2);
        xlabel('ellipse area (a.u.)');
        ylabel('pptn. spectral RFs'); % Proportion of RFs of that colour
        legend;
        title('Spectral');
        set(gca,'FontSize',12);
        %xlim([hist_loop.BinEdges(1) hist_loop.BinEdges(end)]);
%         if exist('hist_loop.Values','var')
%             ylim([0 max(hist_loop.Values)]);
%         end
        annotation('textbox',[.5 .975 0 0],'String','SC -- Gaussian RF: Ellipse Area','FitBoxToText','on','FontSize',12,'FontWeight','bold','LineStyle','none','HorizontalAlignment','center','VerticalAlignment','middle');
        set(gcf,'color','w');
        
        
        %% Major Axis Angle
        % plot angles only for elliptical RFs (not circular ones)
        
        % Spectral RF
        
        figure;
        for i = 1:p.Spectral_Dim
            subplot(2,2,i);
            hist_loop = polarhistogram(SC_Gaus_RF_Angle_Major_Axis_mat(SC_Gaus_Axis_eval_Ratio_mat(:,i)~=1,i),(pi/180)*[-90 -67.5 -45 -22.5 0 22.5 45 67.5 90]);%,'DisplayName',Spectral_Names{i},'Normalization','probability'
            if Hist_Colour_Choice == 2
                set(hist_loop,'FaceColor',colorMap_arr(end,:,i)); % ,'ThetaLim',[-180 180]
            end
            title(Spectral_Names{i});
            set(gca,'ThetaTick',[0 45 90 270 315]);
            set(gca,'ThetaTickLabel',{'0';'45';'90';'-90';'-45'});
            set(gca,'FontSize',12);
            %hist_loop.ThetaLim = [-180 180];
        end
        %legend;
        annotation('textbox',[.5 .975 0 0],'String','SC -- Gaussian RF: Major Axis Angle','FitBoxToText','on','FontSize',12,'FontWeight','bold','LineStyle','none','HorizontalAlignment','center','VerticalAlignment','middle');
        set(gcf,'color','w');
        %set(gca,'ThetaZeroLocation','right');
        
        
        % Full RF
        
        figure;
        subplot(1,2,1);
        polarhistogram(SC_Gaus_FullRF_Angle_Major_Axis_vec(SC_Gaus_FullRF_Axis_eval_Ratio_vec~=1),(pi/180)*[-90 -67.5 -45 -22.5 0 22.5 45 67.5 90]); % ,'FaceColor','m', 'Normalization','probability'
        title('Full'); %,'Interpreter','Latex','FontName','Helvetica','FontWeight','bold'
        legend('Full RF');
        set(gca,'ThetaTick',[0 45 90 270 315]);
        set(gca,'ThetaTickLabel',{'0';'45';'90';'-90';'-45'});
        set(gca,'FontSize',12);
        subplot(1,2,2);
        hist_loop_arr = cell(p.Spectral_Dim,1);
        bin_width_vec = NaN(p.Spectral_Dim,1);
        for i = 1:p.Spectral_Dim%1:p.Spectral_Dim
            hist_loop_arr{i} = polarhistogram(SC_Gaus_RF_Angle_Major_Axis_mat(SC_Gaus_Axis_eval_Ratio_mat(:,i)~=1,i),(pi/180)*[-90 -67.5 -45 -22.5 0 22.5 45 67.5 90],'Normalization','probability','FaceAlpha',0.7,'DisplayName',Spectral_Names{i}); hold on; %,'DisplayName',Spectral_Names{i}
            bin_width_vec(i) = hist_loop_arr{i}.BinWidth;
            if Hist_Colour_Choice == 2
                set(hist_loop_arr{i},'FaceColor',colorMap_arr(end,:,i));
            end
            set(gca,'ThetaTick',[0 45 90 270 315]);
            set(gca,'ThetaTickLabel',{'0';'45';'90';'-90';'-45'});
            set(gca,'FontSize',12);
        end
%         if Hist_Bin_Align_Choice == 1
%             min_bin_width = min(bin_width_vec(bin_width_vec>0));
%             for i = 1:p.Spectral_Dim%1:p.Spectral_Dim
%                 hist_loop_arr{i}.BinWidth = min_bin_width;
%             end
%         end
        title('Spectral');
        legend;
        annotation('textbox',[.5 .975 0 0],'String','SC -- Gaussian RF: Major Axis Angle','FitBoxToText','on','FontSize',12,'FontWeight','bold','LineStyle','none','HorizontalAlignment','center','VerticalAlignment','middle');
        set(gcf,'color','w');
        
        
        %% RF Ellipticity (Axis Eval Ratio)
        
        % Spectral RF
        figure;
        for i = 1:p.Spectral_Dim
            subplot(2,2,i);
            hist_loop = histogram(SC_Gaus_Axis_eval_Ratio_mat(SC_Gaus_Axis_eval_Ratio_mat(:,i)<=20,i)); % I have removed RF ellipticities > 20
            if Hist_Colour_Choice == 2
                set(hist_loop,'FaceColor',colorMap_arr(end,:,i));
            end
            if i>2
                xlabel('RF ellipticity');
            end
            if i==1||i==3
                ylabel('num. cells');
            end
            title(Spectral_Names{i});
            set(gca,'FontSize',12);
            xlim([hist_loop.BinEdges(1) hist_loop.BinEdges(end)]);
            if exist('hist_loop.Values','var')
                ylim([0 max(hist_loop.Values)]);
            end
        end
        annotation('textbox',[.5 .975 0 0],'String','SC -- Gaussian RF: Ellipticity','FitBoxToText','on','FontSize',12,'FontWeight','bold','LineStyle','none','HorizontalAlignment','center','VerticalAlignment','middle');
        set(gcf,'color','w');
        
        % Full RF (use normalisation prob for spectral together)
        figure;
        subplot(1,2,1);
        hist = histogram(SC_Gaus_FullRF_Axis_eval_Ratio_vec(SC_Gaus_FullRF_Axis_eval_Ratio_vec<=20)); % I have removed RF ellipticities > 20
        xlabel('RF ellipticity');
        ylabel('num. cells');
        title('Full');
        set(gca,'FontSize',12);
        xlim([hist.BinEdges(1) hist.BinEdges(end)]);
        ylim([0 max(hist.Values)]);
        subplot(1,2,2);
        hist_loop_arr = cell(p.Spectral_Dim,1);
        bin_width_vec = NaN(p.Spectral_Dim,1);
        for i = 1:p.Spectral_Dim
            hist_loop_arr{i} = histogram(SC_Gaus_Axis_eval_Ratio_mat(SC_Gaus_Axis_eval_Ratio_mat(:,i)<=20,i),'FaceAlpha',0.7,'DisplayName',Spectral_Names{i},'Normalization','probability'); hold on; % I have removed RF ellipticities > 20
            bin_width_vec(i) = hist_loop_arr{i}.BinWidth;
            if Hist_Colour_Choice == 2
                set(hist_loop_arr{i},'FaceColor',colorMap_arr(end,:,i));
            end
        end
        if Hist_Bin_Align_Choice == 1
            min_bin_width = min(bin_width_vec(bin_width_vec>0));
            for i = 1:p.Spectral_Dim
                hist_loop_arr{i}.BinWidth = min_bin_width;
            end
        end
        subplot(1,2,2);
        xlabel('RF ellipticity');
        ylabel('pptn. spectral RFs'); % Proportion of RFs of that colour
        legend;
        title('Spectral');
        set(gca,'FontSize',12);
        %xlim([hist_loop.BinEdges(1) hist_loop.BinEdges(end)]);
%         if exist('hist_loop.Values','var')
%             ylim([0 max(hist_loop.Values)]);
%         end
        annotation('textbox',[.5 .975 0 0],'String','SC -- Gaussian RF: Ellipticity','FitBoxToText','on','FontSize',12,'FontWeight','bold','LineStyle','none','HorizontalAlignment','center','VerticalAlignment','middle');
        set(gcf,'color','w');
        
        
        %% RF Ellipses
        
        % Spectral RF
        ellip_var_vec = linspace(0,2*pi,101); % 101
        figure;
        %tic;
        for j = 1:p.Spectral_Dim
            for i = 1:true_num_cells
                subplot(2,2,j);
                if ~isnan(SC_Gaus_Num_RF_pixels_mat_2(i,j))
                    mu_1_loop = RF_Ident(i).RF_results.SC_Gaus_Gaussian_mean{j}(1);
                    mu_2_loop = RF_Ident(i).RF_results.SC_Gaus_Gaussian_mean{j}(2);
                    covar_loop = RF_Ident(i).RF_results.SC_Gaus_Gaussian_covar{j};
                    eval_descend_vec_loop = RF_Ident(i).RF_results.SC_Gaus_eval_descend_vec{j};
                    evec_descend_mat_loop = RF_Ident(i).RF_results.SC_Gaus_evec_descend_mat{j};
                    sigma_x_loop = sqrt(eval_descend_vec_loop(1));
                    sigma_y_loop = sqrt(eval_descend_vec_loop(2));
                    x_vec_loop = p.RF_SDs*(sigma_x_loop*evec_descend_mat_loop(1,1)*cos(ellip_var_vec) + sigma_y_loop*evec_descend_mat_loop(1,2)*sin(ellip_var_vec)) + mu_1_loop;
                    y_vec_loop = p.RF_SDs*(sigma_x_loop*evec_descend_mat_loop(2,1)*cos(ellip_var_vec) + sigma_y_loop*evec_descend_mat_loop(2,2)*sin(ellip_var_vec)) + mu_2_loop;
                    plot(x_vec_loop,y_vec_loop,'Color',colorMap_arr(end,:,j),'LineWidth',1.5); hold on;
                    %fcontour(RF_Ident{i}.SC_Gaus_gmPDF{j},[0.5 (p.stim_columns+0.5) 0.5 (p.stim_rows+0.5)],'LineColor',colorMap_arr(end,:,j),'LineWidth',1.5,'LevelList',RF_Ident{i}.SC_Gaus_Thresh_height{j}); hold on;
                    % 'MeshDensity' default 71, 'MeshDensity',500
                end
            end
            axis equal;
            xlim([0.5 (p.stim_columns+0.5)]);
            ylim([0.5 (p.stim_columns+0.5)]);
            if j>2
                xlabel('x');
            end
            if j==1||j==3
                ylabel('y');
            end
            title(Spectral_Names{j});
            set(gca,'Ydir','reverse');
            set(gca,'FontSize',12);
        end
        %toc; % Elapsed time is 22.532888 seconds with plotting over all cells in Marvin chicken data 1.
        annotation('textbox',[.5 .975 0 0],'String','SC -- Gaussian RF: SD Contours','FitBoxToText','on','FontSize',12,'FontWeight','bold','LineStyle','none','HorizontalAlignment','center','VerticalAlignment','middle');
        set(gcf,'color','w');
        
        
        % Full RF (and Spectral Together)
        ellip_var_vec = linspace(0,2*pi,101); % 101
        figure;
        %tic;
        for i = 1:true_num_cells
            subplot(1,2,1);
            if ~isnan(SC_Gaus_FullRF_Num_pixels_vec_2(i))
                mu_1_loop  = RF_Ident(i).RF_results.SC_Gaus_FullRF_Gaussian_mean(1);
                mu_2_loop  = RF_Ident(i).RF_results.SC_Gaus_FullRF_Gaussian_mean(2);
                covar_loop = RF_Ident(i).RF_results.SC_Gaus_FullRF_Gaussian_covar;
                eval_descend_vec_loop = RF_Ident(i).RF_results.SC_Gaus_FullRF_eval_descend_vec;
                evec_descend_mat_loop = RF_Ident(i).RF_results.SC_Gaus_FullRF_evec_descend_mat;
                sigma_x_loop = sqrt(eval_descend_vec_loop(1));
                sigma_y_loop = sqrt(eval_descend_vec_loop(2));
                x_vec_loop = p.RF_SDs*(sigma_x_loop*evec_descend_mat_loop(1,1)*cos(ellip_var_vec) + sigma_y_loop*evec_descend_mat_loop(1,2)*sin(ellip_var_vec)) + mu_1_loop;
                y_vec_loop = p.RF_SDs*(sigma_x_loop*evec_descend_mat_loop(2,1)*cos(ellip_var_vec) + sigma_y_loop*evec_descend_mat_loop(2,2)*sin(ellip_var_vec)) + mu_2_loop;
                plot(x_vec_loop,y_vec_loop,'Color','k','LineWidth',1.5); hold on;
            end
            subplot(1,2,2);
            for j = 1:p.Spectral_Dim
                if ~isnan(SC_Gaus_Num_RF_pixels_mat_2(i,j))
                    mu_1_loop  = RF_Ident(i).RF_results.SC_Gaus_Gaussian_mean{j}(1);
                    mu_2_loop  = RF_Ident(i).RF_results.SC_Gaus_Gaussian_mean{j}(2);
                    covar_loop = RF_Ident(i).RF_results.SC_Gaus_Gaussian_covar{j};
                    eval_descend_vec_loop = RF_Ident(i).RF_results.SC_Gaus_eval_descend_vec{j};
                    evec_descend_mat_loop = RF_Ident(i).RF_results.SC_Gaus_evec_descend_mat{j};
                    sigma_x_loop = sqrt(eval_descend_vec_loop(1));
                    sigma_y_loop = sqrt(eval_descend_vec_loop(2));
                    x_vec_loop = p.RF_SDs*(sigma_x_loop*evec_descend_mat_loop(1,1)*cos(ellip_var_vec) + sigma_y_loop*evec_descend_mat_loop(1,2)*sin(ellip_var_vec)) + mu_1_loop;
                    y_vec_loop = p.RF_SDs*(sigma_x_loop*evec_descend_mat_loop(2,1)*cos(ellip_var_vec) + sigma_y_loop*evec_descend_mat_loop(2,2)*sin(ellip_var_vec)) + mu_2_loop;
                    plot(x_vec_loop,y_vec_loop,'Color',colorMap_arr(end,:,j),'LineWidth',1.5); hold on;
                    %fcontour(RF_Ident{i}.SC_Gaus_gmPDF{j},[0.5 (p.stim_columns+0.5) 0.5 (p.stim_rows+0.5)],'LineColor',colorMap_arr(end,:,j),'LineWidth',1.5,'LevelList',RF_Ident{i}.SC_Gaus_Thresh_height{j}); hold on;
                    % 'MeshDensity' default 71, 'MeshDensity',500
                end
                
            end
        end
        %toc; % Elapsed time is 12.102378 seconds with plotting over all cells in Marvin chicken data 1.
        subplot(1,2,1);
        axis equal;
        xlim([0.5 (p.stim_columns+0.5)]);
        ylim([0.5 (p.stim_columns+0.5)]);
        xlabel('x');
        ylabel('y');
        title('Full');
        set(gca,'Ydir','reverse');
        set(gca,'FontSize',12);
        subplot(1,2,2);
        axis equal;
        xlim([0.5 (p.stim_columns+0.5)]);
        ylim([0.5 (p.stim_columns+0.5)]);
        xlabel('x');
        title('Spectral');
        set(gca,'Ydir','reverse');
        set(gca,'FontSize',12);
        annotation('textbox',[.5 .975 0 0],'String','SC -- Gaussian RF: SD Contours','FitBoxToText','on','FontSize',12,'FontWeight','bold','LineStyle','none','HorizontalAlignment','center','VerticalAlignment','middle');
        set(gcf,'color','w');
        % axis equal;  --> equal units (lengths may be different)
        % axis tight;  --> fit window to data range
        % axis square; --> equal axis lengths (units may be different)
        
        
        %% Mean RF Locations
        
        % Spectral RF
        figure;
        %tic; % It's slightly faster to have the for loops this way around
        for j = 1:p.Spectral_Dim
            for i = 1:true_num_cells
                subplot(2,2,j);
                if ~isnan(SC_Gaus_Num_RF_pixels_mat_2(i,j))
                    plot(RF_Ident(i).RF_results.SC_Gaus_Gaussian_mean{j}(1),RF_Ident(i).RF_results.SC_Gaus_Gaussian_mean{j}(2),'Color',colorMap_arr(end,:,j),'Marker','x','LineWidth',1.5); hold on;
                end
            end
            axis equal;
            xlim([0.5 (p.stim_columns+0.5)]);
            ylim([0.5 (p.stim_columns+0.5)]);
            if j>2
                xlabel('x');
            end
            if j==1||j==3
                ylabel('y');
            end
            title(Spectral_Names{j});
            set(gca,'Ydir','reverse');
            set(gca,'FontSize',12);
        end
        %toc; % Elapsed time is 22.773918 seconds with plotting over all cells in Marvin chicken data 1.
        annotation('textbox',[.5 .975 0 0],'String','SC -- Gaussian RF: Means','FitBoxToText','on','FontSize',12,'FontWeight','bold','LineStyle','none','HorizontalAlignment','center','VerticalAlignment','middle');
        set(gcf,'color','w');
        
        
        % Full RF (and Spectral Together)
        figure;
        %tic; % It's slightly faster to have the for loops this way around
        for i = 1:true_num_cells
            subplot(1,2,1);
            if ~isnan(SC_Gaus_FullRF_Num_pixels_vec_2(i))
                plot(RF_Ident(i).RF_results.SC_Gaus_FullRF_Gaussian_mean(1),RF_Ident(i).RF_results.SC_Gaus_FullRF_Gaussian_mean(2),'Color','k','Marker','x','LineWidth',1.5); hold on;
            end
            for j = 1:p.Spectral_Dim
                subplot(1,2,2);
                if ~isnan(SC_Gaus_Num_RF_pixels_mat_2(i,j))
                    plot(RF_Ident(i).RF_results.SC_Gaus_Gaussian_mean{j}(1),RF_Ident(i).RF_results.SC_Gaus_Gaussian_mean{j}(2),'Color',colorMap_arr(end,:,j),'Marker','x','LineWidth',1.5); hold on;
                end
            end
        end
        %toc; %
        subplot(1,2,1);
        axis equal;
        xlim([0.5 (p.stim_columns+0.5)]);
        ylim([0.5 (p.stim_columns+0.5)]);
        xlabel('x');
        ylabel('y');
        title('Full');
        set(gca,'Ydir','reverse');
        set(gca,'FontSize',12);
        subplot(1,2,2);
        axis equal;
        xlim([0.5 (p.stim_columns+0.5)]);
        ylim([0.5 (p.stim_columns+0.5)]);
        xlabel('x');
        title('Spectral');
        set(gca,'Ydir','reverse');
        set(gca,'FontSize',12);
        annotation('textbox',[.5 .975 0 0],'String','SC -- Gaussian RF: Means','FitBoxToText','on','FontSize',12,'FontWeight','bold','LineStyle','none','HorizontalAlignment','center','VerticalAlignment','middle');
        set(gcf,'color','w');
        
        
        %% Mean Location Density (Histograms and Kernel Density Smoothed)
        
        % Spectral RF - Histogram
        figure;
        if Hist_Colour_Choice == 1
            colormap gray;
        end
        for j = 1:p.Spectral_Dim
            ax = subplot(2,2,j);
            histogram2(SC_Gaus_Gaussian_mean_arr(:,j,1),SC_Gaus_Gaussian_mean_arr(:,j,2),...
                [10 10],...
                'XBinLimits',[0.5 (p.stim_columns+0.5)],...
                'YBinLimits',[0.5 (p.stim_columns+0.5)],...
                'Normalization','probability',...
                'DisplayStyle','tile',...
                'ShowEmptyBins','on');
            if Hist_Colour_Choice == 2
                colormap(ax,colorMap_arr(:,:,j));
            end
            if j>2
                xlabel('x');
            end
            if j==1||j==3
                ylabel('y');
            end
            axis equal;
            colorbar;
            title(Spectral_Names{j});
            set(gca,'Ydir','reverse');
            set(gca,'FontSize',12);
        end
        annotation('textbox',[.5 .975 0 0],'String','SC -- Gaussian RF: Mean Dist.','FitBoxToText','on','FontSize',12,'FontWeight','bold','LineStyle','none','HorizontalAlignment','center','VerticalAlignment','middle');
        set(gcf,'color','w');
        
        
        % Full RF - Histogram
        figure;
        colormap gray;
        histogram2(SC_Gaus_FullRF_Gaussian_mean_mat(:,1),SC_Gaus_FullRF_Gaussian_mean_mat(:,2),...
            [10 10],...
            'XBinLimits',[0.5 (p.stim_columns+0.5)],...
            'YBinLimits',[0.5 (p.stim_columns+0.5)],...
            'Normalization','probability',...
            'DisplayStyle','tile',...
            'ShowEmptyBins','on'); hold on;
        xlabel('x');
        ylabel('y');
        axis equal;
        colorbar;
        title('SC -- Gaussian RF: Full Mean Dist.');
        set(gca,'Ydir','reverse');
        set(gca,'FontSize',12);
        set(gcf,'color','w');
        
        
        % Spectral RF - Kernel Density Smoothed
        KernSmth_xGrid = linspace(0.5,(p.stim_columns+0.5),41);
        KernSmth_yGrid = linspace(0.5,(p.stim_columns+0.5),41);
        KernSmth_xGrid_length = length(KernSmth_xGrid);
        KernSmth_yGrid_length = length(KernSmth_yGrid);
        Bwidth_Choice  = 2; % Option 1: default, Option 2: specified
        Bwidth_x       = 1; % 1, 2
        Bwidth_y       = 1; % 1, 2
        [x1,x2]        = meshgrid(KernSmth_xGrid, KernSmth_yGrid);
        x1             = x1(:);
        x2             = x2(:);
        xi             = [x1 x2];
        KernSmth_2D = cell(p.Spectral_Dim,1);
        for j = 1:p.Spectral_Dim
            if Bwidth_Choice == 1
                [KernSmth,~,bw] = ksdensity([SC_Gaus_Gaussian_mean_arr(:,j,1),SC_Gaus_Gaussian_mean_arr(:,j,2)],xi);
            else
                [KernSmth,~,bw] = ksdensity([SC_Gaus_Gaussian_mean_arr(:,j,1),SC_Gaus_Gaussian_mean_arr(:,j,2)],xi,'Bandwidth',[Bwidth_x Bwidth_y]);
            end
            KernSmth_2D{j} = reshape(KernSmth,[KernSmth_xGrid_length,KernSmth_yGrid_length]);
        end
        figure;
        if Hist_Colour_Choice == 1
            colormap gray;
        end
        for j = 1:p.Spectral_Dim
            ax = subplot(2,2,j);
            imagesc(KernSmth_xGrid,KernSmth_yGrid,KernSmth_2D{j});
            %set(gca,'YDir','normal'); % don't reverse, correct as is
            if Hist_Colour_Choice == 2
                colormap(ax,colorMap_arr(:,:,j));
            end
            if j>2
                xlabel('x');
            end
            if j==1||j==3
                ylabel('y');
            end
            axis equal;
            xlim([0.5 (p.stim_columns+0.5)]);
            ylim([0.5 (p.stim_columns+0.5)]);
            colorbar;
            title(Spectral_Names{j});
            set(gca,'FontSize',12);
        end
        annotation('textbox',[.5 .975 0 0],'String','SC -- Gaussian RF: Mean Dist.','FitBoxToText','on','FontSize',12,'FontWeight','bold','LineStyle','none','HorizontalAlignment','center','VerticalAlignment','middle');
        set(gcf,'color','w');
        
        
        % Full RF - Kernel Density Smoothed
        KernSmth_xGrid = linspace(0.5,(p.stim_columns+0.5),41);
        KernSmth_yGrid = linspace(0.5,(p.stim_columns+0.5),41);
        KernSmth_xGrid_length = length(KernSmth_xGrid);
        KernSmth_yGrid_length = length(KernSmth_yGrid);
        Bwidth_Choice  = 2; % Option 1: default, Option 2: specified
        Bwidth_x       = 1; % 1, 2
        Bwidth_y       = 1; % 1, 2
        [x1,x2]        = meshgrid(KernSmth_xGrid, KernSmth_yGrid);
        x1             = x1(:);
        x2             = x2(:);
        xi             = [x1 x2];
        if Bwidth_Choice == 1
            [KernSmth,~,bw] = ksdensity([SC_Gaus_FullRF_Gaussian_mean_mat(:,1),SC_Gaus_FullRF_Gaussian_mean_mat(:,2)],xi);
        else
            [KernSmth,~,bw] = ksdensity([SC_Gaus_FullRF_Gaussian_mean_mat(:,1),SC_Gaus_FullRF_Gaussian_mean_mat(:,2)],xi,'Bandwidth',[Bwidth_x Bwidth_y]);
        end
        KernSmth_2D = reshape(KernSmth,[KernSmth_xGrid_length,KernSmth_yGrid_length]);
        figure;
        colormap gray;
        imagesc(KernSmth_xGrid,KernSmth_yGrid,KernSmth_2D);
        %set(gca,'YDir','normal'); % don't reverse, correct as is
        xlabel('x');
        ylabel('y');
        axis equal;
        xlim([0.5 (p.stim_columns+0.5)]);
        ylim([0.5 (p.stim_columns+0.5)]);
        colorbar;
        title('SC -- Gaussian RF: Full Mean Dist.');
        set(gca,'FontSize',12);
        set(gcf,'color','w');
        
    end
    
end

autoArrangeFigures;

%Ask if next dataset shall be plotted (To avoid that all figures for all datasets
%are plotted all together)

end


out = 1;
end
% figure;
% for i = 1:p.Spectral_Dim
%
%     subplot(2,2,i);
%     hist_loop = histogram(Num_SC_ASP_RF_Pixels_mat(:,i));
%     if Hist_Colour_Choice == 2
%         set(hist_loop,'FaceColor',colorMap_arr(end,:,i));i
%     end
%     if i>2
%         xlabel('num. pixels');
%     end
%     if i==1||i==3
%         ylabel('num. cells');
%     end
%     title(Spectral_Names{i});
%     set(gca,'FontSize',12);
%     xlim([hist_loop.BinEdges(1) hist_loop.BinEdges(end)]);
%     ylim([0 max(hist_loop.Values)]);
%
% end
% annotation('textbox',[.5 .975 0 0],'String','STA-SD ASP -- Num RF Pixels','FitBoxToText','on','FontSize',12,'FontWeight','bold','LineStyle','none','HorizontalAlignment','center','VerticalAlignment','middle');
% set(gcf,'color','w');

