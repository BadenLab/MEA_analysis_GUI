%% STA & SCA Plotting Code - 1
% 26,01,2021 Onwards

function Plot_STA_app (RF_Ident_Plot,settings,RF_type)
%% Plotting Choices

% Choose which cell RF to plot
%Cell_to_Plot  = 19; % 7
%RF_Ident_Plot = RF_Ident{Cell_to_Plot,1};

% Choose whether to use B&W or coloured heatmaps
% 1 = B&W;
% 2 = coloured.
%Heat_Map_Colour_Choice = 2;

% Allow for 2 colour cases: monochromatic and tetrachromatic (RGBUV)



if settings.kernel_new.plot_all_STA
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

if RF_type == 1 % STA-SD method;
    Spectral_Dim = size(RF_Ident_Plot.STA_SD,3);
    figSTA = figure;
    if Spectral_Dim == 1
    elseif Spectral_Dim == 4
        if settings.kernel_new.Heat_gray == 1
            colormap(figSTA,gray(256));
        end
        for i = 1:Spectral_Dim
            for j = First_Bin:Last_Bin%1:p.Num_STE_bins
                loop_plot = subplot_tight(Spectral_Dim,Num_Selected_Bins,Num_Selected_Bins*(i-1)+j);%loop_plot = subplot(p.Spectral_Dim,p.Num_STE_bins,p.Num_STE_bins*(i-1)+j);
                imagesc(RF_Ident_Plot.STA(:,:,j,i)); hold on;
                if settings.kernel_new.Heat_colour == 2
                    colormap(loop_plot,colorMap_arr(:,:,i));
                    if settings.kernel_new.plot_Allpixel == 1 % All Significant Pixels
                        if ~isempty(RF_Ident_Plot.STASD_ASP_Num_RF_pixels{i})
                            scatter(RF_Ident_Plot.STASD_ASP_RF_coords{i}(:,2),RF_Ident_Plot.STASD_ASP_RF_coords{i}(:,1),'wo','LineWidth',1.5);
                        end
                    end
                else % Heat_Map_Colour_Choice == 1
                    if settings.kernel_new.plot_Allpixel == 1 % All Significant Pixels
                        if ~isempty(RF_Ident_Plot.STASD_ASP_Num_RF_pixels{i})
                            scatter(RF_Ident_Plot.STASD_ASP_RF_coords{i}(:,2),RF_Ident_Plot.STASD_ASP_RF_coords{i}(:,1),'go','LineWidth',1.5);
                        end
                    end
                end
                axis equal; axis tight;
                if i==Spectral_Dim
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



% SCA Full

if RF_type == 2 % SCA method;
    Spectral_Dim = size(RF_Ident_Plot.STA_SD,3);
    figSCA = figure;
    if Spectral_Dim == 1
    elseif Spectral_Dim == 4
        if settings.kernel_new.Heat_gray
            colormap(figSCA,gray(256));
        end
        for i = 1:Spectral_Dim
            for j = First_Bin:Last_Bin%1:p.Num_STE_bins
                loop_plot = subplot(Spectral_Dim,Num_Selected_Bins,Num_Selected_Bins*(i-1)+j); %subplot(p.Spectral_Dim,p.Num_STE_bins,p.Num_STE_bins*(i-1)+j);
                imagesc(RF_Ident_Plot.SCA_Stixel_covar(:,:,j,i)); hold on;
                if settings.kernel_new.Heat_colour
                    colormap(loop_plot,colorMap_arr(:,:,i));
                    if settings.kernel_new.plot_Allpixel % All Significant Pixels
                        if ~isempty(RF_Ident_Plot.SC_ASP_Num_RF_pixels{i})
                            scatter(RF_Ident_Plot.SC_ASP_RF_coords{i}(:,2),RF_Ident_Plot.SC_ASP_RF_coords{i}(:,1),'wo','LineWidth',1.5);
                        end
                    end
                else % Heat_Map_Colour_Choice == 1
                    if settings.kernel_new.plot_Allpixel % All Significant Pixels
                        if ~isempty(RF_Ident_Plot.SC_ASP_Num_RF_pixels{i})
                            scatter(RF_Ident_Plot.SC_ASP_RF_coords{i}(:,2),RF_Ident_Plot.SC_ASP_RF_coords{i}(:,1),'go','LineWidth',1.5);
                        end
                    end
                end
                axis equal; axis tight;
                if i==Spectral_Dim
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
if settings.kernel_new.plot_STA_sig_pix  
colour_string = {'r','g','b','m'};  
if RF_type == 1
     Spectral_Dim = size(RF_Ident_Plot.STA_SD,3);
    
 STA = RF_Ident_Plot.STA;
 % get STA for significant pixel
 if settings.kernel_new.plot_Allpixel
    figure
     for ii = 1:RF_Ident_Plot.STASD_ASP_FullRF_Num_pixels
         RF_coords = RF_Ident_Plot.STASD_ASP_FullRF_coords;
         
         traces = squeeze(STA(RF_coords(ii,1),RF_coords(ii,2),:,:));
         
         % Plotting
         
         nr_subplots = numSubplots(RF_Ident_Plot.STASD_ASP_FullRF_Num_pixels);
         
         ax(ii) = subplot(nr_subplots(1),nr_subplots(2),ii);
         for kk = 1:Spectral_Dim
             hold on
             plot(traces(:,kk),'-','Color',colour_string{kk})
             hold off
         end
             
            
     end
 end
     
 
 
  if settings.kernel_new.plot_Box	
    figure
     for ii = 1:RF_Ident_Plot.STASD_Box_FullRF_Num_pixels	
         RF_coords = RF_Ident_Plot.STASD_Box_FullRF_coords;
         
         traces = squeeze(STA(RF_coords(ii,1),RF_coords(ii,2),:,:));
         
         % Plotting
         
         nr_subplots = numSubplots(RF_Ident_Plot.STASD_Box_FullRF_Num_pixels);
         
         ax(ii) = subplot(nr_subplots(1),nr_subplots(2),ii);
         for kk = 1:Spectral_Dim
             hold on
             plot(traces(:,kk),'-','Color',colour_string{kk})
             hold off
         end
             
            
     end
 end
     
if settings.kernel_new.plot_Gauss     
    figure
     for ii = 1:RF_Ident_Plot.STASD_Gaus_FullRF_Num_pixels	
         RF_coords = RF_Ident_Plot.STASD_Gaus_FullRF_coords;
         
         traces = squeeze(STA(RF_coords(ii,1),RF_coords(ii,2),:,:));
         
         % Plotting
         
         nr_subplots = numSubplots(RF_Ident_Plot.STASD_Gaus_FullRF_Num_pixels);
         
         ax(ii) = subplot(nr_subplots(1),nr_subplots(2),ii);
         for kk = 1:Spectral_Dim
             hold on
             plot(traces(:,kk),'-','Color',colour_string{kk})
             hold off
         end
             
            
     end
 end 
     
     
end
 


if RF_type == 2
    
%  STA = RF_Ident_Plot.STA;
%  % get STA for significant pixel
%  if settings.kernel_new.plot_Allpixel
%     figure
%      for ii = 1:RF_Ident_Plot.SC_ASP_FullRF_Num_pixels
%          RF_coords = RF_Ident_Plot.SC_ASP_FullRF_coords;
%          
%          traces = squeeze(STA(RF_coords(ii,1),RF_coords(ii,2),:,:));
%          
%          % Plotting
%          
%          nr_subplots = numSubplots(RF_Ident_Plot.SC_ASP_FullRF_Num_pixels);
%          
%          ax(ii) = subplot(nr_subplots(1),nr_subplots(2),ii);
%          for kk = 1:Spectral_Dim
%              hold on
%              plot(traces(:,kk),'-','Color',colour_string{kk})
%              hold off
%          end
%              
%             
%      end
%  end
%      
%  
%  
%   if settings.kernel_new.plot_Box	
%     figure
%      for ii = 1:RF_Ident_Plot.SC_Box_FullRF_Num_pixels	
%          RF_coords = RF_Ident_Plot.SC_Box_FullRF_coords;
%          
%          traces = squeeze(STA(RF_coords(ii,1),RF_coords(ii,2),:,:));
%          
%          % Plotting
%          
%          nr_subplots = numSubplots(RF_Ident_Plot.SC_Box_FullRF_Num_pixels);
%          
%          ax(ii) = subplot(nr_subplots(1),nr_subplots(2),ii);
%          for kk = 1:Spectral_Dim
%              hold on
%              plot(traces(:,kk),'-','Color',colour_string{kk})
%              hold off
%          end
%              
%             
%      end
%  end
%      
% if settings.kernel_new.plot_Gauss     
%     figure
%      for ii = 1:RF_Ident_Plot.SC_Gaus_FullRF_Num_pixels	
%          RF_coords = RF_Ident_Plot.SC_Gaus_FullRF_coords;
%          
%          traces = squeeze(STA(RF_coords(ii,1),RF_coords(ii,2),:,:));
%          
%          % Plotting
%          
%          nr_subplots = numSubplots(RF_Ident_Plot.SC_Gaus_FullRF_Num_pixels);
%          
%          ax(ii) = subplot(nr_subplots(1),nr_subplots(2),ii);
%          for kk = 1:Spectral_Dim
%              hold on
%              plot(traces(:,kk),'-','Color',colour_string{kk})
%              hold off
%          end
%              
%             
%      end
%  end 
%      
%      
 end
%  
    
    
    
    
end



