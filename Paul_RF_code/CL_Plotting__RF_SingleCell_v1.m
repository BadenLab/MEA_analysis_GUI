%% Plot RF and Associated Summary Statistics - Single Cell

% Choose which cell RF to plot
Cell_to_Plot  = 2;

RF_Ident_Plot = RF_Ident{Cell_to_Plot,1};C



if p.Plot_Choice == 1
    
    %%% RF Plots
    
    if p.RF_Ident_Meth_vec(1) == 1 % STA-SD method;
        
        if p.RF_Type(1) == 1     % Box
            
            fig1 = figure;
            if p.Spectral_Dim == 1
                colormap(fig1,gray(256));
%                 imagesc(RF_Ident_Plot.STA_SD); hold on;
%                 scatter(RF_Ident_Plot.RF_coords_centre(2),RF_Ident_Plot.RF_coords_centre(1),'rx','LineWidth',1.5);
%                 scatter(RF_Ident_Plot.RF_coords_noncentre(:,2),RF_Ident_Plot.RF_coords_noncentre(:,1),'go','LineWidth',1.5);
                axis equal; axis tight;
                xlabel('x');
                ylabel('y');
                title('RF');
                colorbar;
                set(gca,'FontSize',12);
            elseif p.Spectral_Dim == 4
                if Heat_Map_Colour_Choice == 1
                    colormap(fig1,gray(256));
                end
                for i = 1:p.Spectral_Dim
                    loop_plot = subplot(2,2,i);
                    imagesc(RF_Ident_Plot.STA_SD(:,:,i)); hold on;
                    if Heat_Map_Colour_Choice == 2
                        colormap(loop_plot,colorMap_arr(:,:,i));
                    end
                    if ~isempty(RF_Ident_Plot.STASD_Box_Num_RF_pixels{i})
                        scatter(RF_Ident_Plot.STASD_Box_RF_coords_centre{i}(2),RF_Ident_Plot.STASD_Box_RF_coords_centre{i}(1),'rx','LineWidth',1.5);
                        scatter(RF_Ident_Plot.STASD_Box_RF_coords_noncentre{i}(:,2),RF_Ident_Plot.STASD_Box_RF_coords_noncentre{i}(:,1),'go','LineWidth',1.5);
                    end
                    axis equal; axis tight;
                    if i>2
                        xlabel('x');
                    end
                    if i==1||i==3
                        ylabel('y');
                    end
                    title(Spectral_Names{i});
                    colorbar;
                    set(gca,'FontSize',12);
                end
            end
            annotation('textbox',[.5 .975 0 0],'String','STA-SD -- Box RF','FitBoxToText','on','FontSize',12,'FontWeight','bold','LineStyle','none','HorizontalAlignment','center','VerticalAlignment','middle');
            set(gcf,'color','w');
            
        end
        
        if p.RF_Type(2) == 1 % All Significant Pixels
            
            fig2 = figure;
            if p.Spectral_Dim == 1
                colormap(fig2,gray(256));
%                 imagesc(RF_Ident_Plot.STA_SD); hold on;
%                 scatter(RF_Ident_Plot.RF_coords(:,2),RF_Ident_Plot.RF_coords(:,1),'go','LineWidth',1.5);
                axis equal; axis tight;
                xlabel('x');
                ylabel('y');
                title('RF');
                colorbar;
                set(gca,'FontSize',12);
            elseif p.Spectral_Dim == 4
                if Heat_Map_Colour_Choice == 1
                    colormap(fig2,gray(256));
                end
                for i = 1:p.Spectral_Dim
                    loop_plot = subplot(2,2,i);
                    imagesc(RF_Ident_Plot.STA_SD(:,:,i)); hold on;
                    if Heat_Map_Colour_Choice == 2
                        colormap(loop_plot,colorMap_arr(:,:,i));
                        if ~isempty(RF_Ident_Plot.STASD_ASP_Num_RF_pixels{i})
                            scatter(RF_Ident_Plot.STASD_ASP_RF_coords{i}(:,2),RF_Ident_Plot.STASD_ASP_RF_coords{i}(:,1),'wo','LineWidth',1.5);
                        end
                    else % Heat_Map_Colour_Choice == 1
                        if ~isempty(RF_Ident_Plot.STASD_ASP_Num_RF_pixels{i})
                            scatter(RF_Ident_Plot.STASD_ASP_RF_coords{i}(:,2),RF_Ident_Plot.STASD_ASP_RF_coords{i}(:,1),'go','LineWidth',1.5);
                        end
                    end
                    axis equal; axis tight;
                    if i>2
                        xlabel('x');
                    end
                    if i==1||i==3
                        ylabel('y');
                    end
                    title(Spectral_Names{i});
                    colorbar;
                    set(gca,'FontSize',12);
                end
            end
            annotation('textbox',[.5 .975 0 0],'String','STA-SD -- All Sig. Pix. RF','FitBoxToText','on','FontSize',12,'FontWeight','bold','LineStyle','none','HorizontalAlignment','center','VerticalAlignment','middle');
            set(gcf,'color','w');
            
        end
        
        if p.RF_Type(3) == 1 % Gaussian
            
            fig3 = figure;
            if p.Spectral_Dim == 1
                colormap(fig3,gray(256));
%                 imagesc(RF_Ident_Plot.STA_SD,[min(RF_Ident_Plot.STA_SD,[],'all') max(RF_Ident_Plot.STA_SD,[],'all')]); hold on;
%                 scatter(RF_Ident_Plot.STASD_Gaus_RF_coords(:,2),RF_Ident_Plot.STASD_Gaus_RF_coords(:,1),'go','LineWidth',1.5);
%                 plot(RF_Ident_Plot.STASD_Gaus_Gaussian_mean(1),RF_Ident_Plot.STASD_Gaus_Gaussian_mean(2),'rx','LineWidth',1.5); %
%                 fcontour(RF_Ident_Plot.STASD_Gaus_gmPDF,[0.5 (p.stim_columns+0.5) 0.5 (p.stim_rows+0.5)],'r','LineWidth',1.5,'LevelList',RF_Ident_Plot.Thresh_height); % 'LineColor',LightestBlue,'LineWidth',2
                axis equal; axis tight;
                xlabel('x');
                ylabel('y');
                title('RF');
                colorbar;
                set(gca,'FontSize',12);
            elseif p.Spectral_Dim == 4
                if Heat_Map_Colour_Choice == 1
                    colormap(fig3,gray(256));
                end
                for i = 1:p.Spectral_Dim
                    loop_plot = subplot(2,2,i);
                    imagesc(RF_Ident_Plot.STA_SD(:,:,i),[min(RF_Ident_Plot.STA_SD(:,:,i),[],'all') max(RF_Ident_Plot.STA_SD(:,:,i),[],'all')]); hold on;
                    if Heat_Map_Colour_Choice == 2
                        colormap(loop_plot,colorMap_arr(:,:,i));
                        if ~isempty(RF_Ident_Plot.STASD_Gaus_Num_RF_pixels{i})
                            scatter(RF_Ident_Plot.STASD_Gaus_RF_coords{i}(:,2),RF_Ident_Plot.STASD_Gaus_RF_coords{i}(:,1),'wo','LineWidth',1.5);
                            plot(RF_Ident_Plot.STASD_Gaus_Gaussian_mean{i}(1),RF_Ident_Plot.STASD_Gaus_Gaussian_mean{i}(2),'wx','LineWidth',1.5);
                            fcontour(RF_Ident_Plot.STASD_Gaus_gmPDF{i},[0.5 (p.stim_columns+0.5) 0.5 (p.stim_rows+0.5)],'w','LineWidth',1.5,'LevelList',RF_Ident_Plot.STASD_Gaus_Thresh_height{i});
                        end
                    else % Heat_Map_Colour_Choice == 1
                        if ~isempty(RF_Ident_Plot.STASD_Gaus_Num_RF_pixels{i})
                            scatter(RF_Ident_Plot.STASD_Gaus_RF_coords{i}(:,2),RF_Ident_Plot.STASD_Gaus_RF_coords{i}(:,1),'go','LineWidth',1.5);
                            plot(RF_Ident_Plot.STASD_Gaus_Gaussian_mean{i}(1),RF_Ident_Plot.STASD_Gaus_Gaussian_mean{i}(2),'rx','LineWidth',1.5);
                            fcontour(RF_Ident_Plot.STASD_Gaus_gmPDF{i},[0.5 (p.stim_columns+0.5) 0.5 (p.stim_rows+0.5)],'r','LineWidth',1.5,'LevelList',RF_Ident_Plot.STASD_Gaus_Thresh_height{i});
                        end
                    end
                    axis equal; axis tight;
                    if i>2
                        xlabel('x');
                    end
                    if i==1||i==3
                        ylabel('y');
                    end
                    title(Spectral_Names{i});
                    colorbar;
                    set(gca,'FontSize',12);
                end
            end
            annotation('textbox',[.5 .975 0 0],'String','STA-SD -- Gaussian RF','FitBoxToText','on','FontSize',12,'FontWeight','bold','LineStyle','none','HorizontalAlignment','center','VerticalAlignment','middle');
            set(gcf,'color','w');
            
        end
        
    end
    
    if p.RF_Ident_Meth_vec(2) == 1 % LC method (Local Covariance)
        
        if p.RF_Type(1) == 1     % Box
            
            fig4 = figure;
            colormap(fig4,gray(256));
            imagesc(RF_Ident_Plot.Pixel_covar); hold on;
            scatter(RF_Ident_Plot.RF_coords_centre(2),RF_Ident_Plot.RF_coords_centre(1),'rx','LineWidth',1.5);
            scatter(RF_Ident_Plot.RF_coords_noncentre(:,2),RF_Ident_Plot.RF_coords_noncentre(:,1),'go','LineWidth',1.5);
            axis equal; axis tight;
            xlabel('x');
            ylabel('y');
            title('RF');
            colorbar;
            set(gca,'FontSize',12);
            set(gcf,'color','w');
            
        end
        
        if p.RF_Type(2) == 1 % All Significant Pixels
            
            fig5 = figure;
            colormap(fig5,gray(256));
            imagesc(RF_Ident_Plot.Pixel_covar); hold on;
            scatter(RF_Ident_Plot.RF_coords(:,2),RF_Ident_Plot.RF_coords(:,1),'go','LineWidth',1.5);
            axis equal; axis tight;
            xlabel('x');
            ylabel('y');
            title('RF');
            colorbar;
            set(gca,'FontSize',12);
            set(gcf,'color','w');
            
        end
        
        if p.RF_Type(3) == 1 % Gaussian
            
            fig6 = figure;
            colormap(fig6,gray(256));
            imagesc(RF_Ident_Plot.Pixel_covar); hold on;
            scatter(RF_Ident_Plot.RF_coords(:,2),RF_Ident_Plot.RF_coords(:,1),'go','LineWidth',1.5);
            plot(RF_Ident_Plot.Gaussian_mean(1),RF_Ident_Plot.Gaussian_mean(2),'rx','LineWidth',1.5); %
            fcontour(RF_Ident_Plot.gmPDF,[0.5 (p.stim_columns+0.5) 0.5 (p.stim_rows+0.5)],'r','LineWidth',1.5,'LevelList',RF_Ident_Plot.Thresh_height); % 'LineColor',LightestBlue,'LineWidth',2
            axis equal; axis tight;
            xlabel('x');
            ylabel('y');
            title('RF');
            colorbar;
            set(gca,'FontSize',12);
            set(gcf,'color','w');
            
        end
        
    end
    
    if p.RF_Ident_Meth_vec(3) == 1 % MI method (Mutual Information)
        
        if p.RF_Type(1) == 1     % Box
            
            fig7 = figure;
            colormap(fig7,gray(256));
            %imagesc(RF_Ident_Plot.Pixel_covar); hold on; --> Update this line!!
            scatter(RF_Ident_Plot.RF_coords_centre(2),RF_Ident_Plot.RF_coords_centre(1),'rx','LineWidth',1.5);
            scatter(RF_Ident_Plot.RF_coords_noncentre(:,2),RF_Ident_Plot.RF_coords_noncentre(:,1),'go','LineWidth',1.5);
            axis equal; axis tight;
            xlabel('x');
            ylabel('y');
            title('RF');
            colorbar;
            set(gca,'FontSize',12);
            set(gcf,'color','w');
            
        end
        
        if p.RF_Type(2) == 1 % All Significant Pixels
            
            fig8 = figure;
            colormap(fig8,gray(256));
            %imagesc(RF_Ident_Plot.Pixel_covar); hold on; --> Update this line!!
            scatter(RF_Ident_Plot.RF_coords(:,2),RF_Ident_Plot.RF_coords(:,1),'go','LineWidth',1.5);
            axis equal; axis tight;
            xlabel('x');
            ylabel('y');
            title('RF');
            colorbar;
            set(gca,'FontSize',12);
            set(gcf,'color','w');
            
        end
        
        if p.RF_Type(3) == 1 % Gaussian
            
            fig9 = figure;
            colormap(fig9,gray(256));
            %imagesc(RF_Ident_Plot.Pixel_covar); hold on; --> Update this line!!
            scatter(RF_Ident_Plot.RF_coords(:,2),RF_Ident_Plot.RF_coords(:,1),'go','LineWidth',1.5);
            plot(RF_Ident_Plot.Gaussian_mean(1),RF_Ident_Plot.Gaussian_mean(2),'rx','LineWidth',1.5); %
            fcontour(RF_Ident_Plot.gmPDF,[0.5 (p.stim_columns+0.5) 0.5 (p.stim_rows+0.5)],'r','LineWidth',1.5,'LevelList',RF_Ident_Plot.Thresh_height); % 'LineColor',LightestBlue,'LineWidth',2
            axis equal; axis tight;
            xlabel('x');
            ylabel('y');
            title('RF');
            colorbar;
            set(gca,'FontSize',12);
            set(gcf,'color','w');
            
        end
        
    end
    
    
    if p.RF_Ident_Meth_vec(4) == 1 % SC method;
        
        if p.RF_Type(1) == 1     % Box
            
            fig10 = figure;
            if p.Spectral_Dim == 1
                colormap(fig10,gray(256));
                %                 imagesc(RF_Ident_Plot.SCA_Pixel_covar); hold on;
                %                 scatter(RF_Ident_Plot.RF_coords_centre(2),RF_Ident_Plot.RF_coords_centre(1),'rx','LineWidth',1.5);
                %                 scatter(RF_Ident_Plot.RF_coords_noncentre(:,2),RF_Ident_Plot.RF_coords_noncentre(:,1),'go','LineWidth',1.5);
                axis equal; axis tight;
                xlabel('x');
                ylabel('y');
                title('RF');
                colorbar;
                set(gca,'FontSize',12);
            elseif p.Spectral_Dim == 4
                if Heat_Map_Colour_Choice == 1
                    colormap(fig10,gray(256));
                end
                for i = 1:p.Spectral_Dim
                    loop_plot = subplot(2,2,i);
                    imagesc(RF_Ident_Plot.SCA_Pixel_covar(:,:,i)); hold on;
                    if Heat_Map_Colour_Choice == 2
                        colormap(loop_plot,colorMap_arr(:,:,i));
                        if ~isempty(RF_Ident_Plot.SC_Box_Num_RF_pixels{i})
                            scatter(RF_Ident_Plot.SC_Box_RF_coords_centre{i}(2),RF_Ident_Plot.SC_Box_RF_coords_centre{i}(1),'wx','LineWidth',1.5);
                            scatter(RF_Ident_Plot.SC_Box_RF_coords_noncentre{i}(:,2),RF_Ident_Plot.SC_Box_RF_coords_noncentre{i}(:,1),'wo','LineWidth',1.5);
                        end
                    else % Heat_Map_Colour_Choice == 1
                        if ~isempty(RF_Ident_Plot.SC_Box_Num_RF_pixels{i})
                            scatter(RF_Ident_Plot.SC_Box_RF_coords_centre{i}(2),RF_Ident_Plot.SC_Box_RF_coords_centre{i}(1),'rx','LineWidth',1.5);
                            scatter(RF_Ident_Plot.SC_Box_RF_coords_noncentre{i}(:,2),RF_Ident_Plot.SC_Box_RF_coords_noncentre{i}(:,1),'go','LineWidth',1.5);
                        end
                    end
                    axis equal; axis tight;
                    if i>2
                        xlabel('x');
                    end
                    if i==1||i==3
                        ylabel('y');
                    end
                    title(Spectral_Names{i});
                    colorbar;
                    set(gca,'FontSize',12);
                end
            end
            annotation('textbox',[.5 .975 0 0],'String','Self Covar. -- Box RF','FitBoxToText','on','FontSize',12,'FontWeight','bold','LineStyle','none','HorizontalAlignment','center','VerticalAlignment','middle');
            set(gcf,'color','w');
            
        end
        
        if p.RF_Type(2) == 1 % All Significant Pixels
            
            fig11 = figure;
            if p.Spectral_Dim == 1
                colormap(fig11,gray(256));
                %                 imagesc(RF_Ident_Plot.SCA_Pixel_covar); hold on;
                %                 scatter(RF_Ident_Plot.RF_coords(:,2),RF_Ident_Plot.RF_coords(:,1),'go','LineWidth',1.5);
                axis equal; axis tight;
                xlabel('x');
                ylabel('y');
                title('RF');
                colorbar;
                set(gca,'FontSize',12);
            elseif p.Spectral_Dim == 4
                if Heat_Map_Colour_Choice == 1
                    colormap(fig11,gray(256));
                end
                for i = 1:p.Spectral_Dim
                    loop_plot = subplot(2,2,i);
                    imagesc(RF_Ident_Plot.SCA_Pixel_covar(:,:,i)); hold on;
                    if Heat_Map_Colour_Choice == 2
                        colormap(loop_plot,colorMap_arr(:,:,i));
                        if ~isempty(RF_Ident_Plot.SC_ASP_Num_RF_pixels{i})
                            scatter(RF_Ident_Plot.SC_ASP_RF_coords{i}(:,2),RF_Ident_Plot.SC_ASP_RF_coords{i}(:,1),'wo','LineWidth',1.5);
                        end
                    else % Heat_Map_Colour_Choice == 1
                        if ~isempty(RF_Ident_Plot.SC_ASP_Num_RF_pixels{i})
                            scatter(RF_Ident_Plot.SC_ASP_RF_coords{i}(:,2),RF_Ident_Plot.SC_ASP_RF_coords{i}(:,1),'go','LineWidth',1.5);
                        end
                    end
                    axis equal; axis tight;
                    if i>2
                        xlabel('x');
                    end
                    if i==1||i==3
                        ylabel('y');
                    end
                    title(Spectral_Names{i});
                    colorbar;
                    set(gca,'FontSize',12);
                end
            end
            annotation('textbox',[.5 .975 0 0],'String','Self Covar. -- All Sig. Pix. RF','FitBoxToText','on','FontSize',12,'FontWeight','bold','LineStyle','none','HorizontalAlignment','center','VerticalAlignment','middle');
            set(gcf,'color','w');
            
        end
        
        if p.RF_Type(3) == 1 % Gaussian
            
            fig12 = figure;
            if p.Spectral_Dim == 1
                colormap(fig12,gray(256));
%                 imagesc(RF_Ident_Plot.SCA_Pixel_covar,[min(RF_Ident_Plot.SCA_Pixel_covar,[],'all') max(RF_Ident_Plot.SCA_Pixel_covar,[],'all')]); hold on;
%                 scatter(RF_Ident_Plot.SC_Gaus_RF_coords(:,2),RF_Ident_Plot.SC_Gaus_RF_coords(:,1),'go','LineWidth',1.5);
%                 plot(RF_Ident_Plot.SC_Gaus_Gaussian_mean(1),RF_Ident_Plot.SC_Gaus_Gaussian_mean(2),'rx','LineWidth',1.5); %
%                 fcontour(RF_Ident_Plot.SC_Gaus_gmPDF,[0.5 (p.stim_columns+0.5) 0.5 (p.stim_rows+0.5)],'r','LineWidth',1.5,'LevelList',RF_Ident_Plot.Thresh_height); % 'LineColor',LightestBlue,'LineWidth',2
                axis equal; axis tight;
                xlabel('x');
                ylabel('y');
                title('RF');
                colorbar;
                set(gca,'FontSize',12);
            elseif p.Spectral_Dim == 4
                if Heat_Map_Colour_Choice == 1
                    colormap(fig12,gray(256));
                end
                for i = 1:p.Spectral_Dim
                    loop_plot = subplot(2,2,i);
                    imagesc(RF_Ident_Plot.SCA_Pixel_covar(:,:,i),[min(RF_Ident_Plot.SCA_Pixel_covar(:,:,i),[],'all') max(RF_Ident_Plot.SCA_Pixel_covar(:,:,i),[],'all')]); hold on;
                    if Heat_Map_Colour_Choice == 2
                        colormap(loop_plot,colorMap_arr(:,:,i));
                        if ~isempty(RF_Ident_Plot.SC_Gaus_Num_RF_pixels{i})
                            scatter(RF_Ident_Plot.SC_Gaus_RF_coords{i}(:,2),RF_Ident_Plot.SC_Gaus_RF_coords{i}(:,1),'wo','LineWidth',1.5);
                            plot(RF_Ident_Plot.SC_Gaus_Gaussian_mean{i}(1),RF_Ident_Plot.SC_Gaus_Gaussian_mean{i}(2),'wx','LineWidth',1.5);
                            fcontour(RF_Ident_Plot.SC_Gaus_gmPDF{i},[0.5 (p.stim_columns+0.5) 0.5 (p.stim_rows+0.5)],'w','LineWidth',1.5,'LevelList',RF_Ident_Plot.SC_Gaus_Thresh_height{i});
                        end
                    else % Heat_Map_Colour_Choice == 1
                        if ~isempty(RF_Ident_Plot.SC_Gaus_Num_RF_pixels{i})
                            scatter(RF_Ident_Plot.SC_Gaus_RF_coords{i}(:,2),RF_Ident_Plot.SC_Gaus_RF_coords{i}(:,1),'go','LineWidth',1.5);
                            plot(RF_Ident_Plot.SC_Gaus_Gaussian_mean{i}(1),RF_Ident_Plot.SC_Gaus_Gaussian_mean{i}(2),'rx','LineWidth',1.5);
                            fcontour(RF_Ident_Plot.SC_Gaus_gmPDF{i},[0.5 (p.stim_columns+0.5) 0.5 (p.stim_rows+0.5)],'r','LineWidth',1.5,'LevelList',RF_Ident_Plot.SC_Gaus_Thresh_height{i});
                        end
                    end
                    axis equal; axis tight;
                    if i>2
                        xlabel('x');
                    end
                    if i==1||i==3
                        ylabel('y');
                    end
                    title(Spectral_Names{i});
                    colorbar;
                    set(gca,'FontSize',12);
                end
            end
            annotation('textbox',[.5 .975 0 0],'String','Self Covar. -- Gaussian RF','FitBoxToText','on','FontSize',12,'FontWeight','bold','LineStyle','none','HorizontalAlignment','center','VerticalAlignment','middle');
            set(gcf,'color','w');
            
        end
        
    end
    
    
    
    
    
    
    
    
    %%% RF Significance Plots
    
    if p.RF_Ident_Meth_vec(1) == 1 % STA-SD method;
        
        figure;
        for i = 1:p.Spectral_Dim
            subplot(2,2,i);
            hist_loop = histogram(RF_Ident_Plot.STA_SD(:,:,i)); hold on;
            if Hist_Colour_Choice == 2
                set(hist_loop,'FaceColor',colorMap_arr(end,:,i));
            end
            vline_loop_1 = vline(mean(RF_Ident_Plot.STA_SD(:,:,i),'all'),'r'); % mean
            vline_loop_2 = vline(RF_Ident_Plot.STA_SD_Sig_Thresh(i),'g'); % significance thrshold
            vline_loop_3 = vline(hist_loop.BinEdges(end),'b:');   % RH edge of rightmost bin
            if i>2
                xlabel('STA-SD.');
            end
            if i==1||i==3
                ylabel('num. pixels');
            end
            title(Spectral_Names{i});
            set(gca,'FontSize',12);
            set(vline_loop_1,'LineWidth',1.5);
            set(vline_loop_2,'LineWidth',1.5);
            set(vline_loop_3,'LineWidth',1.5);
            xlim([hist_loop.BinEdges(1) max(hist_loop.BinEdges(end),RF_Ident_Plot.STA_SD_Sig_Thresh(i))]);
            ylim([0 max(hist_loop.Values)]);
        end
        annotation('textbox',[.5 .975 0 0],'String','STA-SD','FitBoxToText','on','FontSize',12,'FontWeight','bold','LineStyle','none','HorizontalAlignment','center','VerticalAlignment','middle'); % dim: [x y w h] [.45 .7 .3 .3]
        set(gcf,'color','w');
        
    end
    
    if p.RF_Ident_Meth_vec(2) == 1 % LC method (Local Covariance)
        
        figure;
        for i = 1:p.Spectral_Dim
            subplot(2,2,i);
            hist_loop = histogram(RF_Ident_Plot.LCA_Stixel_covar(:,:,i)); hold on;
            if Hist_Colour_Choice == 2
                set(hist_loop,'FaceColor',colorMap_arr(end,:,i));
            end
            vline_loop_1 = vline(mean(RF_Ident_Plot.LCA_Stixel_covar(:,:,i),'all'),'r'); % mean
            vline_loop_2 = vline(RF_Ident_Plot.LC_Sig_Thresh(i),'g');   % significance thrshold
            vline_loop_3 = vline(hist_loop.BinEdges(end),'b:'); % RH edge of rightmost bin
            if i>2
                xlabel('STA-SD.');
            end
            if i==1||i==3
                ylabel('num. pixels');
            end
            title(Spectral_Names{i});
            set(gca,'FontSize',12);
            set(vline_loop_1,'LineWidth',1.5);
            set(vline_loop_2,'LineWidth',1.5);
            set(vline_loop_3,'LineWidth',1.5);
            xlim([hist_loop.BinEdges(1) max(hist_loop.BinEdges(end),RF_Ident_Plot.LC_Sig_Thresh(i))]);
            ylim([0 max(hist_loop.Values)]);
        end
        annotation('textbox',[.5 .975 0 0],'String','Local Covariance','FitBoxToText','on','FontSize',12,'FontWeight','bold','LineStyle','none','HorizontalAlignment','center','VerticalAlignment','middle');
        set(gcf,'color','w');
        
    end
    
    if p.RF_Ident_Meth_vec(3) == 1 % MI method (Mutual Information)
        
        figure;
        for i = 1:p.Spectral_Dim
            subplot(2,2,i);
            %hist_loop = histogram(); hold on;
            if Hist_Colour_Choice == 2
                set(hist_loop,'FaceColor',colorMap_arr(end,:,i));
            end
            %vline_loop_1 = vline(mean(,'all'),'r'); % mean
            vline_loop_2 = vline(RF_Ident_Plot.MI_Sig_Thresh(i),'g');   % significance thrshold
            vline_loop_3 = vline(hist_loop.BinEdges(end),'b:'); % RH edge of rightmost bin
            if i>2
                xlabel('STA-SD.');
            end
            if i==1||i==3
                ylabel('num. pixels');
            end
            title(Spectral_Names{i});
            set(gca,'FontSize',12);
            set(vline_loop_1,'LineWidth',1.5);
            set(vline_loop_2,'LineWidth',1.5);
            set(vline_loop_3,'LineWidth',1.5);
            xlim([hist_loop.BinEdges(1) max(hist_loop.BinEdges(end),RF_Ident_Plot.MI_Sig_Thresh(i))]);
            ylim([0 max(hist_loop.Values)]);
        end
        annotation('textbox',[.5 .975 0 0],'String','Mutual Information','FitBoxToText','on','FontSize',12,'FontWeight','bold','LineStyle','none','HorizontalAlignment','center','VerticalAlignment','middle');
        set(gcf,'color','w');
        
    end
    
    if p.RF_Ident_Meth_vec(4) == 1 % SC method;
        
        figure;
        for i = 1:p.Spectral_Dim
            subplot(2,2,i);
            hist_loop = histogram(RF_Ident_Plot.SCA_Pixel_covar(:,:,i)); hold on;
            if Hist_Colour_Choice == 2
                set(hist_loop,'FaceColor',colorMap_arr(end,:,i));
            end
            vline_loop_1 = vline(mean(RF_Ident_Plot.SCA_Pixel_covar(:,:,i),'all'),'r'); % mean
            vline_loop_2 = vline(RF_Ident_Plot.SC_Sig_Thresh(i),'g'); % significance thrshold
            vline_loop_3 = vline(hist_loop.BinEdges(1),'b:'); % LH edge of leftmost bin
            if i>2
                xlabel('self-covar.');
            end
            if i==1||i==3
                ylabel('num. pixels');
            end
            title(Spectral_Names{i});
            set(gca,'FontSize',12);
            set(vline_loop_1,'LineWidth',1.5);
            set(vline_loop_2,'LineWidth',1.5);
            set(vline_loop_3,'LineWidth',1.5);
            xlim([min(hist_loop.BinEdges(1),RF_Ident_Plot.SC_Sig_Thresh(i)) hist_loop.BinEdges(end)]);
            ylim([0 max(hist_loop.Values)]);
        end
        annotation('textbox',[.5 .975 0 0],'String','Self Covariance','FitBoxToText','on','FontSize',12,'FontWeight','bold','LineStyle','none','HorizontalAlignment','center','VerticalAlignment','middle');
        set(gcf,'color','w');
        
    end
    
    
    
    
     %%% Pure Heatmap Plots
     
     if p.RF_Ident_Meth_vec(1) == 1 % STA-SD method;
         
         fig13 = figure;
         if p.Spectral_Dim == 1
             colormap(fig13,gray(256));
             %                 imagesc(RF_Ident_Plot.STA_SD); hold on;
             axis equal; axis tight;
             xlabel('x');
             ylabel('y');
             title('RF');
             colorbar;
             set(gca,'FontSize',12);
         elseif p.Spectral_Dim == 4
             if Heat_Map_Colour_Choice == 1
                 colormap(fig13,gray(256));
             end
             for i = 1:p.Spectral_Dim
                 loop_plot = subplot(2,2,i);
                 imagesc(RF_Ident_Plot.STA_SD(:,:,i)); hold on;
                 if Heat_Map_Colour_Choice == 2
                     colormap(loop_plot,colorMap_arr(:,:,i));
                 end
                 axis equal; axis tight;
                 if i>2
                     xlabel('x');
                 end
                 if i==1||i==3
                     ylabel('y');
                 end
                 title(Spectral_Names{i});
                 colorbar;
                 set(gca,'FontSize',12);
             end
         end
         annotation('textbox',[.5 .975 0 0],'String','STA-SD','FitBoxToText','on','FontSize',12,'FontWeight','bold','LineStyle','none','HorizontalAlignment','center','VerticalAlignment','middle');
         set(gcf,'color','w');
         
     end
     
     if p.RF_Ident_Meth_vec(2) == 1 % LC method;
         
         fig14 = figure;
         if p.Spectral_Dim == 1
             colormap(fig14,gray(256));
             %                 imagesc(RF_Ident_Plot.Pixel_covar); hold on;
             axis equal; axis tight;
             xlabel('x');
             ylabel('y');
             title('RF');
             colorbar;
             set(gca,'FontSize',12);
         elseif p.Spectral_Dim == 4
             if Heat_Map_Colour_Choice == 1
                 colormap(fig14,gray(256));
             end
             for i = 1:p.Spectral_Dim
                 loop_plot = subplot(2,2,i);
                 imagesc(RF_Ident_Plot.Pixel_covar(:,:,i)); hold on;
                 if Heat_Map_Colour_Choice == 2
                     colormap(loop_plot,colorMap_arr(:,:,i));
                 end
                 axis equal; axis tight;
                 if i>2
                     xlabel('x');
                 end
                 if i==1||i==3
                     ylabel('y');
                 end
                 title(Spectral_Names{i});
                 colorbar;
                 set(gca,'FontSize',12);
             end
         end
         annotation('textbox',[.5 .975 0 0],'String','Local Covariance','FitBoxToText','on','FontSize',12,'FontWeight','bold','LineStyle','none','HorizontalAlignment','center','VerticalAlignment','middle');
         set(gcf,'color','w');
         
     end
     
     if p.RF_Ident_Meth_vec(3) == 1 % MI method;
         
         fig15 = figure;
         if p.Spectral_Dim == 1
             colormap(fig15,gray(256));
             %                 imagesc(RF_Ident_Plot.Pixel_covar); hold on;
             axis equal; axis tight;
             xlabel('x');
             ylabel('y');
             title('RF');
             colorbar;
             set(gca,'FontSize',12);
         elseif p.Spectral_Dim == 4
             if Heat_Map_Colour_Choice == 1
                 colormap(fig15,gray(256));
             end
             for i = 1:p.Spectral_Dim
                 loop_plot = subplot(2,2,i);
                 %imagesc(RF_Ident_Plot.Pixel_covar(:,:,i)); hold on; Update!!!
                 if Heat_Map_Colour_Choice == 2
                     colormap(loop_plot,colorMap_arr(:,:,i));
                 end
                 axis equal; axis tight;
                 if i>2
                     xlabel('x');
                 end
                 if i==1||i==3
                     ylabel('y');
                 end
                 title(Spectral_Names{i});
                 colorbar;
                 set(gca,'FontSize',12);
             end
         end
         annotation('textbox',[.5 .975 0 0],'String','Mutual Information','FitBoxToText','on','FontSize',12,'FontWeight','bold','LineStyle','none','HorizontalAlignment','center','VerticalAlignment','middle');
         set(gcf,'color','w');
         
     end
     
     if p.RF_Ident_Meth_vec(4) == 1 % SC method;
         
         fig16 = figure;
         if p.Spectral_Dim == 1
             colormap(fig16,gray(256));
             %                 imagesc(RF_Ident_Plot.SCA_Pixel_covar); hold on;
             axis equal; axis tight;
             xlabel('x');
             ylabel('y');
             title('RF');
             colorbar;
             set(gca,'FontSize',12);
         elseif p.Spectral_Dim == 4
             if Heat_Map_Colour_Choice == 1
                 colormap(fig16,gray(256));
             end
             for i = 1:p.Spectral_Dim
                 loop_plot = subplot(2,2,i);
                 imagesc(RF_Ident_Plot.SCA_Pixel_covar(:,:,i)); hold on;
                 if Heat_Map_Colour_Choice == 2
                     colormap(loop_plot,colorMap_arr(:,:,i));
                 end
                 axis equal; axis tight;
                 if i>2
                     xlabel('x');
                 end
                 if i==1||i==3
                     ylabel('y');
                 end
                 title(Spectral_Names{i});
                 colorbar;
                 set(gca,'FontSize',12);
             end
         end
         annotation('textbox',[.5 .975 0 0],'String','Self Covariance','FitBoxToText','on','FontSize',12,'FontWeight','bold','LineStyle','none','HorizontalAlignment','center','VerticalAlignment','middle');
         set(gcf,'color','w');
         
     end
    
    
end

%% Save data

%save data_MChick_1_Cell_7_v_1;






% figure;
% colormap(gray(256));
% subplot(1,5,1);
% imagesc(RF_Ident_Plot.Stixel_covar(:,:,1)); axis equal; axis tight;
% subplot(1,5,2);
% imagesc(RF_Ident_Plot.Stixel_covar(:,:,2)); axis equal; axis tight;
% subplot(1,5,3);
% imagesc(RF_Ident_Plot.Stixel_covar(:,:,3)); axis equal; axis tight;
% subplot(1,5,4);
% imagesc(RF_Ident_Plot.Stixel_covar(:,:,4)); axis equal; axis tight;
% subplot(1,5,5);
% imagesc(RF_Ident_Plot.Stixel_covar(:,:,5)); axis equal; axis tight;
% set(gcf,'color','w');     
%      
%      
%      
% fig2 = figure;
% colormap(fig2,gray(256));
% xslice = [];   
% yslice = [];
% zslice = [1,2,3,4,5];
% slice(RF_Ident_Plot.Stixel_covar,xslice,yslice,zslice);
