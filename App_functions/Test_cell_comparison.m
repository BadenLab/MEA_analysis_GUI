%Create add_info structure for draw gratings function




%%
All_Data = cell(height(plot_table_sorted),1);
for ii = 1:height(plot_table_sorted)
    %extract the indices
    Indices = plot_table_sorted.Indices(ii);
    Indices = Indices{:}';
    %%

    parfor kk = 1:length(Indices)
        
       
       temp = findcell(savepath,Indices(kk),'stim_idx',stim_idx); 
       
       Data(kk) = temp.Stim6;
     
    end
    All_Data{ii,1} = Data;
    
end

%%
%Only the first 10 clusters for now
close all
for ii = 1:10 %height(plot_table_sorted)
  
    Indices = plot_table_sorted.Indices(ii);
    Indices = Indices{1};
    QC_overview = zeros(length(Indices),1);
    
    for kk = 1:length(Indices)
        Data = All_Data{ii};
        QC_overview(kk) = Data(kk).Grating_QC.QC;
        Skewness(kk) = Data(kk).Data_circular
        idx = QC_overview > 0.5;
    end
    a(ii) = figure
    histogram(QC_overview,20,'Facecolor','k')
    title(['Cluster ',num2str(ii)])
    %call plot gratings function from here
    
    plot_fig(ii) = figure;
    plot_fig_raster(ii) = figure;
    add_info.stim_idx = stim_idx;
    add_info.Grating_info.panel_id = plot_fig(ii);
    add_info.Grating_info.panel_id2 = plot_fig_raster(ii);

    
    add_info.Grating_info.Cell = Indices(idx);
    
    out = draw_gratings(savepath,add_info);
    sgtitle(plot_fig(ii),['Cluster ',num2str(ii)]);
    saveas(plot_fig(ii),['Radial_groups',num2str(ii)],'svg')
    saveas(a(ii),['Radial_QC',num2str(ii)],'svg')
end
    %% Compare FFF
FFF_Data = cell(height(plot_table_sorted),1);
for ii = 1:height(plot_table_sorted)
    %extract the indices
    Indices = plot_table_sorted.Indices(ii);
    Indices = Indices{:}';
    %%
    Data = [];
    parfor kk = 1:length(Indices)
        
       
       temp = findcell(savepath,Indices(kk),'stim_idx',FFF_idx); 
       
       Data(kk) = temp.Stim7;
     
    end
    FFF_Data{ii,1} = Data;
    
end

%%

binsize = FFF_Data{1}(1).Binned_data.binsize;
max_bin_nr = length(FFF_Data{1}(1).FFF_average.traces);
max_bin_time = max_bin_nr * binsize;
xtrace = (0:binsize:nr_color*max_bin_time-binsize);

bins_step = max_bin_nr/2;
time_step = bins_step*binsize;

stim_events = (0:time_step:max(xtrace));
color_string_original ={'m','b','g','r'};
color_string = cell(1,length(stim_events));
a = 1;
for cc = 1:length(stim_events)
    
    if mod(cc,2) == 0
        color_string{cc} = 'k';
    else
        color_string{cc} = color_string_original{a};
        a = a+1;
    end
    
end

figure
plot_nr = 10; %height(plot_table_sorted);
for ii = 1:plot_nr
    
   Data_temp = FFF_Data{ii};
   
   trace_length = length(Data_temp(1).FFF_average.traces(:));
   FFF_traces = zeros(length(Data_temp),trace_length);
   
   for kk = 1:length(Data_temp)
       
       FFF_traces(kk,:) = Data_temp(kk).FFF_average.traces(:);
       
       
   end
   
    FFF_mean_traces(ii,:) = mean(FFF_traces,1);
    
    
    fx(ii) = subplot(plot_nr+1,1,ii);
    plot(xtrace,FFF_mean_traces(ii,:),'k')
    
end
fx(ii+1) = subplot(plot_nr+1,1,ii+1);
for kk = 1:length(stim_events)
    
rectangle('Position',[stim_events(kk),0 time_step 1],'FaceColor',color_string{kk});
hold on
end

linkaxes(fx)


