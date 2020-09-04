function out = spiketime_average_FFF_app (savepath, add_info)
%Load stimulus trace
nr_colours = 4;
stim_idx = add_info.stim_idx;

 %% Load saved bins
 S = load(findfile_app(stim_idx,savepath,'Bined_spikes.mat'));
 Bined_spikes = S.Bined_spikes;
 Bins_info = Bined_spikes(1).bins_info;
 
 S = load(savepath,'cell_indices');
 cell_indices = S.cell_indices;
 
 S = load(findfile_app(stim_idx,savepath,'Raster_spikes'));
 Raster_spikes = S.Raster_spikes;
 clear S
 
 max_nr_bins = max([Bins_info.nr_bins]);
 std_nr_bins = std([Bins_info.nr_bins]);
 
 if std_nr_bins > 2
     warning...
         ("Time intervals of stimulus have a high standard devitation, check trigger channel maybe?");
 end
 


 ep_idx = (1:1:nr_colours);
 nr_repeats = ceil(length(Bins_info)/nr_colours);
 ep_idx = repmat(ep_idx,1,nr_repeats); %Get indices for all episodes relative to nr of repeats
 

 nr_ep = length(Bins_info);
 ep_idx = ep_idx(1:nr_ep); %Cut of overlapping episodes
 
 
 
 nr_cells = length(Bined_spikes);
 

 %Average over episodes
for kk = 1:nr_cells
     for ii = 1:nr_colours
          ep_logical = ep_idx == ii;
          nr_repeats = size(Raster_spikes(1).spikes(:,ep_logical),2);
       try
           %Extract bined traces of specific episode for given cell
           cell_bins = Bined_spikes(kk).histcounts(:,ep_logical);
           FFF_average(kk).traces(:,ii) = nanmean(cell_bins,2);
           FFF_average(kk).cell_idx = cell_indices(kk);
           
           %Sort spikes for raster plot 
           FFF_average(kk).spikes(:,1:nr_repeats,ii) = Raster_spikes(kk).spikes(:,ep_logical)...
               +Bins_info(ii).hist_bins(1)-add_info.stim_begin;
           
           
           %Calculate quality criteria for responses
           var_of_mean = nanvar(FFF_average(kk).traces(:,ii));
           mean_of_var = nanmean(nanvar(cell_bins,[],1));
           FFF_average(kk).stats(:,ii) = var_of_mean/mean_of_var;
           FFF_average(kk).stats_mean = nanmean(FFF_average(kk).stats(:,ii));
            
       catch
           continue
       end
        
      end
       
  
    
    
end
 
% Clear empty entries in the structure
a = 1;
for ii = 1:nr_cells
    try
    if isempty(FFF_average(a).cell_idx)
        FFF_average(a) = [];
    else
        a = a+1;
    end
    catch
        continue
    end
end
 nr_cells_new = length(FFF_average);
%Save averaged traces
out = sf_organizer(stim_idx,savepath,'variable_name','FFF_average','variable',FFF_average);
 
 

 
%% Create stat figures
overview_FFF = figure;
bar([FFF_average.cell_idx],[FFF_average.stats_mean],'k');
title("Quality of FFF responses")
ylabel("Quality Index")
xlabel("Cell index")
hline(0.3,'-',"Threshold")
sf_organizer(add_info.stim_idx,savepath,'variable_name','overview_FFF',...
    'variable',overview_FFF);

%Collect quality criteria data
FFF_stats = NaN(nr_colours,length(FFF_average));
for ii = 1:nr_cells_new
    try
        FFF_stats(:,ii) = FFF_average(ii).stats;
    catch
        continue
    end
    
end

detail_stats_FFF = heatbar([FFF_average.cell_idx],FFF_stats);

ynames = {'UV','Blue','Green','Red'};
set(gca,'ytick',[1,2,3,4],'yticklabel',ynames);
title("Quality of responses by colour")
xlabel("Cell Index")

%Change the label of the colourbar (which is basically the ylabel in this
%case)
c = findall(detail_stats_FFF.fig.Children,'type','ColorBar');
c.Label.String = 'Quality index';
sf_organizer(add_info.stim_idx,savepath,'variable_name','detail_stats_FFF',...
    'variable',detail_stats_FFF.fig);




%Plot all traces as heatbar
total_bins = int64(max_nr_bins*nr_colours);
FFF_traces = NaN(nr_cells_new,total_bins);
for ii = 1:nr_cells_new
   try
   FFF_trace = FFF_average(ii).traces;
   FFF_trace = FFF_trace/nanmax(FFF_trace,[],'all');
   FFF_traces(ii,:) = reshape(FFF_trace,[1,total_bins]);
   catch
       continue
   end
end
binsize = Bined_spikes(1).bins_info.binsize;
stim_duration = max_nr_bins*binsize;
x_values = (binsize:binsize:stim_duration*nr_colours);
all_traces_plot = heatbar(x_values,FFF_traces,'Y',[FFF_average.cell_idx],...
    'gap', false);

title("All traces heatmap")
xlabel("Time in s")
ylabel("Cell Index")
c = findall(all_traces_plot.fig.Children,'type','ColorBar');
c.Label.String = 'Spikes per bin';
sf_organizer(add_info.stim_idx,savepath,'variable_name','all_traces_plot',...
    'variable',all_traces_plot.fig);

 
 
 
 
end