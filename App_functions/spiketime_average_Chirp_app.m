function out = spiketime_average_Chirp_app (savepath, add_info)
%Load stimulus trace
nr_colours = 2;
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
 
 %Find the two different nr of bins
 nr_bins_unique = unique([Bins_info.nr_bins]);
 
 
 ep_idx = (1:1:nr_colours);
 nr_repeats = ceil(length(Bins_info)/nr_colours);
 ep_idx = repmat(ep_idx,1,nr_repeats); %Get indices for all episodes relative to nr of repeats
 

 nr_ep = length(Bins_info);
 ep_idx = ep_idx(1:nr_ep); %Cut of overlapping episodes
  
 nr_cells = length(Bined_spikes);
 
 %Preallocate
 for ii = 1:nr_cells
     Chirp_average(ii).traces = zeros(int64(nr_bins_unique(2)),nr_colours);
     
 end
 
 
 %Average over episodes
for kk = 1:nr_cells
     for ii = 1:nr_colours
          ep_logical = ep_idx == ii;
       
       try
           nr_repeats = size(Raster_spikes(1).spikes(:,ep_logical),2);
           %Extract bined traces of specific episode for given cell
           cell_bins = Bined_spikes(kk).histcounts(1:nr_bins_unique(ii),ep_logical);
           Chirp_average(kk).traces(1:nr_bins_unique(ii),ii) = nanmean(cell_bins,2);
           Chirp_average(kk).cell_idx = cell_indices(kk);
           
           %Sort spikes for raster plot 
           Chirp_average(kk).spikes(:,1:nr_repeats,ii) = Raster_spikes(kk).spikes(:,ep_logical)...
               +Bins_info(ii).hist_bins(1)-add_info.stim_begin;
           
           
           %Calculate quality criteria for responses
           var_of_mean = nanvar(Chirp_average(kk).traces(:,ii));
           mean_of_var = nanmean(nanvar(cell_bins,[],1));
           Chirp_average(kk).stats(:,ii) = var_of_mean/mean_of_var;
           Chirp_average(kk).stats_max = nanmax(Chirp_average(kk).stats(:,ii));
            
       catch
           continue
       end
        
      end
       
  
    
    
end
 
% Clear empty entries in the structure
a = 1;
for ii = 1:nr_cells
    try
    if isempty(Chirp_average(a).cell_idx)
        Chirp_average(a) = [];
    else
        a = a+1;
    end
    catch
        continue
    end
end
 nr_cells_new = length(Chirp_average);
%Save averaged traces
out = sf_organizer(stim_idx,savepath,'variable_name','Chirp_average','variable',Chirp_average);
 
 

 
%% Create stat figures
overview_FFF = figure;
bar([Chirp_average.cell_idx],[Chirp_average.stats_max],'k');
title("Quality of FFF responses")
ylabel("Quality Index")
xlabel("Cell index")
hline(0.3,'-',"Threshold")
sf_organizer(add_info.stim_idx,savepath,'variable_name','overview_Chirp',...
    'variable',overview_FFF);

%Collect quality criteria data
FFF_stats = NaN(nr_colours,length(Chirp_average));
for ii = 1:nr_cells_new
    try
        FFF_stats(:,ii) = Chirp_average(ii).stats;
    catch
        continue
    end
    
end

detail_stats_FFF = heatbar([Chirp_average.cell_idx],FFF_stats);

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
total_bins = int64(nr_bins_unique(2));
Chirp_traces = NaN(nr_cells_new,total_bins);
for ii = 1:nr_cells_new
   try
   Chirp_trace = Chirp_average(ii).traces(:,2);
   Chirp_trace = Chirp_trace/nanmax(Chirp_trace,[],'all');
   Chirp_traces(ii,:) = reshape(Chirp_trace,[1,total_bins]);
   catch
       continue
   end
end
binsize = Bined_spikes(1).bins_info.binsize;
stim_duration = nr_bins_unique(2)*binsize;
x_values = (binsize:binsize:stim_duration);
all_traces_plot = heatbar(x_values,Chirp_traces,'Y',[Chirp_average.cell_idx],...
    'gap', false);

title("All traces heatmap")
xlabel("Time in s")
ylabel("Cell Index")
c = findall(all_traces_plot.fig.Children,'type','ColorBar');
c.Label.String = 'Spikes per bin';
sf_organizer(add_info.stim_idx,savepath,'variable_name','all_traces_plot',...
    'variable',all_traces_plot.fig);

 
 
 
 
end