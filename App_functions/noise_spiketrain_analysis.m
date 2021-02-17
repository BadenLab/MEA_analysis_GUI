function out = noise_spiketrain_analysis (savepath, add_info)

nr_stimuli = length(add_info.stim_idx);

if nr_stimuli > 1
    error("This script only works for a single stimulus at a time")
end

S = load(savepath,'Ch');
trigger_ch = S.Ch.Ch01_02;

cell_idx = add_info.kernel_raster_cell;

noise_begin_fr = add_info.stim_begin*S.Ch.SamplingFrequency;
noise_end_fr = add_info.stim_end*S.Ch.SamplingFrequency;

trigger_ch = trigger_ch(noise_begin_fr:noise_end_fr);

noise_trigger_norm = trigger_ch(1,:) > 2500;
%find the peaks in the differences (which correspond to the beginning of a
%trigger event
[~,locs] = findpeaks(gather(double(noise_trigger_norm)),'MinPeakProminence',1,'MinPeakDistance',178);
locs_temp = 1;
locs_temp(2:length(locs)+1) = locs;
locs = locs_temp;
locs_diff = diff(locs);
%Find locs which indicate beginning of chunk of frozen noise
locs_diff_norm = ceil05(locs_diff,0.05);
locs_begin = locs_diff_norm > 5*mean(locs_diff_norm);
locs_end_idx = locs(locs_begin);
locs_end_idx(end+1) = locs(end);
locs_end_s = (locs_end_idx+noise_begin_fr)/S.Ch.SamplingFrequency;

locs_idx = find(locs_begin);
locs_idx = locs_idx + 1;
locs_begin_idx_temp = locs(locs_idx);
locs_begin_idx = ones(1,length(locs_begin_idx_temp)+1);
locs_begin_idx(2:end) = locs_begin_idx_temp;
%Get it in seconds
locs_begin_s = (locs_begin_idx+noise_begin_fr)/S.Ch.SamplingFrequency;
locs_s = locs/S.Ch.SamplingFrequency;
%% Load actual stimulus sequence (to see correlations between different colours)
%Check how many frames per repeat
begin_second_locs = find(locs_begin,1,'first');
nr_frames = begin_second_locs +1;
stimulus_arr = load_noise_from_hdf5(string(add_info.settings.location.noise),...
    true,1,double(nr_frames));
nr_pixel = size(stimulus_arr,2);

%lets estimate the overall luminance (nr of active pixel per trigger)
nr_active_pixel = squeeze(sum(stimulus_arr,2));
nr_active_pixel_sm = movmean(nr_active_pixel,10);
nr_active_pixel_sm_relative = nr_active_pixel_sm/nr_pixel;
%Calculate stimulus frequency
noise_interval = median(locs_diff_norm)/S.Ch.SamplingFrequency;







%% Load spikes 


%Load spikes according to noise chunks
chunk_info.stim_begin = locs_begin_s;
chunk_info.stim_end = locs_end_s;
spiketimestamps = load_spiketimestamps_app(savepath,chunk_info,[cell_idx cell_idx]);
%spiketimestamps = spiketimestamps - add_info.stim_begin;


locs_idx_temp = 1;
locs_idx_temp(2:length(locs_idx)+1) = locs_idx;
locs_idx_temp(end+1) = length(locs_s)+1;
spike_cell = squeeze(spiketimestamps(:,1,:));
for ii = 1:size(spike_cell,2)
    spike_cell(:,ii) = spike_cell(:,ii)-locs_begin_s(ii);
    locs_subset(:,ii) = locs_s(locs_idx_temp(ii):locs_idx_temp(ii+1)-1);
    locs_subset(:,ii) = locs_subset(:,ii) - locs_subset(1,ii);
    
end


for ii = 1:size(spike_cell,2)
   locs_subset(begin_second_locs+1,ii) = locs_subset(end,ii)+noise_interval; 
end



%% Alignment of spikes (horror show)
spikes_aligned = NaN(length(locs_subset),size(spike_cell,2));
for ii = 1:size(spike_cell,2)
   spikes_temp = spike_cell(:,ii); 
   a = 1;
   for kk = 1:begin_second_locs
        spikes_match = logical((spikes_temp > locs_subset(kk,ii)).*(...
            spikes_temp<= locs_subset(kk+1,ii))); 
        try
        spikes_aligned(a,ii) = spikes_temp(spikes_match)-locs_subset(kk,ii)+kk;
        a = a+1;
        catch
            continue
        end
   end
    
end




%% Bin Data
%Transform data to a series of 0 and 1 at a given bin size
binsize = 0.001; %To create a matrix with 0 and 1 (spike or no spike)
maxspikes = nanmax(spikes_aligned,[],'all');
if isnan(maxspikes)
    %continue
end
nr_repeats = size(spikes_aligned,2);
nr_bins = ceil(maxspikes/binsize);

max_bin_time = binsize*nr_bins; %Last edge of the histcount

edges = (binsize:binsize:max_bin_time);

cell_histcount = zeros(nr_bins-1,nr_repeats);

%Bin all traces
for ii = 1:nr_repeats
    try
    cell_histcount(:,ii) = histcounts(spikes_aligned(:,ii),edges);
    catch
        continue
    end
end



%% convolve with gaussian kernel
%each spike falls into a certain time-probability window determined by a
%gaussian kernel

w = gausswin(2500);
spikes_smoothed = zeros(size(edges,2)-1,nr_repeats);
for ii = 1:nr_repeats
    
    spikes_smoothed(:,ii) = conv(cell_histcount(:,ii),w,'same');
    spikes_smoothed(:,ii) = spikes_smoothed(:,ii)/max(spikes_smoothed(:,ii));
   
end

%% Plot overview
% Calculate average trace (old stuff)
% point wise multiplication
spikes_smoothed1 = spikes_smoothed+1;
trace_average = [];
trace_average(:,1) = spikes_smoothed1(:,1);
for ii =1:size(spikes_smoothed1,2)-1
    trace_average = trace_average(:).*spikes_smoothed1(:,ii+1);
end
trace_average = trace_average/2-1;
%Normalize average trace
%trace_average = trace_average/max(trace_average);


figure
    ax = [];
    ax(1) = subplot(nr_repeats+2,1,1);
    plot_raster(spikes_aligned,1);

for ii = 1:nr_repeats
   ax(ii+1) = subplot(nr_repeats+2,1,ii+1);
   
   plot(edges(1:end-1),spikes_smoothed(:,ii),'k');
   
end

ax(nr_repeats+2) = subplot(nr_repeats+2,1,nr_repeats+2);
plot(edges(1:end-1),trace_average,'k')
linkaxes(ax,'x')



%% Plot
figure;
ax1 = subplot(3,1,1);
plot_raster(spike_cell,'ax',ax1);
title(['Cell ',num2str(cell_idx)])
ylabel("Repeats")
set(ax1,'xticklabel',[])
ax2 = subplot(3,1,2);
hold on
title("Stimulus trace")
xvalues = (0:1/S.Ch.SamplingFrequency:(locs_end_s(1)-locs_begin_s(1)));
trigger_plot = trigger_ch(locs_begin_idx(1):locs_end_idx(1));
%Sometimes there is a slight difference between x values due to how the 
%series is created above, this is fixed in this if statement
if length(xvalues) < length(trigger_plot)
    tr_p_length = length(trigger_plot);
    xvalues(end:tr_p_length) = xvalues(end);
end

plot(xvalues,trigger_plot,'k')
ylabel("Trigger signal, a.u.")


ax3 = subplot(3,1,3);
hold on
title("Stimulus mean luminance")
xvalues = (0:noise_interval:noise_interval*nr_frames-noise_interval);
colour = {'red','green','blue','magenta'};
for ii = 1:size(nr_active_pixel_sm,2)
    plot(xvalues,nr_active_pixel_sm_relative(:,ii),'Color',colour{1,ii})
end
ylabel("Proportion of active channel by colour, a.u.")


xlabel("Time in s")



linkaxes([ax1,ax2,ax3],'x')
out = 1;


end