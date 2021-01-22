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



spike_cell = squeeze(spiketimestamps(:,1,:));
for ii = 1:size(spike_cell,2)
    spike_cell(:,ii) = spike_cell(:,ii)-locs_begin_s(ii);
        
end

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