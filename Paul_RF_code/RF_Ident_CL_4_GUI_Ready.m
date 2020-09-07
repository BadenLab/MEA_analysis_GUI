%% RF Indentification Command Line 4 - GUI Ready
% 04,06,2020 Onwards
% Works with RF_Ident_fn_v4.m

clear all;
clc;

%% Choose whether to use parallel processing and number of cores
Parpool   = 2; % 1: Yes, 2: No.
Num_Cores = 6;

%% Make RF Identification Choices

% Choose RF Identification Method
% If multiple options are chosen
% p.RF_Ident_Meth_vec = [a,b,c,d];
% a: STA-SD method;
% b: LC method (Local Covariance);
% c: MI method (Mutual Information).
% d: SC method (Self Covariance).
p.RF_Ident_Meth_vec = [1,0,0,1];

% Choose whether to calculate STE_Full for use with LC/MI/SC % PAR Mod 27,08,2020
% 1: Yes,
% 2: No.
p.STE_Full_Choice = 1;

% Choose RF Quality Control (QC) and Thresholds
if p.RF_Ident_Meth_vec(1) == 1
    
    % Choose STA-SD QC type
    % 1. SD threshold;
    % 2. CI (Confidence Interval).
    p.STA_SD_QC_Type = 1;
    
    if p.STA_SD_QC_Type == 1
        % Choose SD Threshold
        p.STA_SD_Thresh = 5; % MEA Data 1 (50px) -- 5 seems good. 
    elseif p.STA_SD_QC_Type == 2
        % Choose Upper Percentile
        p.STA_SD_CI_Per = 1; % default: 1 (%)
    end
         
end

if p.RF_Ident_Meth_vec(2)==1
    
    % Choose LC QC type
    % 1. SD threshold;
    % 2. CI (Confidence Interval).
    p.LC_QC_Type = 1;
    
    if p.LC_QC_Type == 1
        % Choose SD Threshold
        p.LC_Thresh = 3; % default: 6 (SDs)
    elseif p.LC_QC_Type == 2
        % Choose Upper Percentile
        p.LC_CI_Per = 1; % default: 1 (%)
    end
    
end

if p.RF_Ident_Meth_vec(3)==1
    
    % Choose MI QC type
    % 1. SD threshold;
    % 2. CI (Confidence Interval).
    p.MI_QC_Type = 1;
    
    if p.MI_QC_Type == 1
        % Choose SD Threshold
        p.MI_Thresh = 3; % default: 6 (SDs)
    elseif p.MI_QC_Type == 2
        % Choose Upper Percentile
        p.MI_CI_Per = 1; % default: 1 (%)
    end
    
end

if p.RF_Ident_Meth_vec(4) == 1
    
    % Choose SC QC type
    % 1. SD threshold;
    % 2. CI (Confidence Interval).
    p.SC_QC_Type = 1;
    
    if p.SC_QC_Type == 1
        % Choose SD Threshold
        p.SC_Thresh = 7; % MEA Data 1 (50px): 7 seems good.
    elseif p.SC_QC_Type == 2
        % Choose Upper Percentile
        p.SC_CI_Per = 1; % default: 1 (%)
    end
         
end

% Choose RF Type
% p.RF_Type = [a,b,c];
% a. Box;
% b. All Significant Pixels;
% c. Gaussian.
p.RF_Type = [1,1,1];

% Choose RF Parameters
if p.RF_Type(1) == 1   % Box
    % Choose number of rings of pixels around central RF pixel
    p.RF_layers = 1;
    % Choose number of rings of pixels around significant Full RF pixels
    p.FullRF_layers = 1;
end

if p.RF_Type(3) == 1   % Gaussian
    % Choose number of standard deviations at which to set RF contour
    p.RF_SDs = 2; % 2 default
end

%% Post RF Calc Options

% Choose whether to plot RF indentification results
% 1. Yes;
% 2. No.
p.Plot_Choice = 1;

%% Make Data Related Choices

%%% Choose MEA cell(s)
% Choose all cells or a subset
% 1. All cells;
% 2. Subset of cells.
Cell_Choice = 2;
if Cell_Choice == 2
    First_Cell  = 1;
    Last_Cell   = 10; % 10, 1000, 5000
    %Num_Cells   = Last_Cell - First_Cell + 1; % Never used
    Cell_Choice_vec = First_Cell:Last_Cell; % '5000' cells in 'Colour_noise50px_data.mat'
end

% Choose whether to work with time bins or define own time grid
% 1. work with stimulus frames;
% 2. define own time grid.
p.Time_Choice = 1;

% Choose whether to subtract the mean raw stimulus in the STA calculation
% 1. don't subtract mean raw stim;
% 2. subtract mean raw stim.
p.STA_Choice  = 2;

% Choose whether to calc mean raw stim from stimulus frames (much quicker)
% or own time grid (much slower but a bit more accurate)
% 1. Calc from stim frames;
% 2. Calc from own time grid.
if p.Time_Choice == 2 && p.STA_Choice == 2
    p.Mean_Stim_Choice  = 1;
end

% Choose whether to use local covariance or local covariance difference
% 1. use local covariance directly;
% 2. use local covariance difference.
if p.RF_Ident_Meth_vec(2) == 1
    LocCovar_Choice = 1;
end

% Choose whether to calc raw local covar from stimulus frames (much quicker)
% or own time grid (much slower but a bit more accurate)
% 1. Calc from stim frames;
% 2. Calc from own time grid.
if p.Time_Choice == 2 && p.RF_Ident_Meth_vec(2) == 1 && LocCovar_Choice == 2
    RawLocCovar_Choice  = 1;
end

% Choose length of time window in which to find the STEs
if p.Time_Choice == 2
    STE_int      = 2;  % sec (default 2) (if create own interpolated values)
end

% Choose stimulus resolution (number of points sampled from stimulus time
% window)
p.Num_STE_bins = 10; % Default: 5 (Tom data, Time_Choice 1) or 10 (Tom data,Time_Choice 2); 20 (Marvin data)


%% Load Data

%load Marvin_Data_Chicken_24_08_2019_PARMod.mat;
%load Marvin_Data_Chicken_06_12_2019_PARMod.mat;

%load Colour_noise50px_data_PARMod.mat;
%load Colour_noise50px_data_PARMod2.mat;

load Colour_noise50px_data.mat;
%load 30pxNoise_subset.mat; % There is something wrong with this data set, it gives no RF.
%load 20pxNoise_subset.mat;

if Cell_Choice == 1 % All cells
    spike_times_mat    = spiketimestamps;
    Cell_Choice_vec    = 1:1:size(spiketimestamps,2);
    True_cell_index_vec = Cell_Choice_vec(any(~isnan(spike_times_mat)));
else % Cell_Choice == 2 % Subset of cells
    spike_times_mat    = spiketimestamps(:,Cell_Choice_vec);
    % Remove NaN columns (non-responsive cells)
    True_cell_index_vec = Cell_Choice_vec(any(~isnan(spike_times_mat))); % Record indices of remaining cells
end

spike_times_mat(:,~any(~isnan(spike_times_mat))) = [];
True_Num_Cells     = size(spike_times_mat,2);

trigCh_vec         = Ch_new.trigger_ch;
min_trigCh_vec     = min(trigCh_vec);
max_trigCh_vec     = max(trigCh_vec);
trigThreshFac      = 0.05; % 0.1
trigHigh_vec       = double(trigCh_vec > min_trigCh_vec + trigThreshFac*(max_trigCh_vec-min_trigCh_vec));
[~,trig_index_vec] = findpeaks(trigHigh_vec); % 'MinPeakProminence',#,'MinPeakDistance',#

sampling_freq      = Ch_new.SamplingFrequency;
sampling_int       = 1/sampling_freq;

trig_times_vec_temp = (trig_index_vec-1)*sampling_int;
% I take it that the rising edge starts at the timepoint directly
% before the jump upwards - hence '(trig_index_vec-1)' rather than 'trig_index_vec'
trig_times_vec      = trig_times_vec_temp - trig_times_vec_temp(1); % set first frame time to 0 s.

clear spiketimestamps Ch_new; % Remove original data files as no longer needed.

% Check that triggers are correctly identified.
% figure; plot(trigCh_vec); hold on;
% plot(trig_index_vec,4085*ones(1,length(trig_index_vec)),'ro');

stimulus_arr(:,:,1) = importdata('Red_Noise.txt',',')';
stimulus_arr(:,:,2) = importdata('Green_Noise.txt',',')';
stimulus_arr(:,:,3) = importdata('Blue_Noise.txt',',')';
stimulus_arr(:,:,4) = importdata('UV_Noise.txt',',')';
stimulus_arr(:,1,:) = []; % stimulus_arr(1,:,:) before transposed above
stimulus_arr        = reshape(stimulus_arr,[40,40,6000,4]);

p.Spectral_Dim = size(stimulus_arr,4);

% if size(stimulus_arr,3) ~= length(trig_times_vec)
%     disp('NB: There is a mismatch between the number of stimulus frames and the number of triggers!');
% end

p.stim_frames = size(stimulus_arr,3);   % PAR Mod 27,08,2020 (moved from below, was stim_frames)
p.Num_trigs   = length(trig_times_vec); % PAR Mod 27,08,2020

if p.Num_trigs > p.stim_frames % Frozen noise repeated case % PAR Mod 27,08,2020 (whole if statement)
    
    Num_FNoise_rep      = p.Num_trigs/p.stim_frames;
    p.Num_FNoise_rep_ceil = ceil(Num_FNoise_rep);
    max_trig_int        = max(diff(trig_times_vec(1:p.stim_frames)));
    inter_noise_int_vec = NaN(p.Num_FNoise_rep_ceil-1,1);
    for i = 1:p.Num_FNoise_rep_ceil-1
        inter_noise_int_vec(i) = trig_times_vec(i*p.stim_frames+1)-trig_times_vec(i*p.stim_frames);
    end
    
end


noise_length         = min(p.Num_trigs,p.stim_frames); % PAR Mod 27,08,2020 Do this rather than just stim_frames incase only part of one chunk used
trig_times_vec_trunc = trig_times_vec(1:noise_length); % PAR Mod 27,08,2020 (truncated to first noise chunk)

%%% Set Parameters

% Stimulus ppties
stim_freq = 20;          % Hz
stim_int  = 1/stim_freq; % sec

if p.Time_Choice == 1 % stimulus frames
    
    p.stim_timesample_vec = linspace(-stim_int*p.Num_STE_bins,-stim_int,p.Num_STE_bins);
    
    % Mimimum time intervals after first and last triggers in which to consider
    % spikes (to ensure spikes occur following stimulus presentation across
    % full STE window)
    min_start_int = trig_times_vec(p.Num_STE_bins+1);  % sec (was 'stim_int*p.Num_STE_bins', this gave 0 indices when calc STA/STE_Full)
    min_end_int   = 2*stim_int;               % sec
    
elseif p.Time_Choice == 2 % define own time grid
    
    % Pre-spike stimulus segment time vector
    p.stim_timesample_vec = linspace(-STE_int,0,p.Num_STE_bins);
    
    % Mimimum time intervals after first and last triggers in which to consider
    % spikes (to ensure spikes occur following stimulus presentation across
    % full STE window)
    min_start_int = STE_int;  % sec
    min_end_int   = stim_int; % sec
    
end


% Trigger times are the times each new frame in the stimulus starts
% The first frame starts at trig_times_vec(1) and the last frame starts
% at trig_times_vec(end) and hence ends at trig_times_vec(end) + stim_int.
% I had made last_trig_time  = trig_times_vec(end) + stim_int; but I'm
% changing it to last_trig_time  = trig_times_vec(end); since the
% '+min_end_int' in spike_times_vec should account for spiking after stim
% end.
first_trig_time       = trig_times_vec(1);
last_trig_time        = trig_times_vec(end); % trig_times_vec(end) + stim_int
last_trig_time_trunc  = trig_times_vec_trunc(end); % PAR Mod 27,08,2020

p.stim_rows        = size(stimulus_arr,1);
p.stim_columns     = size(stimulus_arr,2);
p.stim_pixels      = p.stim_rows*p.stim_columns;
%stim_frames        = size(stimulus_arr,3); % PAR Mod 27,08,2020 (moved above)

Num_stim_spat_dim  = p.stim_rows*p.stim_columns;

Num_STE_Stixels = Num_stim_spat_dim*p.Num_STE_bins;

if p.STA_Choice  == 2 % subtract mean raw stim
    if p.Time_Choice == 1 || p.Mean_Stim_Choice  == 1 % stimulus frames || Calc from stim frames
        p.Num_Raw_Stim = p.stim_frames - p.Num_STE_bins + 1;
        % Total number of raw stimuli displayed
    elseif p.Time_Choice == 2                       % define own time grid
        p.Num_Raw_Stim = ceil((last_trig_time_trunc - first_trig_time + stim_int - STE_int)/(STE_int/p.Num_STE_bins)) + 1; % PAR Mod 27,08,2020 'last_trig_time' --> 'last_trig_time_trunc'
        % This is the total time over which stimuli are played 'last_trig_time-first_trig_time+stim_int'
        % minus the length of an STE 'STE_int' divided by the time between time
        % samples 'STE_int/p.Num_STE_bins' rounded upwards to the nearest integer, can times 10 to try to
        % make sure don't miss any stimulus combinations and plus 1 to count
        % the first stimulus before shifting the window.
        % 111.208542 seconds for model 1 without *10 so do this o/w takes too
        % long.
        p.Num_Raw_Stim_t_vec = linspace(first_trig_time + STE_int,last_trig_time_trunc + stim_int,p.Num_Raw_Stim); % PAR Mod 27,08,2020 'last_trig_time' --> 'last_trig_time_trunc'
    end
end

%%% Preliminary Plots

% Find and plot the distribution of trigger intervals
actual_stim_int = diff(trig_times_vec);
mean_stim_int   = mean(actual_stim_int);
stddev_stim_int = std(actual_stim_int);
figure;
histogram(actual_stim_int); hold on;
v1 = vline(mean_stim_int,'r');
v2 = vline(mean_stim_int-stddev_stim_int,'g');
v3 = vline(mean_stim_int+stddev_stim_int,'g');
xlabel('trigger int.');
ylabel('freq.');
title('trigger int. dist.');
set(v1,'LineWidth',1.5);
set(v2,'LineWidth',1.5);
set(v3,'LineWidth',1.5);
set(gca,'FontSize',12);
set(gcf,'color','w');


%% Call RF Identification Function
% MEA data: stimulus_array := 40 rows, 40 columns,6000 frames, 4 colours

RF_Ident = cell(True_Num_Cells,1);

if p.RF_Ident_Meth_vec(1) == 1 % STA-SD method
    if p.STA_Choice  == 2 % subtract average stim
        p.mean_raw_stim_arr = NaN(p.stim_rows,p.stim_columns,p.Num_STE_bins,p.Spectral_Dim);
        for i = 1:p.Spectral_Dim
            p.mean_raw_stim_arr(:,:,:,i) = mean_raw_stim_SpaceTime_fn_v2(stimulus_arr(:,:,:,i),trig_times_vec_trunc,p); % PAR Mod 27,08,2020 was 'mean_raw_stim_SpaceTime_fn' and 'trig_times_vec'
        end
    end
end

tic;
if Parpool == 1 % Parpool on
    
    parpool(Num_Cores);
    
    parfor i = 1:True_Num_Cells % for/parfor
        
        spike_times_vec_loop      = spike_times_mat(:,i);
        spike_times_vec_loop(isnan(spike_times_vec_loop)) = []; % spiketimestamps(~isnan(spiketimestamps(:,Cell_Choice)),Cell_Choice);
        spike_times_vec_loop      = spike_times_vec_loop(spike_times_vec_loop>(first_trig_time + min_start_int));
        spike_times_vec_loop      = spike_times_vec_loop(spike_times_vec_loop<(last_trig_time + min_end_int));
        
        if p.Num_trigs > p.stim_frames % Frozen noise repeated case % PAR Mod 27,08,2020 (whole if statement)
            
            if sum(inter_noise_int_vec <= max_trig_int) == p.Num_FNoise_rep_ceil-1 % If the frozen noise repeates w/o a gap
                
                % No need to remove further spikes in this case
                
                for j = 2:p.Num_FNoise_rep_ceil % map spikes in repeated stimulus chunks to original chunk
                    
                    if j==p.Num_FNoise_rep_ceil
                        frame_end_loop = p.Num_trigs - p.stim_frames*(p.Num_FNoise_rep_ceil-1);
                    else
                        frame_end_loop = p.stim_frames;
                    end
                    
                    for k = 1:frame_end_loop
                        if (k < frame_end_loop) || (j < p.Num_FNoise_rep_ceil)
                            indices_loop = find((spike_times_vec_loop>trig_times_vec((j-1)*p.stim_frames + k))&&(spike_times_vec_loop<trig_times_vec((j-1)*p.stim_frames + k + 1)));
                            a_loop      = trig_times_vec(k);
                            b_loop      = trig_times_vec(k+1);
                            a_dash_loop = trig_times_vec((j-1)*p.stim_frames + k);
                            b_dash_loop = trig_times_vec((j-1)*p.stim_frames + k + 1);
                            spike_times_vec_loop(indices_loop) = a_loop + (spike_times_vec_loop(indices_loop) - a_dash_loop)*(b_loop - a_loop)/(b_dash_loop - a_dash_loop);
                        else % (k == frame_end_loop) && (j == p.Num_FNoise_rep_ceil) % No end time is given for the last frame so have to specify differently
                            indices_loop = find((spike_times_vec_loop>trig_times_vec((j-1)*p.stim_frames + k))&&(spike_times_vec_loop<trig_times_vec((j-1)*p.stim_frames + k) + stim_int));
                            a_loop      = trig_times_vec(k);
                            b_loop      = trig_times_vec(k+1);
                            a_dash_loop = trig_times_vec((j-1)*p.stim_frames + k);
                            b_dash_loop = trig_times_vec((j-1)*p.stim_frames + k) + stim_int;
                            spike_times_vec_loop(indices_loop) = a_loop + (spike_times_vec_loop(indices_loop) - a_dash_loop)*(b_loop - a_loop)/(b_dash_loop - a_dash_loop);
                        end
                    end
                    
                end
                
            else % If there are gaps between (any) frozen noise chuncks
                
                for j = 1:p.Num_FNoise_rep_ceil-1 % remove spikes with stimulus windows that fall between the noise chunks
                    spike_times_vec_loop((spike_times_vec_loop>(trig_times_vec(j*p.stim_frames) + min_end_int))&(spike_times_vec_loop<trig_times_vec(j*p.stim_frames+p.Num_STE_bins+1))) = [];
                end
                
                for j = 2:p.Num_FNoise_rep_ceil % map spikes in repeated stimulus chunks to original chunk
                    
                    if j==p.Num_FNoise_rep_ceil
                        frame_end_loop = p.Num_trigs - p.stim_frames*(p.Num_FNoise_rep_ceil-1);
                    else
                        frame_end_loop = p.stim_frames;
                    end
                    
                    for k = 1:frame_end_loop
                        if k<frame_end_loop % frame_end_loop<p.stim_frames
                            indices_loop = find((spike_times_vec_loop>trig_times_vec((j-1)*p.stim_frames + k))&&(spike_times_vec_loop<trig_times_vec((j-1)*p.stim_frames + k + 1)));
                            a_loop      = trig_times_vec(k);
                            b_loop      = trig_times_vec(k+1);
                            a_dash_loop = trig_times_vec((j-1)*p.stim_frames + k);
                            b_dash_loop = trig_times_vec((j-1)*p.stim_frames + k + 1);
                            spike_times_vec_loop(indices_loop) = a_loop + (spike_times_vec_loop(indices_loop) - a_dash_loop)*(b_loop - a_loop)/(b_dash_loop - a_dash_loop);
                        else % k==frame_end_loop % frame_end_loop==p.stim_frames % No end time is given for the last frame so have to specify differently
                            indices_loop = find((spike_times_vec_loop>trig_times_vec((j-1)*p.stim_frames + k))&&(spike_times_vec_loop<trig_times_vec((j-1)*p.stim_frames + k) + stim_int));
                            a_loop      = trig_times_vec(k);
                            b_loop      = trig_times_vec(k+1);
                            a_dash_loop = trig_times_vec((j-1)*p.stim_frames + k);
                            b_dash_loop = trig_times_vec((j-1)*p.stim_frames + k) + stim_int;
                            spike_times_vec_loop(indices_loop) = a_loop + (spike_times_vec_loop(indices_loop) - a_dash_loop)*(b_loop - a_loop)/(b_dash_loop - a_dash_loop);
                        end
                    end
                    
                end
                
            end
            
        end
        
        %p.length_spike_times      = length(spike_times_vec_loop);
        length_spike_times_loop   = length(spike_times_vec_loop); % To allow parfor
        
        RF_Ident{i} = RF_Ident_fn_v5(stimulus_arr,trig_times_vec_trunc,spike_times_vec_loop,length_spike_times_loop,p); % PAR Mod 27,08,2020 --> was 'RF_Ident_fn_v4' and 'trig_times_vec'
        
        disp(i);
        
    end
    
    delete(gcp('nocreate'));
    
else %  Parpool == 2 % Parpool off
    
    for i = 1:True_Num_Cells
        
        spike_times_vec_loop      = spike_times_mat(:,i);
        spike_times_vec_loop(isnan(spike_times_vec_loop)) = []; % spiketimestamps(~isnan(spiketimestamps(:,Cell_Choice)),Cell_Choice);
        spike_times_vec_loop      = spike_times_vec_loop(spike_times_vec_loop>(first_trig_time + min_start_int));
        spike_times_vec_loop      = spike_times_vec_loop(spike_times_vec_loop<(last_trig_time + min_end_int));
        
        if p.Num_trigs > p.stim_frames % Frozen noise repeated case % PAR Mod 27,08,2020 (whole if statement)
            
            if sum(inter_noise_int_vec <= max_trig_int) == p.Num_FNoise_rep_ceil-1 % If the frozen noise repeates w/o a gap
                
                % No need to remove further spikes in this case
                
                for j = 2:p.Num_FNoise_rep_ceil % map spikes in repeated stimulus chunks to original chunk
                    
                    if j==p.Num_FNoise_rep_ceil
                        frame_end_loop = p.Num_trigs - p.stim_frames*(p.Num_FNoise_rep_ceil-1);
                    else
                        frame_end_loop = p.stim_frames;
                    end
                    
                    for k = 1:frame_end_loop
                        if (k < frame_end_loop) || (j < p.Num_FNoise_rep_ceil)
                            indices_loop = find((spike_times_vec_loop>trig_times_vec((j-1)*p.stim_frames + k))&(spike_times_vec_loop<trig_times_vec((j-1)*p.stim_frames + k + 1)));
                            a_loop      = trig_times_vec(k);
                            b_loop      = trig_times_vec(k+1);
                            a_dash_loop = trig_times_vec((j-1)*p.stim_frames + k);
                            b_dash_loop = trig_times_vec((j-1)*p.stim_frames + k + 1);
                            spike_times_vec_loop(indices_loop) = a_loop + (spike_times_vec_loop(indices_loop) - a_dash_loop)*(b_loop - a_loop)/(b_dash_loop - a_dash_loop);
                        else % (k == frame_end_loop) && (j == p.Num_FNoise_rep_ceil) % No end time is given for the last frame so have to specify differently
                            indices_loop = find((spike_times_vec_loop>trig_times_vec((j-1)*p.stim_frames + k))&(spike_times_vec_loop<trig_times_vec((j-1)*p.stim_frames + k) + stim_int));
                            a_loop      = trig_times_vec(k);
                            b_loop      = trig_times_vec(k+1);
                            a_dash_loop = trig_times_vec((j-1)*p.stim_frames + k);
                            b_dash_loop = trig_times_vec((j-1)*p.stim_frames + k) + stim_int;
                            spike_times_vec_loop(indices_loop) = a_loop + (spike_times_vec_loop(indices_loop) - a_dash_loop)*(b_loop - a_loop)/(b_dash_loop - a_dash_loop);
                        end
                    end
                    
                end
                
            else % If there are gaps between (any) frozen noise chuncks
                
                for j = 1:p.Num_FNoise_rep_ceil-1 % remove spikes with stimulus windows that fall between the noise chunks
                    spike_times_vec_loop((spike_times_vec_loop>(trig_times_vec(j*p.stim_frames) + min_end_int))&(spike_times_vec_loop<trig_times_vec(j*p.stim_frames+p.Num_STE_bins+1))) = [];
                end
                
                for j = 2:p.Num_FNoise_rep_ceil % map spikes in repeated stimulus chunks to original chunk
                    
                    if j==p.Num_FNoise_rep_ceil
                        frame_end_loop = p.Num_trigs - p.stim_frames*(p.Num_FNoise_rep_ceil-1);
                    else
                        frame_end_loop = p.stim_frames;
                    end
                    
                    for k = 1:frame_end_loop
                        if k<frame_end_loop % frame_end_loop<p.stim_frames
                            indices_loop = find((spike_times_vec_loop>trig_times_vec((j-1)*p.stim_frames + k))&&(spike_times_vec_loop<trig_times_vec((j-1)*p.stim_frames + k + 1)));
                            a_loop      = trig_times_vec(k);
                            b_loop      = trig_times_vec(k+1);
                            a_dash_loop = trig_times_vec((j-1)*p.stim_frames + k);
                            b_dash_loop = trig_times_vec((j-1)*p.stim_frames + k + 1);
                            spike_times_vec_loop(indices_loop) = a_loop + (spike_times_vec_loop(indices_loop) - a_dash_loop)*(b_loop - a_loop)/(b_dash_loop - a_dash_loop);
                        else % k==frame_end_loop % frame_end_loop==p.stim_frames % No end time is given for the last frame so have to specify differently
                            indices_loop = find((spike_times_vec_loop>trig_times_vec((j-1)*p.stim_frames + k))&&(spike_times_vec_loop<trig_times_vec((j-1)*p.stim_frames + k) + stim_int));
                            a_loop      = trig_times_vec(k);
                            b_loop      = trig_times_vec(k+1);
                            a_dash_loop = trig_times_vec((j-1)*p.stim_frames + k);
                            b_dash_loop = trig_times_vec((j-1)*p.stim_frames + k) + stim_int;
                            spike_times_vec_loop(indices_loop) = a_loop + (spike_times_vec_loop(indices_loop) - a_dash_loop)*(b_loop - a_loop)/(b_dash_loop - a_dash_loop);
                        end
                    end
                    
                end
                
            end
            
        end
        
        %p.length_spike_times      = length(spike_times_vec_loop);
        length_spike_times_loop   = length(spike_times_vec_loop); % To allow parfor
        
        RF_Ident{i} = RF_Ident_fn_v5(stimulus_arr,trig_times_vec_trunc,spike_times_vec_loop,length_spike_times_loop,p); % PAR Mod 27,08,2020 --> was 'RF_Ident_fn_v4' and 'trig_times_vec'
        
        disp(i);
        
    end
    
end
toc;

%% Save Data
%save('Save_Name.mat','-v7.3')

%save data_Test_#;
