%% RF Indentification Command Line 4 - GUI Ready
% 04,06,2020 Onwards
% Works with RF_Ident_fn_v4.m
function out = RF_Ident_CL_4_GUI_Ready (savepath, add_info)
%% Get inputs
kernel_settings = add_info.settings.kernel_new;
M = matfile(savepath,'Writable',false);
cell_indices = M.cell_indices;

%% Choose whether to use parallel processing and number of cores
Parpool = add_info.settings.parpro; %Decided in the gui
%Num_Cores = 4; %implement later

%% Make RF Identification Choices

% Choose RF Identification Method
% If multiple options are chosen
% p.RF_Ident_Meth_vec = [a,b,c,d];
% a: STA-SD method;
% b: LC method (Local Covariance);
% c: MI method (Mutual Information).
% d: SC method (Self Covariance).
p.RF_Ident_Meth_vec = [kernel_settings.SS,0,0,kernel_settings.SC];

% Choose whether to calculate STE_Full for use with LC/MI/SC % PAR Mod 27,08,2020
% 1: Yes,
% 2: No.
p.STE_Full_Choice = 1;%2;

% Choose RF Quality Control (QC) and Thresholds
if p.RF_Ident_Meth_vec(1) == 1
    
    % Choose STA-SD QC type
    % 1. SD threshold;
    % 2. CI (Confidence Interval).
    p.STA_SD_QC_Type = 1;
    
    if p.STA_SD_QC_Type == 1
        % Choose SD Threshold
        p.STA_SD_Thresh = kernel_settings.SS_STD; % MEA Data 1 (50px) -- 5 seems good. 
    elseif p.STA_SD_QC_Type == 2
        % Choose Upper Percentile
        p.STA_SD_CI_Per = kernel_settings.SS_CI_Upper; % default: 1 (%)
    end
         
end

if p.RF_Ident_Meth_vec(2)==1 %Remeains to be implemented @mars
    
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

if p.RF_Ident_Meth_vec(3)==1 %Remains to be implemented @mars
    
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
        p.SC_Thresh = kernel_settings.SC_STD;  % MEA Data 1 (50px): 7 seems good.
    elseif p.SC_QC_Type == 2
        % Choose Upper Percentile
        p.SC_CI_Per = kernel_settings.SC_CI_Lower; % default: 1 (%)
    end
         
end

% Choose RF Type
% p.RF_Type = [a,b,c];
% a. Box;
% b. All Significant Pixels;
% c. Gaussian.
p.RF_Type = [kernel_settings.RF_Box,kernel_settings.RF_AS,kernel_settings.RF_Gaus];

% Choose RF Parameters
if p.RF_Type(1) == 1   % Box
    % Choose number of rings of pixels around central RF pixel
    p.RF_layers = kernel_settings.RF_RBoxes;
    % Choose number of rings of pixels around significant Full RF pixels
    p.FullRF_layers = 1; %mars: What does this mean?
end

if p.RF_Type(3) == 1   % Gaussian
    % Choose number of standard deviations at which to set RF contour
    p.RF_SDs = kernel_settings.RF_SDGaus; % 2 default
end

%% Post RF Calc Options

% Choose whether to plot RF indentification results
% 1. Yes;
% 2. No.
p.Plot_Choice = kernel_settings.plot_RF_overview;
p.Time_Plot_Choice = 2;
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
p.Time_Choice = 1; %Should be in the gui

% Choose whether to subtract the mean raw stimulus in the STA calculation
% 1. don't subtract mean raw stim;
% 2. subtract mean raw stim.
p.STA_Choice  = 2; %Should be in the gui

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
    RawLocCovar_Choice  = 1; %@mars This is unused?
end

% Choose length of time window in which to find the STEs
if p.Time_Choice == 2
    STE_int      = 2;  % sec (default 2) (if create own interpolated values)
else
    STE_int = 1;
end

% Choose stimulus resolution (number of points sampled from stimulus time
% window)
p.Num_STE_bins = kernel_settings.Num_STE_bins; % Default: 5 (Tom data, Time_Choice 1) or 10 (Tom data,Time_Choice 2); 20 (Marvin data)


%% Load Data
%mars: I added a new way to load data. It depends on which stimulus the
%user selects in the GUI. The spiketimestamps get loaded in batches, to
%prevent the memory from overfloating. It calls functions which come with
%the GUI and are used in other scripts of the GUI as well. 


%% Memory handling  @mars
%Check how many spiketimestamps are in recording
M = matfile(savepath,'Writable',false); %Using matfile here so we dont have 
%to load the actual file
spikesize = numel(M.spiketimestamps);
stx_before = size(M.spiketimestamps,2);

%Check how much memory is available
%One double entry is 8 bytes, we dont want to use more than one fourth of the
%available memory (So that returned data has place as well)
memory_available = memory;
memory_available = memory_available.MaxPossibleArrayBytes/4;
array_size_available = memory_available/8;

%Calculate how many batches we need to process all the data. If only one
%batch is need, all data will be processed as one set.

%Calculate batch size
batch_nr = ceil(spikesize/array_size_available);
batch_stx = ceil(stx_before/batch_nr);
batch_begins = (1:batch_stx:stx_before);
batch_ends = batch_begins(2:end)-1;
batch_ends(end+1) = stx_before;


%% Specifiy trigger channel @mars
%This has to be done only once for all batches, as they are from the same
%stimulus normally. I had to swap some of the variable declarations around
%to make them fit to the new loop structure. Trigger channel stuff is
%outside the loop, spiketimestamp stuff inside the loop (of batches).

%mars: Load stim time details
stim_idx = add_info.stim_idx;
stim_begin = add_info.stim_begin;
stim_end = add_info.stim_end;
ch = load(savepath,'-mat','Ch');
ch = ch.Ch;
Ch01_02 = ch.Ch01_02;
SamplingFrequency = ch.SamplingFrequency;
stimulus_nr = numel(stim_begin); %Nr of selected noise stimuli

for nn = 1:stimulus_nr %Loop over selected stimuli

%% Load Data modified @mars
%% Trigger channel for each stimulus    
noise_begin = stim_begin(nn)*SamplingFrequency;
noise_end = stim_end(nn)*SamplingFrequency;
noise_begin_s = stim_begin(nn);
%noise_end_s = stim_end(nn);

%@mars: Cell selection, not sure this should stay in the code
% if Cell_Choice == 1 % All cells
%     spike_times_mat    = spiketimestamps;
%     Cell_Choice_vec    = 1:1:size(spiketimestamps,2);
%     True_cell_index_vec = Cell_Choice_vec(any(~isnan(spike_times_mat)));
% else % Cell_Choice == 2 % Subset of cells
%     spike_times_mat    = spiketimestamps(:,Cell_Choice_vec);
%     % Remove NaN columns (non-responsive cells)
%     True_cell_index_vec = Cell_Choice_vec(any(~isnan(spike_times_mat))); % Record indices of remaining cells
% end

trigCh_vec(1,:) = Ch01_02(noise_begin:noise_end); %This cuts the
%trigger trace at stimulus begin and end

min_trigCh_vec     = min(trigCh_vec);
max_trigCh_vec     = max(trigCh_vec);
trigThreshFac      = 0.05; % 0.1
trigHigh_vec       = double(trigCh_vec > min_trigCh_vec +...
    trigThreshFac*(max_trigCh_vec-min_trigCh_vec));
[~,trig_index_vec_temp] = findpeaks(trigHigh_vec);% 'MinPeakProminence',#,'MinPeakDistance',#
trig_index_vec = zeros(1,size(trig_index_vec_temp,2)+1);
trig_index_vec(1,1) = 1;
trig_index_vec(1,2:end) = trig_index_vec_temp;
%trig_index_vec(1,end) = size(trigHigh_vec,2);
sampling_freq      = SamplingFrequency;
sampling_int       = 1/sampling_freq;

trig_times_vec_temp = (trig_index_vec-1)*sampling_int;
% I take it that the rising edge starts at the timepoint directly
% before the jump upwards - hence '(trig_index_vec-1)' rather than 'trig_index_vec'
trig_times_vec      = trig_times_vec_temp - trig_times_vec_temp(1); % set first frame time to 0 s.
%@mars: this somehow didnt work correctly, the first trigger is at time
%0, so it must be (1) in trig_times_vec.

%% Check for frozen repeats @mars
%@mars: We need to find out how often the stimulus sequence was repeated. 
trig_times_diff = ceil05(diff(trig_times_vec),0.05); %@mars: find beginning and
%end of noise repeats
noise_begins = trig_times_diff > 5*mean(trig_times_diff);
nr_frozen = nnz(noise_begins)+1;%@mars: +1 because diff finds the difference
noise_begins_idx = ones(1,nr_frozen);
noise_begins_idx(2:end) = find(noise_begins); %@mars: index in trig_index_vec_temp 
p.noise_begins_idx = noise_begins_idx;
% of the begining of each frozen noise repeat


if nr_frozen > 1
    %Check if the last repeat of the noise is complete
    nr_trig_first_repeat = noise_begins_idx(2);
    nr_trig_last_repeat = length(trig_index_vec)-noise_begins_idx(end)-1;
    
    if nr_trig_first_repeat == nr_trig_last_repeat
        trig_per_frozen = int64(numel(trig_times_vec)/nr_frozen);
    else %@mars: Only possibility is that the last repeat is shorter
        trig_per_frozen = int64(noise_begins_idx(2));
    end
   
                
else %If effectively no frozen noise is played but only on sequence
    trig_per_frozen = length(trig_times_vec);
    warning('No noise repeats detected, assuming only 1 noise sequence and no frozen repeats');
    nr_trig_last_repeat = 0;
    
end
 p.Gap_Ind = 2;
    
%@mars: Now we can load the sequence according to what was played
%@mars: Need to add part that deals with case of no gaps between noise
%chunks


p.stim_frames = double(trig_per_frozen);   % PAR Mod 27,08,2020 (moved from below, was stim_frames)
p.Num_trigs   = length(trig_times_vec); % PAR Mod 27,08,2020
if p.Num_trigs < p.Num_STE_bins+1 % PAR Mod 09,10,2020 (whole if statement)
    disp('Number of STE frames exceeds number of triggers; therefore, RF identification is not possible');
end
p.noise_length         = min(p.Num_trigs,p.stim_frames); % PAR Mod 27,08,2020 Do this rather than just stim_frames incase only part of one chunk used
trig_times_vec_trunc = trig_times_vec(1:p.noise_length);

Num_FNoise_rep      = nr_frozen; 
p.Num_FNoise_rep_ceil = ceil(Num_FNoise_rep);
max_trig_int        = max(diff(trig_times_vec(1:p.stim_frames)));
inter_noise_int_vec = NaN(p.Num_FNoise_rep_ceil-1,1);
for i = 1:p.Num_FNoise_rep_ceil-1
    inter_noise_int_vec(i) = trig_times_vec(i*p.stim_frames+1)-...
        trig_times_vec(i*p.stim_frames);
end

%% Noise sequence for each stimulus @mars
%Noise sequence is either loaded from HDF5 file or text file
noise_file = add_info.settings.location.noise;
%Check for file ending 
file_ending = get_fileformat(noise_file);

if strcmp(file_ending,'.h5')
stimulus_arr = load_noise_from_hdf5(string(noise_file),true,1,double(trig_per_frozen));
nr_boxes = int64(sqrt(size(stimulus_arr,2)));
nr_colours = int64(size(stimulus_arr,3));
stimulus_arr = permute(stimulus_arr,[2 1 3]); % PAR Mod 01,02,2021
stimulus_arr = reshape(stimulus_arr,[nr_boxes,nr_boxes,trig_per_frozen,nr_colours]);

elseif strcmp(file_ending,'.txt') %Old version of getting the sequence
%Get folder location
folder_ending = max(strfind(noise_file,'\'));
folder_dir = noise_file(1:folder_ending);

stimulus_arr(:,:,1) = importdata([folder_dir,'Red_Noise.txt'],',')';
stimulus_arr(:,:,2) = importdata([folder_dir,'Green_Noise.txt'],',')';
stimulus_arr(:,:,3) = importdata([folder_dir,'Blue_Noise.txt'],',')';
stimulus_arr(:,:,4) = importdata([folder_dir,'UV_Noise.txt'],',')';
stimulus_arr(:,1,:) = []; % stimulus_arr(1,:,:) before transposed above
stimulus_arr        = reshape(stimulus_arr,[40,40,6000,4]);
nr_colours = 4; %mars: has to be in the gui
end

p.Spectral_Dim = nr_colours;

fprintf(['Loaded ', num2str(trig_per_frozen),' frames from noise file']);



%% Old stimulus loading
% % Check that triggers are correctly identified.
% %figure; plot(trigCh_vec); hold on;
% %plot(trig_index_vec,4085*ones(1,length(trig_index_vec)),'ro');
% 
% stimulus_arr(:,:,1) = importdata('Red_Noise.txt',',')';
% stimulus_arr(:,:,2) = importdata('Green_Noise.txt',',')';
% stimulus_arr(:,:,3) = importdata('Blue_Noise.txt',',')';
% stimulus_arr(:,:,4) = importdata('UV_Noise.txt',',')';
% stimulus_arr(:,1,:) = []; % stimulus_arr(1,:,:) before transposed above
% stimulus_arr        = reshape(stimulus_arr,[40,40,6000,4]);
% 
% p.Spectral_Dim = size(stimulus_arr,4);
% 
% % if size(stimulus_arr,3) ~= length(trig_times_vec)
% %     disp('NB: There is a mismatch between the number of stimulus frames and the number of triggers!');
% % end
% 
% 
% % p.stim_frames = size(stimulus_arr,3);   % PAR Mod 27,08,2020 (moved from below, was stim_frames)
% % p.Num_trigs   = length(trig_times_vec); % PAR Mod 27,08,2020
% %@mars: This will not work and is kind of circular. How do we know how many
% %stim frames we need to load? Only based on the stimulus channel, but than
% %we can never have more trigger signals than frames by definition.


%% Set and check parameters

% Stimulus ppties
% stim_freq = 20;        % Hz @mars: This should be defined programatically and not hard coded

stim_freq_check = 1/ceil05(nanmean(diff(trig_times_vec(1:p.noise_length))),0.05); %@mars: This should work in most cases... but if
%there are too many outliers we might want to ask the user...

%Compare calculated frequency to the user input and through error if they
%dont match
stim_freq = add_info.settings.kernel_new.Hz;

if stim_freq ~= round(stim_freq_check)
    warning('User input stimulus frequency does not match trigger frequency')
end

stim_int  = 1/stim_freq; % sec
p.stim_int = stim_int;



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

%Num_stim_spat_dim  = p.stim_rows*p.stim_columns;

%Num_STE_Stixels = Num_stim_spat_dim*p.Num_STE_bins;

if p.STA_Choice  == 2 % subtract mean raw stim
    if p.Time_Choice == 1 || p.Mean_Stim_Choice  == 1 % stimulus frames || Calc from stim frames
        p.Num_Raw_Stim = p.stim_frames - p.Num_STE_bins + 1;
        % Total number of raw stimuli displayed
    elseif p.Time_Choice == 2                       % define own time grid
        p.Num_Raw_Stim = ceil((last_trig_time_trunc -...
            first_trig_time + stim_int - STE_int)/(STE_int/p.Num_STE_bins)) + 1; % PAR Mod 27,08,2020 'last_trig_time' --> 'last_trig_time_trunc'
        % This is the total time over which stimuli are played 'last_trig_time-first_trig_time+stim_int'
        % minus the length of an STE 'STE_int' divided by the time between time
        % samples 'STE_int/p.Num_STE_bins' rounded upwards to the nearest integer, can times 10 to try to
        % make sure don't miss any stimulus combinations and plus 1 to count
        % the first stimulus before shifting the window.
        % 111.208542 seconds for model 1 without *10 so do this o/w takes too
        % long.
        p.Num_Raw_Stim_t_vec = linspace(first_trig_time +...
            STE_int,last_trig_time_trunc + stim_int,p.Num_Raw_Stim); % PAR Mod 27,08,2020 'last_trig_time' --> 'last_trig_time_trunc'
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

%% Loop over batch
%Load spiketimes according to the batch size defined above and set the save
%parameters for the function output

h = waitbar(0,'Batch nr: ');
% for bb = 1:batch_nr
for bb = 1:batch_nr
    waitbar(bb/batch_nr,h,[num2str(bb),' of ', batch_nr]);
%Initiate the batch    
spiketimestamps = load_spiketimestamps_app (savepath, add_info,...
    'cell_subset',[batch_begins(bb),batch_ends(bb)]);
cell_indices_temp = cell_indices(batch_begins(bb):batch_ends(bb));

%Quality Check %Possible usefull?? maybe make this Gui element

% right_indices = quality_check_spikes(spiketimestamps,10); %This needs to be 
% %tested further;
% spiketimestamps = spiketimestamps(:,right_indices);
%set times relative to stimulus begin
spiketimestamps = spiketimestamps - noise_begin_s;
%spiketimestamps(:,~any(~isnan(spiketimestamps))) = [];

stx = size(spiketimestamps,2);
%sty = size(spiketimestamps,1);  

%cell_indices_temp = cell_indices_temp(1,right_indices);
    
%% Call RF Identification Function
% MEA data: stimulus_array := 40 rows, 40 columns,6000 frames, 4 colours
%a = 1;
%RF_Ident = cell(stx,1);

if p.RF_Ident_Meth_vec(1) == 1 % STA-SD method
    if p.STA_Choice  == 2 % subtract average stim
        p.mean_raw_stim_arr = NaN(p.stim_rows,p.stim_columns,p.Num_STE_bins,p.Spectral_Dim);
        for i = 1:p.Spectral_Dim
            p.mean_raw_stim_arr(:,:,:,i) = mean_raw_stim_SpaceTime_fn_v2(stimulus_arr(:,:,:,i),trig_times_vec_trunc,p); % PAR Mod 27,08,2020 was 'mean_raw_stim_SpaceTime_fn' and 'trig_times_vec'
        end
      
    end
end

tic;
if Parpool == 1 % Parpool on %@mars: This doesnt work if a parpool is already active
    
    %parpool(Num_Cores);
    
    parfor i = 1:stx % for/parfor
        i
        
        settings = add_info;  
        p_par = p;
        trig_times_vec_par = trig_times_vec; %@mars: added this declaration
        %to avoid uneccessary communication between the workers in the parloop.
        spike_times_vec_loop      = spiketimestamps(:,i);
        %Check if all values are NAN
        test_for_NaN = find(all(isnan(spike_times_vec_loop),1));
        if test_for_NaN
            RF_Ident(i).RF_results = NaN;
            RF_Ident(i).cell_idx = cell_indices_temp(i);
             
            RF_Ident(i).add_info = settings;
            
            continue
        end
        
        spike_times_vec_loop(isnan(spike_times_vec_loop)) = []; % spiketimestamps(~isnan(spiketimestamps(:,Cell_Choice)),Cell_Choice);
        spike_times_vec_loop      = spike_times_vec_loop(spike_times_vec_loop>(first_trig_time + min_start_int));
        %spike_times_vec_loop      = spike_times_vec_loop(spike_times_vec_loop<(last_trig_time + min_end_int));
        if (nr_trig_last_repeat >= p.Num_STE_bins + 1) || (nr_trig_last_repeat == 0)  
            % PAR Mod 17,09,2020 (whole if statement, replaces 'spike_times_vec_loop      = spike_times_vec_loop(spike_times_vec_loop<(last_trig_time + min_end_int));')
            spike_times_vec_loop      = spike_times_vec_loop(spike_times_vec_loop<(last_trig_time + min_end_int));
        else % Num_Trig_Final_NoiseChunk < p.Num_STE_bins + 1
            spike_times_vec_loop      = spike_times_vec_loop(spike_times_vec_loop<...
                (trig_times_vec(p.stim_frames*(p.Num_FNoise_rep_ceil-1)) + min_end_int));
        end
        
        
        
        
        
        if p_par.Num_trigs > p_par.stim_frames % Frozen noise repeated case % PAR Mod 27,08,2020 (whole if statement)
            
            if sum(inter_noise_int_vec <= max_trig_int) == p_par.Num_FNoise_rep_ceil-1 % If the frozen noise repeates w/o a gap
                
                % No need to remove further spikes in this case
                
                for j = 2:p_par.Num_FNoise_rep_ceil % map spikes in repeated stimulus chunks to original chunk
                    
                    if j==p_par.Num_FNoise_rep_ceil
                        frame_end_loop = p_par.Num_trigs - p_par.stim_frames*(p_par.Num_FNoise_rep_ceil-1);
                    else
                        frame_end_loop = p_par.stim_frames;
                    end
                    
                    for k = 1:frame_end_loop
                        if (k < frame_end_loop) || (j < p_par.Num_FNoise_rep_ceil)
                            indices_loop = find((spike_times_vec_loop>...
                                trig_times_vec_par((j-1)*p_par.stim_frames + k))&&...
                                (spike_times_vec_loop<trig_times_vec_par((j-1)*p_par.stim_frames + k + 1)));
                            a_loop      = trig_times_vec_par(k);
                            b_loop      = trig_times_vec_par(k+1);
                            a_dash_loop = trig_times_vec_par((j-1)*p_par.stim_frames + k);
                            b_dash_loop = trig_times_vec_par((j-1)*p_par.stim_frames + k + 1);
                            %spike_times_vec_loop(indices_loop) = a_loop + (spike_times_vec_loop(indices_loop) - a_dash_loop)*(b_loop - a_loop)/(b_dash_loop - a_dash_loop);
                        else % (k == frame_end_loop) && (j == p.Num_FNoise_rep_ceil) % No end time is given for the last frame so have to specify differently
                            indices_loop = find((spike_times_vec_loop>...
                                trig_times_vec_par((j-1)*p_par.stim_frames + k))&&...
                                (spike_times_vec_loop<trig_times_vec_par((j-1)*p_par.stim_frames + k) + min_end_int));
                            a_loop      = trig_times_vec_par(k);
                            b_loop      = trig_times_vec_par(k+1);
                            a_dash_loop = trig_times_vec_par((j-1)*p_par.stim_frames + k);
                            b_dash_loop = trig_times_vec_par((j-1)*p_par.stim_frames + k) + min_end_int;
                            %spike_times_vec_loop(indices_loop) = a_loop + (spike_times_vec_loop(indices_loop) - a_dash_loop)*(b_loop - a_loop)/(b_dash_loop - a_dash_loop);
                        end
                        spike_times_vec_loop(indices_loop) = a_loop +...
                            (spike_times_vec_loop(indices_loop) - a_dash_loop)*...
                            (b_loop - a_loop)/(b_dash_loop - a_dash_loop); % PAR Mod 17,09,2020 (moved from inside if statement above)
                    end
                    
                end
                
                
                 % Map spike in interval after last stimulus window % PAR Mod 17,09,2020 (everything from here to the next 'else' is new)
                indices_loop = find((spike_times_vec_loop>trig_times_vec(p.Num_trigs) + stim_int)...
                                   &(spike_times_vec_loop<trig_times_vec(p.Num_trigs) + min_end_int));
                a_loop      = trig_times_vec(1);
                b_loop      = trig_times_vec(2);
                a_dash_loop = trig_times_vec(p.Num_trigs) + stim_int;
                b_dash_loop = trig_times_vec(p.Num_trigs) + min_end_int;
                spike_times_vec_loop(indices_loop) = a_loop +...
                    (spike_times_vec_loop(indices_loop) - a_dash_loop)*...
                    (b_loop - a_loop)/(b_dash_loop - a_dash_loop);
                
            else % If there are gaps between (any) frozen noise chuncks
                
                 if (nr_trig_last_repeat >= p.Num_STE_bins + 1) || (Num_Trig_Final_NoiseChunk == 0) % PAR Mod 17,09,2020 (if statement is new, but first for loop it contains was here originally)
                        for j = 1:p.Num_FNoise_rep_ceil-1 % remove spikes with stimulus windows that fall between the noise chunks
                            if p.Time_Choice == 1 % stimulus frames % PAR Mod 09,10,2020 (if statement is new, but first for loop it contains was here originally)
                                spike_times_vec_loop((spike_times_vec_loop>...
                                    (trig_times_vec(j*p.stim_frames) + min_end_int))&...
                                    (spike_times_vec_loop<trig_times_vec(j*p.stim_frames+p.Num_STE_bins+1))) = [];
                            else % p.Time_Choice == 2 % own time grid
                                spike_times_vec_loop((spike_times_vec_loop>...
                                    trig_times_vec(j*p.stim_frames) + min_end_int)&...
                                    (spike_times_vec_loop<trig_times_vec(j*p.stim_frames+1)+STE_int)) = [];
                            end
                        end
                elseif p.Num_FNoise_rep_ceil>2 % still need to remove spike between earlier chunks, priveded there were more than 2 originally
                     for j = 1:p.Num_FNoise_rep_ceil-2 % remove spikes with stimulus windows that fall between the noise chunks
                            if p.Time_Choice == 1 % stimulus frames
                                spike_times_vec_loop((spike_times_vec_loop>...
                                    (trig_times_vec(j*p.stim_frames) + min_end_int))&...
                                    (spike_times_vec_loop<trig_times_vec(j*p.stim_frames+p.Num_STE_bins+1))) = [];
                            else % p.Time_Choice == 2 % own time grid
                                spike_times_vec_loop((spike_times_vec_loop>...
                                    trig_times_vec(j*p.stim_frames) + min_end_int)&...
                                    (spike_times_vec_loop<trig_times_vec(j*p.stim_frames+1)+STE_int)) = [];
                            end
                      end
                end
                
                for j = 2:p_par.Num_FNoise_rep_ceil % map spikes in repeated stimulus chunks to original chunk
                    
                    if j==p_par.Num_FNoise_rep_ceil
                        frame_end_loop = p_par.Num_trigs - p_par.stim_frames*(p_par.Num_FNoise_rep_ceil-1);
                    else
                        frame_end_loop = p_par.stim_frames;
                    end
                    
                    for k = 1:frame_end_loop+1
                        if k<frame_end_loop % frame_end_loop<p.stim_frames
                            
                            indices_loop = find((spike_times_vec_loop>=trig_times_vec((j-1)*p_par.stim_frames + k)).*...
                                (spike_times_vec_loop<trig_times_vec((j-1)*p_par.stim_frames + k + 1)));
                            a_loop      = trig_times_vec(k);
                            b_loop      = trig_times_vec(k+1);
                            a_dash_loop = trig_times_vec((j-1)*p_par.stim_frames + k);
                            b_dash_loop = trig_times_vec((j-1)*p_par.stim_frames + k + 1);
                            %spike_times_vec_loop(indices_loop) = a_loop + (spike_times_vec_loop(indices_loop) - a_dash_loop)*(b_loop - a_loop)/(b_dash_loop - a_dash_loop);
                       elseif  k==frame_end_loop % No end time is given for the last frame so have to specify differently % PAR Mod 17,09,2020 (was 'else')
                            indices_loop = find((spike_times_vec_loop>trig_times_vec((j-1)*p.stim_frames + k))&...
                                (spike_times_vec_loop<trig_times_vec((j-1)*p.stim_frames + k) + stim_int)); % PAR Mod 17,09,2020 (& not &&)
                            a_loop      = trig_times_vec(k);
                            b_loop      = trig_times_vec(k) + stim_int; % PAR Mod 17,09,2020 (trig_times_vec(k) + stim_int not trig_times_vec(k+1))
                            a_dash_loop = trig_times_vec((j-1)*p_par.stim_frames + k);
                            b_dash_loop = trig_times_vec((j-1)*p_par.stim_frames + k) + min_end_int;
                            %spike_times_vec_loop(indices_loop) = a_loop + (spike_times_vec_loop(indices_loop) - a_dash_loop)*(b_loop - a_loop)/(b_dash_loop - a_dash_loop);
                        
                        else % k==frame_end_loop+1 % PAR Mod 17,09,2020 (everything between this else statement and the end of the if statement is new)
                            indices_loop = find((spike_times_vec_loop>trig_times_vec((j-1)*p.stim_frames + k-1) + stim_int)...
                                               &(spike_times_vec_loop<trig_times_vec((j-1)*p.stim_frames + k-1) + min_end_int));
                            if frame_end_loop == p.stim_frames
                                a_loop      = trig_times_vec(p.stim_frames) + stim_int;
                                b_loop      = trig_times_vec(p.stim_frames) + min_end_int;
                                a_dash_loop = trig_times_vec(j*p.stim_frames) + stim_int;
                                b_dash_loop = trig_times_vec(j*p.stim_frames) + min_end_int;
                            else %frame_end_loop < p.stim_frames
                                a_loop      = trig_times_vec(k);
                                b_loop      = trig_times_vec(k+1);
                                a_dash_loop = trig_times_vec((j-1)*p.stim_frames + k-1) + stim_int;
                                b_dash_loop = trig_times_vec((j-1)*p.stim_frames + k-1) + min_end_int;
                            end
                                                     
                        end    
                        spike_times_vec_loop(indices_loop) = a_loop +...
                            (spike_times_vec_loop(indices_loop) - a_dash_loop)*...
                            (b_loop - a_loop)/(b_dash_loop - a_dash_loop); % PAR Mod 17,09,2020 (moved from inside if statement above)               
                        
                    end
                    
                end
                
            end
            
        end
        
        %p.length_spike_times      = length(spike_times_vec_loop);
        length_spike_times_loop   = length(spike_times_vec_loop); % To allow parfor
        
        %RF_Ident{i} = RF_Ident_fn_v5(stimulus_arr,trig_times_vec_trunc,spike_times_vec_loop,length_spike_times_loop,p_par); % PAR Mod 27,08,2020 --> was 'RF_Ident_fn_v4' and 'trig_times_vec'
        RF_Ident(i).RF_results =...
            RF_Ident_fn_v8(stimulus_arr,trig_times_vec_trunc,spike_times_vec_loop,length_spike_times_loop,p); % PAR Mod 17,09,2020 (RF_Ident_fn_v5 -->  RF_Ident_fn_v6) %%% PAR Mod 27,08,2020 --> was 'RF_Ident_fn_v4' and 'trig_times_vec'
        RF_Ident(i).cell_idx = cell_indices_temp(i);
        
        RF_Ident(i).add_info = settings.settings.kernel_new;
        
        
               
        
        disp(i);
        
    end
    
    delete(gcp('nocreate'));
    
else %  Parpool == 2 % Parpool off
    
    for i = 1:True_Num_Cells
        
        spike_times_vec_loop      = spiketimestamps(:,i);
        
        test_for_NaN = find(all(isnan(spike_times_vec_loop),1));
        if test_for_NaN
            RF_Ident(i).RF_results = NaN;
            RF_Ident(i).cell_idx = cell_indices_temp(i);
            if i == 1&& bb == 1 
                RF_Ident(i).add_info = add_info.settings.kernel_new;
            end
            continue
        end
        
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
                                indices_loop = find((spike_times_vec_loop>=trig_times_vec((j-1)*p.stim_frames + k))&...
                                    (spike_times_vec_loop<trig_times_vec((j-1)*p.stim_frames + k + 1))); % PAR Mod 09,10,2020 (changed first inequality '>' to '>=')
                                a_loop      = trig_times_vec(k);
                                b_loop      = trig_times_vec(k+1);
                                a_dash_loop = trig_times_vec((j-1)*p.stim_frames + k);
                                b_dash_loop = trig_times_vec((j-1)*p.stim_frames + k + 1);
                                %spike_times_vec_loop(indices_loop) = a_loop + (spike_times_vec_loop(indices_loop) - a_dash_loop)*(b_loop - a_loop)/(b_dash_loop - a_dash_loop); % PAR Mod 17,09,2020 (moved outside of if statement)
                            else % (k == frame_end_loop) && (j == p.Num_FNoise_rep_ceil) % No end time is given for the last frame so have to specify differently
                                indices_loop = find((spike_times_vec_loop>=trig_times_vec((j-1)*p.stim_frames + k))&...
                                    (spike_times_vec_loop<trig_times_vec((j-1)*p.stim_frames + k) + stim_int)); % PAR Mod 09,10,2020 (changed first inequality '>' to '>=')
                                a_loop      = trig_times_vec(k);
                                b_loop      = trig_times_vec(k+1);
                                a_dash_loop = trig_times_vec((j-1)*p.stim_frames + k);
                                b_dash_loop = trig_times_vec((j-1)*p.stim_frames + k) + stim_int;
                                %spike_times_vec_loop(indices_loop) = a_loop + (spike_times_vec_loop(indices_loop) - a_dash_loop)*(b_loop - a_loop)/(b_dash_loop - a_dash_loop); % PAR Mod 17,09,2020 (moved outside of if statement)
                            end
                            spike_times_vec_loop(indices_loop) = a_loop + (spike_times_vec_loop(indices_loop) -...
                                a_dash_loop)*(b_loop - a_loop)/(b_dash_loop - a_dash_loop); % PAR Mod 17,09,2020 (moved from inside if statement above)
                        end
                    end
                    % Map spike in interval after last stimulus window % PAR Mod 17,09,2020 (everything from here to the next 'else' is new)
                    indices_loop = find((spike_times_vec_loop>=trig_times_vec(p.Num_trigs) + stim_int)&...
                        (spike_times_vec_loop<trig_times_vec(p.Num_trigs) + min_end_int)); % PAR Mod 09,10,2020 (changed first inequality '>' to '>=')
                    a_loop      = trig_times_vec(Num_Trig_Final_NoiseChunk+1); % (NB: If last chunck is full then Num_Trig_Final_NoiseChunk = 0)
                    b_loop      = trig_times_vec(Num_Trig_Final_NoiseChunk+2); % (NB: If last chunck is full then Num_Trig_Final_NoiseChunk = 0)
                    a_dash_loop = trig_times_vec(p.Num_trigs) + stim_int;
                    b_dash_loop = trig_times_vec(p.Num_trigs) + min_end_int;
                    spike_times_vec_loop(indices_loop) = a_loop +...
                        (spike_times_vec_loop(indices_loop) - a_dash_loop)*(b_loop - a_loop)/(b_dash_loop - a_dash_loop);
                    
                else % If there are gaps between (any) frozen noise chuncks
                    
                    if (nr_trig_last_repeat >= p.Num_STE_bins + 1) || (Num_Trig_Final_NoiseChunk == 0) % PAR Mod 17,09,2020 (if statement is new, but first for loop it contains was here originally)
                        for j = 1:p.Num_FNoise_rep_ceil-1 % remove spikes with stimulus windows that fall between the noise chunks
                            if p.Time_Choice == 1 % stimulus frames % PAR Mod 09,10,2020 (if statement is new, but first for loop it contains was here originally)
                                spike_times_vec_loop((spike_times_vec_loop>...
                                    (trig_times_vec(j*p.stim_frames) + min_end_int))&...
                                    (spike_times_vec_loop<trig_times_vec(j*p.stim_frames+p.Num_STE_bins+1))) = [];
                            else % p.Time_Choice == 2 % own time grid
                                spike_times_vec_loop((spike_times_vec_loop>...
                                    trig_times_vec(j*p.stim_frames) + min_end_int)&...
                                    (spike_times_vec_loop<trig_times_vec(j*p.stim_frames+1)+STE_int)) = [];
                            end
                        end
                    elseif p.Num_FNoise_rep_ceil>2 % still need to remove spike between earlier chunks, priveded there were more than 2 originally
                        for j = 1:p.Num_FNoise_rep_ceil-2 % remove spikes with stimulus windows that fall between the noise chunks
                            if p.Time_Choice == 1 % stimulus frames
                                spike_times_vec_loop((spike_times_vec_loop>...
                                    (trig_times_vec(j*p.stim_frames) + min_end_int))&...
                                    (spike_times_vec_loop<trig_times_vec(j*p.stim_frames+p.Num_STE_bins+1))) = [];
                            else % p.Time_Choice == 2 % own time grid
                                spike_times_vec_loop((spike_times_vec_loop>...
                                    trig_times_vec(j*p.stim_frames) + min_end_int)&...
                                    (spike_times_vec_loop<trig_times_vec(j*p.stim_frames+1)+STE_int)) = [];
                            end
                        end
                    end
                    
                    for j = 2:p.Num_FNoise_rep_ceil % map spikes in repeated stimulus chunks to original chunk
                        
                        if j==p.Num_FNoise_rep_ceil
                            frame_end_loop = p.Num_trigs - p.stim_frames*(p.Num_FNoise_rep_ceil-1);
                        else
                            frame_end_loop = p.stim_frames;
                        end
                        
                        for k = 1:frame_end_loop+1 % PAR Mod 17,09,2020 (was 1:frame_end_loop)
                            if k<frame_end_loop % frame_end_loop<p.stim_frames
                                indices_loop = find((spike_times_vec_loop>=...
                                    trig_times_vec((j-1)*p.stim_frames + k))&...
                                    (spike_times_vec_loop<trig_times_vec((j-1)*p.stim_frames + k + 1))); % PAR Mod 09,10,2020 (changed first inequality '>' to '>=') % PAR Mod 17,09,2020 (& not &&)
                                a_loop      = trig_times_vec(k);
                                b_loop      = trig_times_vec(k+1);
                                a_dash_loop = trig_times_vec((j-1)*p.stim_frames + k);
                                b_dash_loop = trig_times_vec((j-1)*p.stim_frames + k + 1);
                                %spike_times_vec_loop(indices_loop) = a_loop + (spike_times_vec_loop(indices_loop) - a_dash_loop)*(b_loop - a_loop)/(b_dash_loop - a_dash_loop); % PAR Mod 17,09,2020 (moved outside of if statement)
                            elseif  k==frame_end_loop % No end time is given for the last frame so have to specify differently % PAR Mod 17,09,2020 (was 'else')
                                indices_loop = find((spike_times_vec_loop>=...
                                    trig_times_vec((j-1)*p.stim_frames + k))&...
                                    (spike_times_vec_loop<trig_times_vec((j-1)*p.stim_frames + k) + stim_int)); % PAR Mod 09,10,2020 (changed first inequality '>' to '>=') % PAR Mod 17,09,2020 (& not &&)
                                a_loop      = trig_times_vec(k);
                                b_loop      = trig_times_vec(k) + stim_int; % PAR Mod 17,09,2020 (trig_times_vec(k) + stim_int not trig_times_vec(k+1))
                                a_dash_loop = trig_times_vec((j-1)*p.stim_frames + k);
                                b_dash_loop = trig_times_vec((j-1)*p.stim_frames + k) + stim_int;
                                %spike_times_vec_loop(indices_loop) = a_loop + (spike_times_vec_loop(indices_loop) - a_dash_loop)*(b_loop - a_loop)/(b_dash_loop - a_dash_loop); % PAR Mod 17,09,2020 (moved outside of if statement)
                            else % k==frame_end_loop+1 % PAR Mod 17,09,2020 (everything between this else statement and the end of the if statement is new)
                                indices_loop = find((spike_times_vec_loop>=...t
                                    rig_times_vec((j-1)*p.stim_frames + k-1) + stim_int)&...
                                    (spike_times_vec_loop<trig_times_vec((j-1)*p.stim_frames + k-1) + min_end_int)); % PAR Mod 09,10,2020 (changed first inequality '>' to '>=')
                                if frame_end_loop == p.stim_frames
                                    a_loop      = trig_times_vec(p.stim_frames) + stim_int;
                                    b_loop      = trig_times_vec(p.stim_frames) + min_end_int;
                                    a_dash_loop = trig_times_vec(j*p.stim_frames) + stim_int;
                                    b_dash_loop = trig_times_vec(j*p.stim_frames) + min_end_int;
                                else %frame_end_loop < p.stim_frames
                                    a_loop      = trig_times_vec(k);
                                    b_loop      = trig_times_vec(k+1);
                                    a_dash_loop = trig_times_vec((j-1)*p.stim_frames + k-1) + stim_int;
                                    b_dash_loop = trig_times_vec((j-1)*p.stim_frames + k-1) + min_end_int;
                                end
                            end
                            spike_times_vec_loop(indices_loop) = a_loop +...
                                (spike_times_vec_loop(indices_loop) - a_dash_loop)*...
                                (b_loop - a_loop)/(b_dash_loop - a_dash_loop); % PAR Mod 17,09,2020 (moved from inside if statement above)
                        end
                        
                    end
                    
                end
                
        end
        
        %p.length_spike_times      = length(spike_times_vec_loop);
        length_spike_times_loop   = length(spike_times_vec_loop); % To allow parfor
        
        RF_Ident(i).RF_results = RF_Ident_fn_v8(stimulus_arr,trig_times_vec_trunc,spike_times_vec_loop,length_spike_times_loop,p); % PAR Mod 17,09,2020 (RF_Ident_fn_v5 -->  RF_Ident_fn_v6) %%% PAR Mod 27,08,2020 --> was 'RF_Ident_fn_v4' and 'trig_times_vec'
        RF_Ident(i).cell_idx = cell_indices_temp(i);
        
        %Add info if first cell
        
        RF_Ident(i).add_info = add_info.settings.kernel_new;
        
        
           
        
    end
    
end

%Add some extra info
%Check if RF was found


        
        
        
%% Data organization and saving
%We have 2 x 2 conditions:
%STASD - SDThreshold, STASD - Confidence Interval
%SelfCovariance - SDThreshold, SelfCovaricance - Confidence Interval

%We need to save these independently, so it doesnt get overwritten when the
%repective other QC Type is choosen.


for rf = 1:length(RF_Ident)
   RF_overview(rf).p = p; 
    
end


if add_info.settings.kernel_new.SS_Std 
    
    %Saving the RF info
    subfoldername = ['RF_Ident_SS_Std_Thr_',num2str(add_info.settings.kernel_new.SS_STD)];
    save_folder = sf_organizer(stim_idx,savepath,'subfoldername',subfoldername);
    f = waitbar(0,'Saving files, file:');
    
    for rf = 1:length(RF_Ident)
        if isstruct(RF_Ident(rf).RF_results)
         waitbar(rf/length(RF_Ident),f,['Saving files, file: ',num2str(rf)])
         RF_overview(rf).file = strcat(save_folder,'\RF_Ident_',num2str(cell_indices_temp(rf)));
         fileID = fopen(strcat(save_folder,'\RF_Ident_',num2str(cell_indices_temp(rf)),'.bin'),'w');
         fwrite(fileID,hlp_serialize(RF_Ident(rf)));
         fclose(fileID);
        end
%          save_location = sf_organizer(stim_idx,savepath,'subfoldername',subfoldername,'variable',RF_Ident(rf),...
%                 'variable_name',['RF_Ident_',num2str(cell_indices_temp(rf))]);
    end
    try
        close(f)
    catch
        continue
    end
    [~] = sf_organizer(stim_idx,savepath,'subfoldername',subfoldername,'collect_files',true);
    
    %Check for Gaussian RF
    for i = 1:length(RF_Ident) 
        try
            if RF_Ident(i).RF_results.STASD_Gaus_FullRF_Num_pixels >= 1
                RF_overview(i).STASD_Gaus_RF = true;

            else

                RF_overview(i).STASD_Gaus_RF = false;
            end
        catch 
            RF_overview(i).STASD_Gaus_RF = false;
        end
    %Check for ASP RF
        try
            if RF_Ident(i).RF_results.STASD_ASP_FullRF_Num_pixels >= 1
                RF_overview(i).STASD_ASP_RF = true;

            else

                RF_overview(i).STASD_ASP_RF = false;
            end
        catch 
            RF_overview(i).STASD_ASP_RF = false;
        end
     %Check for Box RF
        try
            if RF_Ident(i).RF_results.STASD_Box_FullRF_Num_pixels >= 1
                RF_overview(i).STASD_Box_RF = true;

            else

                RF_overview(i).STASD_Box_RF = false;
            end
        catch 
            RF_overview(i).STASD_Box_RF = false;
        end
        
        RF_overview(i).cell_idx = cell_indices_temp(i);
        
        
    end
    
    if bb == 1
        save_location = sf_organizer(stim_idx,savepath,'subfoldername',subfoldername,'variable',RF_overview,...
            'variable_name','RF_overview');
       
    else
            %Create matfile and attach files to existing batch
        M = matfile(save_location,'Writable',true);
        M.RF_overview(1,batch_begins(bb):batch_ends(bb)) = RF_overview;


    end
    
   
    
elseif add_info.settings.kernel_new.SS_CI
    
    %Saving the RF info
    subfoldername = ['RF_Ident_SS_CI_Thr_',num2str(add_info.settings.kernel_new.SS_CI_Upper)];
    save_folder = sf_organizer(stim_idx,savepath,'subfoldername',subfoldername);
    f = waitbar(0,'Saving files, file:');
    
    for rf = 1:length(RF_Ident)
        if isstruct(RF_Ident(rf).RF_results)
         waitbar(rf/length(RF_Ident),f,['Saving files, file: ',num2str(rf)])
         RF_overview(rf).file = strcat(save_folder,'\RF_Ident_',num2str(cell_indices_temp(rf)));
         fileID = fopen(strcat(save_folder,'\RF_Ident_',num2str(cell_indices_temp(rf)),'.bin'),'w');
         fwrite(fileID,hlp_serialize(RF_Ident(rf)));
         fclose(fileID);
        end
%          save_location = sf_organizer(stim_idx,savepath,'subfoldername',subfoldername,'variable',RF_Ident(rf),...
%                 'variable_name',['RF_Ident_',num2str(cell_indices_temp(rf))]);
    end
    save_folder = sf_organizer(stim_idx,savepath,'subfoldername',subfoldername,'collect_files',true);
    
    
    %Check for Gaussian RF
    for i = 1:stx 
        try
            if RF_Ident(i).RF_results.STASD_Gaus_FullRF_Num_pixels >= 1
                RF_overview(i).STASD_Gaus_RF = true;

            else

                RF_overview(i).STASD_Gaus_RF = false;
            end
        catch 
            RF_overview(i).STASD_Gaus_RF = false;
        end
    %Check for ASP RF
        try
            if RF_Ident(i).RF_results.STASD_ASP_FullRF_Num_pixels >= 1
                RF_overview(i).STASD_ASP_RF = true;

            else

                RF_overview(i).STASD_ASP_RF = false;
            end
        catch 
            RF_overview(i).STASD_ASP_RF = false;
        end
     %Check for Box RF
        try
            if RF_Ident(i).RF_results.STASD_Box_FullRF_Num_pixels >= 1
                RF_overview(i).STASD_Box_RF = true;

            else

                RF_overview(i).STASD_Box_RF = false;
            end
        catch 
            RF_overview(i).STASD_Box_RF = false;
        end
        RF_overview(i).cell_idx = cell_indices_temp(i);
        
        
    end
        
    
    if bb == 1
        save_location = sf_organizer(stim_idx,savepath,'subfoldername',subfoldername,'variable',RF_overview,...
            'variable_name','RF_overview');
       
    else
            %Create matfile and attach files to existing batch
        M = matfile(save_location,'Writable',true);
        M.RF_overview(1,batch_begins(bb):batch_ends(bb)) = RF_overview;


    end

    
elseif add_info.settings.kernel_new.SC_Std
    
    subfoldername = ['RF_Ident_SC_Std_Thr_',num2str(add_info.settings.kernel_new.SC_STD)];
    save_folder = sf_organizer(stim_idx,savepath,'subfoldername',subfoldername);
    f = waitbar(0,'Saving files, file:');
    
    for rf = 1:length(RF_Ident)
        if isstruct(RF_Ident(rf).RF_results)
         waitbar(rf/length(RF_Ident),f,['Saving files, file: ',num2str(rf)])
         RF_overview(rf).file = strcat(save_folder,'\RF_Ident_',num2str(cell_indices_temp(rf)));
         fileID = fopen(strcat(save_folder,'\RF_Ident_',num2str(cell_indices_temp(rf)),'.bin'),'w');
         fwrite(fileID,hlp_serialize(RF_Ident(rf)));
         fclose(fileID);
        end
%          save_location = sf_organizer(stim_idx,savepath,'subfoldername',subfoldername,'variable',RF_Ident(rf),...
%                 'variable_name',['RF_Ident_',num2str(cell_indices_temp(rf))]);
    end
    save_folder = sf_organizer(stim_idx,savepath,'subfoldername',subfoldername,'collect_files',true);
    
    for i = 1:stx 
        %Check for Gaussian RF
        try
            if RF_Ident(i).RF_results.SC_Gaus_FullRF_Num_pixels >= 1
                RF_overview(i).SC_Gaus_RF = true;
            else
                RF_overview(i).SC_Gaus_RF = false;
            end
        catch 
            RF_overview(i).SC_Gaus_RF = false;
        end
        %Check for ASP RF
        try
            if RF_Ident(i).RF_results.SC_ASP_FullRF_Num_pixels >= 1
                RF_overview(i).SC_ASP_RF = true;
            else
                RF_overview(i).SC_ASP_RF = false;
            end
        catch 
            RF_overview(i).SC_ASP_RF = false;
        end
        %Check for Box RF
         try
            if RF_Ident(i).RF_results.SC_Box_FullRF_Num_pixels >= 1
                RF_overview(i).SC_Box_RF = true;
            else
                RF_overview(i).SC_Box_RF = false;
            end
        catch 
            RF_overview(i).SC_Box_RF = false;
        end
        RF_overview(i).cell_idx = cell_indices_temp(i);
    end
    
    if bb == 1
        save_location = sf_organizer(stim_idx,savepath,'subfoldername',subfoldername,'variable',RF_overview,...
            'variable_name','RF_overview');
       
    else
            %Create matfile and attach files to existing batch
        M = matfile(save_location,'Writable',true);
        M.RF_overview(1,batch_begins(bb):batch_ends(bb)) = RF_overview;


    end

elseif add_info.settings.kernel_new.SC_CI
    
    subfoldername = ['RF_Ident_SC_CI_Thr_',num2str(add_info.settings.kernel_new.SC_CI_Lower)];
    save_folder = sf_organizer(stim_idx,savepath,'subfoldername',subfoldername);
    f = waitbar(0,'Saving files, file:');
    
    for rf = 1:length(RF_Ident)
        if isstruct(RF_Ident(rf).RF_results)
         waitbar(rf/length(RF_Ident),f,['Saving files, file: ',num2str(rf)])
         RF_overview(rf).file = strcat(save_folder,'\RF_Ident_',num2str(cell_indices_temp(rf)));
         fileID = fopen(strcat(save_folder,'\RF_Ident_',num2str(cell_indices_temp(rf)),'.bin'),'w');
         fwrite(fileID,hlp_serialize(RF_Ident(rf)));
         fclose(fileID);
        end
%          save_location = sf_organizer(stim_idx,savepath,'subfoldername',subfoldername,'variable',RF_Ident(rf),...
%                 'variable_name',['RF_Ident_',num2str(cell_indices_temp(rf))]);
    end
    save_folder = sf_organizer(stim_idx,savepath,'subfoldername',subfoldername,'collect_files',true);
    
    for i = 1:stx 
        %Check for Gaussian RF
        try
            if RF_Ident(i).RF_results.SC_Gaus_FullRF_Num_pixels >= 1
                RF_overview(i).SC_Gaus_RF = true;
            else
                RF_overview(i).SC_Gaus_RF = false;
            end
        catch 
            RF_overview(i).SC_Gaus_RF = false;
        end
        %Check for ASP RF
        try
            if RF_Ident(i).RF_results.SC_ASP_FullRF_Num_pixels >= 1
                RF_overview(i).SC_ASP_RF = true;
            else
                RF_overview(i).SC_ASP_RF = false;
            end
        catch 
            RF_overview(i).SC_ASP_RF = false;
        end
        %Check for Box RF
         try
            if RF_Ident(i).RF_results.SC_Box_FullRF_Num_pixels >= 1
                RF_overview(i).SC_Box_RF = true;
            else
                RF_overview(i).SC_Box_RF = false;
            end
        catch 
            RF_overview(i).SC_Box_RF = false;
        end
        RF_overview(i).cell_idx = cell_indices_temp(i);
    end
    
    if bb == 1
        save_location = sf_organizer(stim_idx,savepath,'subfoldername',subfoldername,'variable',RF_overview,...
            'variable_name','RF_overview');
       
    else
            %Create matfile and attach files to existing batch
        M = matfile(save_location,'Writable',true);
        M.RF_overview(1,batch_begins(bb):batch_ends(bb)) = RF_overview;


    end
end
   

close(f)

clear RF_Ident
clear RF_overview
end
end
toc;



out = 1;
end





%% Save Data
%save('Save_Name.mat','-v7.3')

%save data_Test_#;
