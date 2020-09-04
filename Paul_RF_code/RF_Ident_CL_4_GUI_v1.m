%% RF Indentification Command Line 4 - GUI Ready
% 04,06,2020 Onwards
% modified by marvin seifert to make work with MEA analysis GUI
% Works with RF_Ident_fn_v4.m
function out = RF_Ident_CL_4_GUI_v1 (savepath, add_info)



%% Get GUI input
kernel_settings = add_info.settings.kernel_new;
Parpool = add_info.settings.parpro;
Num_Cores = 4;

%% Choose whether to use parallel processing and number of cores
% Parpool   = 1; % 1: Yes, 2: No.
% Num_Cores = 4;

%mars: Will add this later 

%% Make RF Identification Choices

% Choose RF Identification Method
% If multiple options are chosen
% p.RF_Ident_Meth_vec = [a,b,c,d];
% a: STA-SD method;
% b: LC method (Local Covariance);
% c: MI method (Mutual Information).
% d: SC method (Self Covariance).
p.RF_Ident_Meth_vec = [kernel_settings.SS,0,0,kernel_settings.CI];

% Choose RF Quality Control (QC) and Thresholds
if p.RF_Ident_Meth_vec(1) == 1
    
    % Choose STA-SD QC type
    % 1. SD threshold;
    % 2. CI (Confidence Interval).
    p.STA_SD_QC_Type = 1;
    
    if p.STA_SD_QC_Type == 1
        % Choose SD Threshold
        p.STA_SD_Thresh = kernel_settings.SS_Std; % MEA Data 1 (50px) -- 5 seems good. 
    elseif p.STA_SD_QC_Type == 2
        % Choose Upper Percentile
        p.STA_SD_CI_Per = kernel_settings.QC_CI_Upper; % default: 1 (%)
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
        p.SC_Thresh = kernel_settings.QC_Std; % MEA Data 1 (50px): 7 seems good.
    elseif p.SC_QC_Type == 2
        % Choose Upper Percentile
        p.SC_CI_Per = kernel_settings.QC_CI_Upper; % default: 1 (%)
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
p.Plot_Choice = 1; %mars: Should this be a tickbox in the GUI?

%% Make Data Related Choices

%%% Choose MEA cell(s)
% Choose all cells or a subset
% 1. All cells;
% 2. Subset of cells.
Cell_Choice = 1; %mars: Not sure we need this option
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
p.Num_STE_bins = kernel_settings.Num_STE_bins; % Default: 5 (Tom data, Time_Choice 1) or 10 (Tom data,Time_Choice 2); 20 (Marvin data)


%% Load Data
%mars: I added a new way to load data. It depends on which stimulus the
%user selects in the GUI. The spiketimestamps get loaded in batches, to
%prevent the memory from overfloating. It calls functions which come with
%the GUI and are used in other scripts of the GUI as well. 

%% mars: Memory handling 
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

%% mars: Specifiy trigger channel 
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
    
%% Trigger channel for each stimulus    
noise_begin = stim_begin(nn)*SamplingFrequency;
noise_end = stim_end(nn)*SamplingFrequency;

trigCh_vec(1,:) = Ch01_02(noise_begin:noise_end); %This cuts the
%trigger trace at stimulus begin and end

min_trigCh_vec     = min(trigCh_vec);
max_trigCh_vec     = max(trigCh_vec);
trigThreshFac      = 0.05; % 0.1
trigHigh_vec       = double(trigCh_vec > min_trigCh_vec + trigThreshFac*(max_trigCh_vec-min_trigCh_vec));
[~,trig_index_vec_temp] = findpeaks(trigHigh_vec);% 'MinPeakProminence',#,'MinPeakDistance',#
trig_index_vec = zeros(1,size(trig_index_vec_temp,2)+2);
trig_index_vec(1,1) = 1;
trig_index_vec(1,2:end-1) = trig_index_vec_temp;
trig_index_vec(1,end) = size(trigHigh_vec,2);
sampling_freq      = SamplingFrequency;
sampling_int       = 1/sampling_freq;

trig_times_vec_temp = (trig_index_vec-1)*sampling_int;
% I take it that the rising edge starts at the timepoint directly
% before the jump upwards - hence '(trig_index_vec-1)' rather than 'trig_index_vec'
trig_times_vec      = trig_times_vec_temp - trig_times_vec_temp(1); % set first frame time to 0 s.
nr_noise_frames = length(trig_times_vec);

%% Noise sequence for each stimulus 
%Noise sequence is either loaded from HDF5 file or text file

noise_file = add_info.settings.location.noise;
%Check for file ending 
dot_ending = strfind(noise_file,'.');
file_ending = noise_file(dot_ending+1:end);

if strcmp(file_ending,'h5')
stimulus_arr = load_noise_from_hdf5(string(noise_file),true,1,nr_noise_frames);
nr_boxes = sqrt(size(stimulus_arr,2));
nr_colours = size(stimulus_arr,3);
stimulus_arr = reshape(stimulus_arr,[nr_boxes,nr_boxes,nr_noise_frames,nr_colours]);

elseif strcmp(file_ending,'txt') %Old version of getting the sequence
%Get folder location
folder_ending = max(strfind(noise_file,'\'));
folder_dir = noise_file(1:folder_ending);

stimulus_arr(:,:,1) = importdata([folder_dir,'Red_Noise.txt'],',')';
stimulus_arr(:,:,2) = importdata([folder_dir,'Green_Noise.txt'],',')';
stimulus_arr(:,:,3) = importdata([folder_dir,'Blue_Noise.txt'],',')';
stimulus_arr(:,:,4) = importdata([folder_dir,'UV_Noise.txt'],',')';
stimulus_arr(:,1,:) = []; % stimulus_arr(1,:,:) before transposed above
stimulus_arr        = reshape(stimulus_arr,[40,40,6000,4]);
nr_colours = 4;
end

p.Spectral_Dim = nr_colours;

if size(stimulus_arr,3) ~= nr_noise_frames %mars This will not work if we 
    %only load a subset of frames depending on how many trigger signals we
    %find in the trigger channel (kinda circular in that case)
    disp('NB: There is a mismatch between the number of stimulus frames and the number of triggers!');
end

% if Cell_Choice == 1 % All cells
%     spike_times_mat    = spiketimestamps;
%     Cell_Choice_vec    = 1:1:size(spiketimestamps,2);
%     True_cell_index_vec = Cell_Choice_vec(any(~isnan(spike_times_mat)));
% else % Cell_Choice == 2 % Subset of cells
%     spike_times_mat    = spiketimestamps(:,Cell_Choice_vec);
%     % Remove NaN columns (non-responsive cells)
%     True_cell_index_vec = Cell_Choice_vec(any(~isnan(spike_times_mat))); % Record indices of remaining cells
% end



%% mars Loop over batches
%Load supset of cells into memory
h = waitbar(0,'Batch nr: ');
for bb = 1:batch_nr
waitbar(bb/batch_nr,h,[num2str(bb),' of ', batch_nr]);
spike_times_mat = load_spiketimestamps_app (savepath, add_info,...
    'cell_subset',[batch_begins(bb),batch_ends(bb)]);
Cell_Choice_vec    = 1:1:size(spike_times_mat,2);
True_cell_index_vec = Cell_Choice_vec(any(~isnan(spike_times_mat)));

%cell_indices_temp = cell_indices(batch_begins(bb):batch_ends(bb)); %Load indices of the cells

spike_times_mat(:,~any(~isnan(spike_times_mat))) = [];
True_Num_Cells     = size(spike_times_mat,2);

% trigCh_vec         = Ch_new.trigger_ch;
% Check that triggers are correctly identified.
% figure; plot(trigCh_vec); hold on;
% plot(trig_index_vec,4085*ones(1,length(trig_index_vec)),'ro');




%%% Set Parameters

% Stimulus ppties
stim_freq = kernel_settings.Hz;          % Hz %mars: This has to be not hard coded
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
first_trig_time = trig_times_vec(1);
last_trig_time  = trig_times_vec(end); % trig_times_vec(end) + stim_int

p.stim_rows        = size(stimulus_arr,1);
p.stim_columns     = size(stimulus_arr,2);
p.stim_pixels      = p.stim_rows*p.stim_columns;
stim_frames        = size(stimulus_arr,3);

Num_stim_spat_dim  = p.stim_rows*p.stim_columns;

Num_STE_Stixels = Num_stim_spat_dim*p.Num_STE_bins;

if p.STA_Choice  == 2 % subtract mean raw stim
    if p.Time_Choice == 1 || p.Mean_Stim_Choice  == 1 % stimulus frames || Calc from stim frames
        p.Num_Raw_Stim = stim_frames - p.Num_STE_bins + 1;
        % Total number of raw stimuli displayed
    elseif p.Time_Choice == 2                       % define own time grid
        p.Num_Raw_Stim = ceil((last_trig_time - first_trig_time + stim_int - STE_int)/(STE_int/p.Num_STE_bins)) + 1;%10*ceil((last_trig_time - first_trig_time + stim_int - STE_int)/(STE_int/p.Num_STE_bins)) + 1
        % This is the total time over which stimuli are played 'last_trig_time-first_trig_time+stim_int'
        % minus the length of an STE 'STE_int' divided by the time between time
        % samples 'STE_int/p.Num_STE_bins' rounded upwards to the nearest integer, can times 10 to try to
        % make sure don't miss any stimulus combinations and plus 1 to count
        % the first stimulus before shifting the window.
        % 111.208542 seconds for model 1 without *10 so do this o/w takes too
        % long.
        p.Num_Raw_Stim_t_vec = linspace(first_trig_time + STE_int,last_trig_time + stim_int,p.Num_Raw_Stim);
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
            p.mean_raw_stim_arr(:,:,:,i) = mean_raw_stim_SpaceTime_fn(stimulus_arr(:,:,:,i),trig_times_vec,p);
        end
    end
end




if Parpool % Parpool on
    try
    parpool(Num_Cores);
    catch
    delete(gcp('nocreate'));
    parpool(Num_Cores);
    end
    
    parfor i = 1:True_Num_Cells % for/parfor
        i
        spike_times_vec_loop      = spike_times_mat(:,i);
        spike_times_vec_loop(isnan(spike_times_vec_loop)) = []; % spiketimestamps(~isnan(spiketimestamps(:,Cell_Choice)),Cell_Choice);
        spike_times_vec_loop      = spike_times_vec_loop(spike_times_vec_loop>(first_trig_time + min_start_int));
        spike_times_vec_loop      = spike_times_vec_loop(spike_times_vec_loop<(last_trig_time + min_end_int));
        %p.length_spike_times      = length(spike_times_vec_loop);
        length_spike_times_loop   = length(spike_times_vec_loop); % To allow parfor
        
        RF_Ident(i) = RF_Ident_fn_v4(stimulus_arr,trig_times_vec,spike_times_vec_loop,length_spike_times_loop,p);
        %mars: Does not work at the moment because some variables are not
        %preallocated and therefore nonexistant if no value is assigned...
        
        disp(i);
        
    end
    
    delete(gcp('nocreate'));
    
else %  Parpool == 2 % Parpool off
    
    for i = [62,124,123,186,248,247]%1:True_Num_Cells
        
        spike_times_vec_loop      = spike_times_mat(:,i);
        spike_times_vec_loop(isnan(spike_times_vec_loop)) = []; % spiketimestamps(~isnan(spiketimestamps(:,Cell_Choice)),Cell_Choice);
        spike_times_vec_loop      = spike_times_vec_loop(spike_times_vec_loop>(first_trig_time + min_start_int));
        spike_times_vec_loop      = spike_times_vec_loop(spike_times_vec_loop<(last_trig_time + min_end_int));
        %p.length_spike_times      = length(spike_times_vec_loop);
        length_spike_times_loop   = length(spike_times_vec_loop); % To allow parfor
        
        RF_Ident(i) = RF_Ident_fn_v4(stimulus_arr,trig_times_vec,spike_times_vec_loop,length_spike_times_loop,p);
        %mars: Does not work at the moment because some variables are not
        %preallocated and therefore nonexistant if no value is assigned...
        
        disp(i);
        
    end
    
end
end
end
end

%% Save Data
%save('Save_Name.mat','-v7.3')
