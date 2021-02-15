function [stim_begin, stim_end, filename, Ch] = findstimtimes 
%This function find the stimulus begins and ends in the trigger channel
%recorded with the Biocam MEA. Stimulus begins and ends are identified by
%first searching for all peaks (trigger signals) in the trigger channel and
%than identifying the trigger signals which have unique intervals (The idea
%is that trigger intervals remain constant if they are part of a repetitive
%stimulus but vary when the user starts a new stimulus) 
%Trigger channel can be loaded as .brw file or as .mat file. @MSeifert 2020

%% Load the file
%Load the stimulus file (Channel 1,2)
[stimulus_file,stimulus_path] = uigetfile({'*.mat';'*.brw'},'Select stimulus file');
cd (stimulus_path);
format = get_fileformat(stimulus_file);
filename = [stimulus_path,stimulus_file];
 
 %Either .mat files or .brw files can be loaded
 if strcmp('.mat',format)
    Ch = load(stimulus_file);
  
 elseif strcmp('.brw',format) %This loads the data directly from brw file
   Ch.Ch01_02 = h5read(stimulus_file,'/3BData/Raw'); 
   Ch.Ch01_02 = Ch.Ch01_02'; %Tranformation is necessary because matfile and brw file
   %use different array dimensions
   Ch.SamplingFrequency = h5read(stimulus_file,'/3BRecInfo/3BRecVars/SamplingRate');
     
 end
 
 stim_Ch = Ch.Ch01_02;
 
 %if cuda driver exists save variable into gpu Array (makes it easier to plot)
 if gpuDeviceCount
     stim_Ch = gpuArray(stim_Ch);
 end
 
 %plot(stim_Ch)
 
 %% Find trigger signals
 %Estimate the right threshold for the trigger detection 
 [value,~] = findpeaks(gather(double(stim_Ch)),'MinPeakProminence',1000,'MinPeakDistance',178);
 
 detection_threshold = mean(value)*0.8;
 %Check if the threshold is safe
 threshold_std = std(value);
 if threshold_std > 1
     error(['couldnt find trigger size, check if trigger signals are of equal'...
         ' power'])
 end
 
 %transfer the trigger channel into a logical matrix so that all events
 %bigger than the detection threshold are true, all other are false
 stim_Ch_norm = stim_Ch > detection_threshold;
 
 % identify the peaks in the trigger channel
 [~,locs] = findpeaks(gather(double(stim_Ch_norm)),'MinPeakProminence',1,'MinPeakDistance',178);
 %find the beginning and end of stimuli 
 locs_diff = diff(locs);
 
 %Check if time difference is lager than a second
 locs_diff_true = locs_diff>Ch.SamplingFrequency;
 position_true = find(locs_diff_true);
 

 %Check for unique trigger times which are different from all other non
 %unique trigger times
 [unique_locs, position] = unique(locs_diff(locs_diff_true));
 position = position_true(position);
%% Identify unique trigger signals 
 true_unique = logical(unique_test(unique_locs,0.01));
 
 %position = position';
 locs_idx = position(1,true_unique);
 locs_idx = sort(locs_idx);
 stim_end = locs(1,locs_idx);
 
 %Find the last end of the last stimulus (which is not found yet)
 last_peak = locs(1,end);
 Ch_last = stim_Ch(1,last_peak:end);
 Ch_last_norm = Ch_last > detection_threshold;
 last_peak_end = gather(find(Ch_last_norm==0, 1, 'first'));
 stim_end(1,end+1) = last_peak+last_peak_end;
 
 
 %find the beginning of the stimuli
 %locs_idx is the index of the last peak in each stimulus, thus, locs_idx+1
 %is the beginning index of the previous stimulus (except for the first
 %stimulus)
 locs_idx_b = locs_idx + 1;
 stim_begin = locs(1,locs_idx_b);
 
 %The begin of the first stimulus is the first peak detected
 stim_begin_temp = zeros(1,length(stim_begin)+1);
 stim_begin_temp(1,2:end) = stim_begin;
 stim_begin_temp(1,1) = locs(1,1);
 stim_begin = stim_begin_temp;
 
 
 %% Add the overhanging stimulus end
 %the last trigger is the end of the stimulus, but the time window has be
 %be extended beyond that so that we get the responding spikes.
 
 for ii = 1:length(stim_begin)
    try
    [~,locs_stim] = findpeaks(gather(double(stim_Ch_norm(stim_begin(ii):stim_end(ii))),'MinPeakProminence',1,'MinPeakDistance',178));
    if isempty(locs_stim)
        continue
    else
    stim_diff = diff(locs_stim);
    trig_max = max(stim_diff);
    stim_end(ii) = stim_end(ii)+trig_max; %Problem if the trigger channel is shorter...
    end
    catch
        continue
        %Error can be thrown when the stim begin and stim end are too close
        %together
    end
         
 end
 
 %Check if last stim end is outside the trigger channel size, if thats the
 %case the trigger channel gets extended with zeros;
 
 
 if stim_end(end) > size(stim_Ch,2)

    stim_Ch(end:stim_end(end)) = 0;
    Ch.Ch01_02 = stim_Ch;
    
 end
 
 
    
 
 
 
 
 
    

 
 
 
 
 
 
 
 
 
 


end