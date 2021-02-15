function out = extract_dataset (savepath,add_info)
%This script extracts the spiketimes and trigger times for a given
%stimulus, or a list of stimuli and returns a new matfile with the
%respective information.
%@MSeifert 2020

%% Error check
%Check if stim_begin and stim_end exists
if isempty(add_info.stim_idx)
    error('No stimulus selected, select stimulus first')
end

%Check if matfile exists
if exist(savepath,'file') == 0
    error("Savepath refers to nonexisting file, check file location")
end

%% Main
stim_idx = add_info.stim_idx;
load(savepath,'stim_begin','stim_end');
stim_begin = stim_begin(stim_idx);
stim_end = stim_end(stim_idx);

nr_stimuli = numel(stim_begin);


%Ask if all cells or only subset of cells shall be exported
prompt = {'Enter index of first cell:','Enter index of last cell:'};
dlgtitle = 'Export subset of cells? If no click cancel';
dims = [1 35];
definput = {'1','100'};
cell_subset = inputdlg(prompt,dlgtitle,dims,definput);
if isempty(cell_subset)
    cell_subset = 0;
else
    cell_subset = str2double(cell_subset');
end





for kk = 1:nr_stimuli
    
   %Get the stim times for the stimulus
    st_b = stim_begin(kk);
    st_e = stim_end(kk);
    
    data.stim_begin = st_b;
    data.stim_end = st_e;
%Load the spiketimestamps
    if cell_subset == 0
        spiketimestamps = load_spiketimestamps_app (savepath, data);
    else
        spiketimestamps = load_spiketimestamps_app (savepath, data, 'cell_subset',cell_subset);
    end
    
%This resets the spiketimestamps of each stimulus back to the stimulus
%begin so that spiketimestamps start relative to the beginning.
spiketimestamps = spiketimestamps-st_b;


%Next load the trigger times, for this we need to load the Ch variable

S = load(savepath,'Ch');
Ch = S.Ch;
clear S

trigger_begin = st_b*Ch.SamplingFrequency-100;
trigger_end = st_e * Ch.SamplingFrequency+100;

Ch_new.trigger_ch = Ch.Ch01_02(trigger_begin:trigger_end);
Ch_new.SamplingFrequency = Ch.SamplingFrequency;

%Ask if also analysis folders belonging to the stimulus shall be extracted

%To be finished


%% Save the data in a mat file
title = ['Save data for stimulus ',num2str(stim_idx(kk))]; 
filter = {'*.mat';'*.m';'*.slx';'*.*'}; 
[file, path] = uiputfile(filter,title);
cd (path)

save(file,'spiketimestamps','Ch_new');
end




    
    
    
        
    
    
    
    
    
    
   out = 1; 
end
















