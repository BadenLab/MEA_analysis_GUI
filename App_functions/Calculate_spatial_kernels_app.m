function out = Calculate_spatial_kernels_app (savepath, add_info)

%Check how many clusters we have
M = matfile(savepath,'Writable',false);
spikesize = numel(M.spiketimestamps);
stx_before = size(M.spiketimestamps,2);

%One double entry is 8 bytes, we dont want to use more than one fourth of the
%available memory
memory_available = memory;
memory_available = memory_available.MaxPossibleArrayBytes/4;
array_size_available = memory_available/8;

%If loading all spiketimestamps requires more than half the memory
%available, batch mode is activated.

    %Calculate batch size
    batch_nr = ceil(spikesize/array_size_available);
    batch_stx = ceil(stx_before/batch_nr);
    batch_begins = (1:batch_stx:stx_before);
    batch_ends = batch_begins(2:end)-1;
    batch_ends(end+1) = stx_before;

    


%Extract the important information from input
stim_idx = add_info.stim_idx;
stim_begin = add_info.stim_begin;
stim_end = add_info.stim_end;
ch = load(savepath,'-mat','Ch');
ch = ch.Ch;
Ch01_02 = ch.Ch01_02;


loadpath = load(savepath,'pathname'); 
pathname = loadpath.pathname;

SamplingFrequency = ch.SamplingFrequency;

%Load cell indices
S = load(savepath,'-mat','cell_indices');
cell_indices = S.cell_indices;

%Check how many stimulus repeas are in the file
noise_repeats = numel(stim_begin);

for nn = 1:noise_repeats
%Define path for saving the output
% stim_idx = add_info.stim_idx(nn);
% stimulus_str = string(['Stimulus_',num2str(stim_idx),'_Kernels']);
% cd(pathname)
% mkdir (stimulus_str)
% savepath_Kernel = strcat(pathname,stimulus_str,'\');    
    
noise_begin = stim_begin(nn)*SamplingFrequency;
noise_end = stim_end(nn)*SamplingFrequency;
noise_begin_s = stim_begin(nn);
noise_end_s = stim_end(nn);

noise_trigger_trace(1,:) = Ch01_02(noise_begin:noise_end);
noise_trigger_trace(2,:) = (noise_begin:1:noise_end);
noise_diff_trigger = diff(noise_trigger_trace(1,:));
noise_trigger_norm = noise_trigger_trace(1,:) > 2500;
%find the peaks in the differences (which correspond to the beginning of a
%trigger event
[~,locs] = findpeaks(gather(double(noise_trigger_norm)),'MinPeakProminence',1,'MinPeakDistance',178);
%add the time from the stimulus begin
locs_realtime = locs + noise_begin;
%plot trigger times in histogram
locs_diff = diff(locs_realtime)/SamplingFrequency;
locs_diff_norm = ceil05(locs_diff,0.05);
locs_s = zeros(1,length(locs(1,:)));
for i = 1:length(locs(1,:))
    if i == 1
        locs_s(i) = noise_begin_s;
    else
        locs_s(i) = locs_s(i-1)+locs_diff_norm(i-1);
    end
end
figure
hold on
histogram(locs_diff)
title ('Noise interval histogram')
xlabel('time in [s]')
ylabel('counts')
% noise_end_s - noise_begin_s;
% sum(locs_diff) + mean(locs_diff);
locs_s = locs_realtime / SamplingFrequency;
%add the start time to the array
locs_s_temp = NaN(1,length(locs_s)+1);
locs_s_temp(1,2:end) = locs_s;
locs_s_temp(1,1) = noise_begin_s;
locs_s = locs_s_temp;
clear locs_s_temp

%% Load stimulus files

colour_noise_path = 'D:\Stimuli_Data\Lightcrafter\colour_noise_spatial';
cd (colour_noise_path)
Red_name = 'Red_Noise';
Green_name = 'Green_Noise';
Blue_name = 'Blue_Noise';
UV_name = 'UV_Noise';

Colour_noise(:,:,4) = importdata([colour_noise_path, '\', Red_name,'.txt'],',');
Colour_noise(:,:,3) = importdata([colour_noise_path, '\', Green_name,'.txt'],',');
Colour_noise(:,:,2) = importdata([colour_noise_path, '\', Blue_name,'.txt'],',');
Colour_noise(:,:,1) = importdata([colour_noise_path, '\', UV_name,'.txt'],',');

Colour_noise(1,:,:) = [];


h = waitbar(0,'Batch nr: ');
% for bb = 1:batch_nr
Kernel_location = strings(stx_before,2);
for bb = 1:batch_nr
    waitbar(bb/batch_nr,h,[num2str(bb),' of ', batch_nr]);
%Initiate the batch    
spiketimestamps = load_spiketimestamps_app (savepath, add_info,...
    'cell_subset',[batch_begins(bb),batch_ends(bb)]);
cell_indices_temp = cell_indices(batch_begins(bb):batch_ends(bb));

%Quality Check

right_indices = quality_check_spikes(spiketimestamps,10);
spiketimestamps = spiketimestamps(:,right_indices);
stx = size(spiketimestamps,2);
sty = size(spiketimestamps,1);  

cell_indices_temp = cell_indices_temp(1,right_indices);


    
%% Extracting the right spiketimes
real_spikes = NaN(sty,stx);
for ss = 1:stx
    spikes = full(spiketimestamps(:,ss));
    real_spikes_idx = logical((spikes>noise_begin_s).*(spikes<noise_end_s));
    spikes_l = sum(real_spikes_idx);
    real_spikes(1:spikes_l,ss) = spikes(real_spikes_idx);
end
%cut the NaN values
[~, real_spikes] = lastNaN(real_spikes);

stx_rs = length(real_spikes(1,:));
sty_rs = length(real_spikes(:,1));

%% calculate the kernels

%initilizing parfor loop
%parfor rs = 1:stx_rs
kbinsize = 0.001;
kernel_window = 0.5;
bins = kernel_window/kbinsize;
%create matfile to save the kernels
%first create placeholde

cd (pathname)
% x = 'kernels';
% save('Kernels.mat','x');
% Kernels = matfile('Kernels.mat','Writable',true);

%Check if the subfolder for the stimulus exists
out = sf_organizer(stim_idx,savepath,'subfoldername',"Kernel");


waitmessage = parfor_wait(stx_rs,'Waitbar',true);
Kernel_location_temp = strings(stx_rs,1);

parfor rs = 1:stx_rs
    waitmessage.Send;
    rs %= Best_cells(rs1);
    %cd (savepath_Kernel)
    %Kernels = matfile(sprintf('output%d.mat',rs),'Writable',true)
    Colour_noise_l = Colour_noise;
    %name = num2str(rs);
    %the following variables are defined new for use in the parfoor loop
    spikes = real_spikes(:,rs);
    locs_sl = locs_s;
    nnz_spikes = nnz(~isnan(spikes));
    kernel_window = 0.6;
    pre_window = 0.5*kernel_window;
    post_window = pre_window;
    stim_idx_temp = stim_idx;
    
     
    kernels_temp = zeros(bins,length(Colour_noise_l(1,:,1)),4);
    for ii = 1:nnz_spikes
        
        
        spiketime = spikes(ii,1);
        %%find the right trigger time
        delta_time = locs_sl - spiketime;
        %look for the largest negative value, a is the difference between
        %the spike and the last trigger event
        [~,delta_idx] = max(delta_time(delta_time<=0));
        %calculate the range for averaging the stimulus
        delta_time_1 = locs_sl - (spiketime-pre_window);
        [a_1,delta_idx_1] = max(delta_time_1(delta_time_1<=0));
        
        delta_time_2 = locs_sl - (spiketime + post_window);
        [a_2,delta_idx_2] = max(delta_time_2(delta_time_2<=0));
      
        %delta_1 indicates at which sequence we need to start, delta 2
        %indicates where to stop. Next we need to load the respective
        %sequence
        if delta_idx <= 100
            %Skip the first 5 trigger
        elseif delta_idx_2 >= length(locs_sl(1,:))
            break
            %skip spikes outside of the sequence
        else
        sequence = Colour_noise_l(delta_idx_1:delta_idx_2,:,:);
        sequence_x = length(sequence(:,1,1));
        sequence_y = length(sequence(1,:,1));
        sequence_z = length(sequence(1,1,:));
        locs_window = locs_sl(1,(delta_idx_1:delta_idx_2+1));
        
        
        
%         noise_freq = locs_window/(delta_idx_2 - delta_idx_1+1);
%          upsample_factor = 50; %(noise_freq/kbinsize);
         %sequence_upsample = zeros(sequence_x,sequence_y, sequence_z);
         %Check the trigger channel
         test_diff = max(diff(locs_window));
         if test_diff > 1
             continue
         end
         
         
        [sequence_upsample, trigger_times] = sequence_complete(sequence,locs_window);
         sequence_timing = locs_window(1) + (0:length(sequence_upsample(:,1,1))-1)*kbinsize;
        
            
        
        delta_time_3 = sequence_timing - spiketime;
        [a_3,delta_idx_3] = max(delta_time_3(delta_time_3<=0));
        %delta_idx_3 is the bin closest to the actual spike
        
        
%         sequence_upsample = repelem(sequence,50,1,1);
%         a = 1;
%         for kk = 1:length(sequence(:,1,1))
%             for tt = 1:upsample_factor
%                 sequence_upsample(a,:,:) = sequence(kk,:,:);
%                 a = a+1;
%             end
%         end
        
        %Next we have to identify the right range for the time window we
        % have set
        %a1 is the difference between the first trigger and the beginning
        %of the time window, a2 is the time from the last trigger to the
        %end of the time window
                       
%         if bins_delete <=0
%             break
%         end
                  
        sequence_upsample(delta_idx_3+bins/2:end,:,:) = [];
        bins_o = length(sequence_upsample(:,1,1))-bins;
        
        try
            sequence_upsample = sequence_upsample(bins_o+1:end,:,:);
        catch
            sequence_upsample = zeros(bins,length(sequence_upsample(1,:,1)),4);
        
        end
        
%         if length(sequence_upsample(:,1,1)) > bins
%         delete_a1 = floor((-1*a_1)/kbinsize);
%         delete_a2 = ceil(upsample_factor - (-1*a_2)/kbinsize);
%         if delete_a2 < 0 
%             delete_a2 = 0;
%         end
%         sequence_upsample =  sequence_upsample(1:end-delete_a2,:,:);
%         sequence_upsample(1:delete_a1,:,:) = [];
%         end
%         if length(sequence_upsample(:,1,1)) < bins
%             diff_seq = bins - length(sequence_upsample(:,1,1));
%             sequence_upsample(end:end+diff_seq,1,1) = zeros;
%         end
%         if length(sequence_upsample(:,1,1))>bins
%            diff_seq = length(sequence_upsample(:,1,1));
%            sequence_upsample(bins+1:diff_seq,:,:) = [];
%         end
        kernels_temp = kernels_temp + sequence_upsample;
       
        end
        
    end
    %kernels_save = zeros(1,600,1600,4);
    %kernels_save(1,:,:,:) = kernels_temp;
    %Create save name
    %file_name = ['Kernels_', num2str(rs),'.mat'];
%     file_name = 'test';
%     partmatfile(sprintf('output%d.mat',rs),kernels_temp);
%     Kernel_location(rs,1) = strcat(savepath_Kernel,sprintf('output%d.mat',rs));
    %Create name for the kernel file
    Kernel_name = ['Kernel_Cell_',num2str(cell_indices_temp(rs))];
    Kernel_location_temp(rs,1) = sf_organizer(stim_idx_temp,savepath,'variable_name',...
        'Kernels','filename',Kernel_name,'subfoldername',"Kernel",'variable',kernels_temp,...
        'overwrite',true,'update_data',false)
    
    
    %Kernels.Kernels = kernels_temp;
    %Save data to existing mat file
    
end
waitmessage.Destroy;


%Add cell indices to Kernel_location
kernel_idx = size(Kernel_location_temp,1);

if isempty(Kernel_location)
    Kernel_location(1:kernel_idx,1) = Kernel_location_temp;
    Kernel_location(1:kernel_idx,2) = cell_indices_temp; 
else
    
    Kernel_location(end:end+kernel_idx-1,1) = Kernel_location_temp;
    Kernel_location(end:end+kernel_idx-1,2) = cell_indices_temp; 
end


[~] = sf_organizer(stim_idx, savepath, 'subfoldername', "Kernel", 'collect_files',true);
[~] = sf_organizer(stim_idx,savepath,'variable_name','Kernel_location','variable',Kernel_location);
end
% S = matfile(savepath, 'Writable', true);
% S.Kernel_location = Kernel_location';


    

out = add_info;


%%
%out = sf_organizer(3,savepath,'variable_name','Test','variable',10,'filename','Test1')
end


end

        






     