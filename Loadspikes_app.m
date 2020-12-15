function out = Loadspikes_app
%This function loads spiketimes from the hdf5 file which is the output of
%the spikesorting (herdingspikes) algorithm. The script checks which
%spiketimes are empty and which cells pass the natural spiking frequency
%threshold. The stimulus traces is loaded and the stimulus names are
%assigned based on user input (function "findstimtimes"). The information
%is than stored in an array which is written to a analysis main file in
%the respective folder or based on user inoput at a location of choice. 
%@MSeifert 2019

%% variables
hz_threshold.min_distance = 0.005;
hz_threshold.hz_threshold = 0.3;

%% Get the Hdf5 file
%Load the sorted spikes
[hdfname,pathname] = uigetfile('*.hdf5','Select hdf5 file');
%pathname = app.working_path;
% want_spike_trace = 0
HdfFile=strcat(pathname,hdfname);
centres_temp = double(h5read(HdfFile,'/centres'));   % get the cluster localisation
centres_temp = centres_temp';
% columns 1&2 refer to x and y coordinates respectively
cluster_id =double(h5read(HdfFile,'/cluster_id'));  % the cluter id's for every spike
times = double(h5read(HdfFile,'/times'));  % the time ( in data point)
%for each spike (note cluster_id and times has the same length, becomes important later)
Sampling = double(h5read(HdfFile,'/Sampling'));
Channels = double(h5read(HdfFile,'/ch'));
%Shapes = double(h5read(HdfFile,'/shapes'));
time_in_s = times / Sampling; %Converts the information about frames into time in s

units = double(tabulate(cluster_id));
nunits = numel(units(:,1));
cell_indices = (1:1:nunits);
maxspikes = max(units(:,2));





%% Plot overview over what the hdf5 file contains
hdf_fig = figure;
hold on
bar(cell_indices, units(:,2),'k');
title("Number of spikes per cluster")
xlabel("units/cells")
ylabel("Nr of spikes")
set(gca, 'YScale', 'log')

%% Ask if subset of clusters should be loaded

answer = questdlg('Do you wannt to load the full set or a subset of clusters?', ...
	'Loading choice', ...
	'Full set','Subset','Full set');
% Handle response
switch answer
    case 'Full set'
        disp([answer 'Loading full set'])
        choice = 1;
    case 'Subset'
        disp([answer 'Loading subset'])
        choice = 2;
end
 
if choice == 2
    prompt = {'Enter index of first cell you want to load',...
        'Enter the index of the last cell you want to load'};
    dlgtitle = 'Input';
    dims = [1 40];
    definput = {'1',num2str(nunits)};
    answer = inputdlg(prompt,dlgtitle,dims,definput);
    batch_begin = str2double(answer{1,1});
    batch_end = str2double(answer{2,1});
    batch_size = batch_end-batch_begin+1;
    
    %Sort out the spikes we dont want
    nunits = batch_size;
    units = units(batch_begin:batch_end,:);
    cell_indices = cell_indices(1,batch_begin:batch_end);
     
     
    %Plot overview over subset selection
    hdf_fig = figure;
    hold on
    bar(cell_indices, units(:,2));
    title("Number of spikes per cluster")
    xlabel("cell index")
    ylabel("Nr of spikes")
   

    
end
  

%% Assigne spikes to array
% some things to do before we can strt. Here I generate a big matrix with
%all spikes - each column stands for a cluster (feel free to organise it your way)

%spiketimestamps = sparse(maxspikes,nunits);
spiketimestamps = spalloc(maxspikes,nunits,maxspikes);
spiketimestamps_empty = zeros(nunits,1);
%spiketimestamps = tall(zeros(maxspikes,nunits));
a = 1;

for i = 1:nunits
    i
    spiketimestamps(1:units(i,2),a)=(times(cluster_id==units(i,1)));
    empty_cluster = nnz(spiketimestamps(:,a));

    %This only saves the centers into the matrix which refer to a not
    %empty cell
    centres(a,:) = centres_temp(i,:);
    a = a+1;

end
%Transfer timestamps from frames to seconds
spiketimestamps = spiketimestamps / Sampling;
clear centres_temp
disp ("Data Loaded");


%% Old version

% for i = 1:nunits
%     i
%     spiketimestamps(1:units(i,2),a)=(times(cluster_id==units(i,1)));
%     empty_cluster = nnz(spiketimestamps(:,a));
%     if empty_cluster <100
%         %This deletes almost empty clusters and adds a true value into the
%         %matrix for empty cells
%         spiketimestamps(:,a) = [];
%         spiketimestamps_empty(i,1) = 1;
%         %delet the indices of rejected cells
%         cell_indices(:,a) = [];
% 
%     else
%         %This only saves the centers into the matrix which refer to a not
%         %empty cell
%         centres(a,:) = centres_temp(i,:);
%         a = a+1;
%     end
% 
% end





%% Spiketime validation
%This rejects all cells which spike at a too high rate (10% or more of the
%cells have a inter spike intervall of less than 10 milliseconds)
test = 1;
if test == 0
stx = length(spiketimestamps(1,:));
sty = length(spiketimestamps(:,1));
st_validated = NaN(1,stx);

if nunits < 200 %For few clusters no parfor loop is established

    for i=1:length(spiketimestamps(1,:))

    min_distance = hz_threshold.min_distance; %This is the maximal firing rate that is considered
    %min_distance = 0.0005;
    ISP_temp = spiketimestamps(:,i);
    sty = nnz(ISP_temp);
    ISP = diff(ISP_temp(1:sty));
    ISP_neg = logical(ISP<min_distance);
    ISP_neg_temp = nnz(ISP_neg);
%     nnz_temp = nnz(ISP_temp==0);
%     nr_spikes = sty-nnz_temp;
%     false_spikes = nr_spikes/ISP_neg_temp;
    false_spikes = ISP_neg_temp/sty;
    st_validated(1,i) = false_spikes;
    end
else
    
    parfor i=1:length(spiketimestamps(1,:))

        min_distance = hz_threshold.min_distance; %This is the maximal firing rate that is considered
        %min_distance = 0.0005;
        ISP_temp = spiketimestamps(:,i);
        sty = nnz(ISP_temp);
        ISP = diff(ISP_temp(1:sty));
        ISP_neg = logical(ISP<min_distance);
        ISP_neg_temp = nnz(ISP_neg);
    %     nnz_temp = nnz(ISP_temp==0);
    %     nr_spikes = sty-nnz_temp;
    %     false_spikes = nr_spikes/ISP_neg_temp;
        false_spikes = ISP_neg_temp/sty;
        st_validated(1,i) = false_spikes;
    end
end

% Plot results of the test
max_hz = 1/hz_threshold.min_distance;
freq_fig = figure;
hold on
bar(cell_indices,st_validated,'k')
title(["Number of cells that show response frequencies above ",max_hz," Hz"]);
ylabel("Nr of cells in %")
xlabel("Cell index")
hline(hz_threshold.hz_threshold,'-','Threshold')

dim = [.2 .5 .3 .3];
str = 'High number of cells above the threshold indicates wrong clustering bandwidth in spikesorting';
annotation('textbox',dim,'String',str,'FitBoxToText','on');


%Clusters with more than 10% of wrong spikes will be sorted out
wrong_spikes = logical(st_validated>hz_threshold.hz_threshold);
spiketimestamps(:,wrong_spikes) = [];
cell_indices(wrong_spikes) = [];
centres(wrong_spikes,:) = [];
spiketimestamps_empty(wrong_spikes,1) = 2;
stx = stx-nnz(wrong_spikes);
%find last nonzero element and reduce matrix size accordingly
nnz_row = NaN(1,stx);
for i = 1:stx
    [nnz_row(i),~] = find(spiketimestamps(:,i),1,'last');
end
st_end = nanmax(nnz_row);
spiketimestamps = spiketimestamps(1:st_end,:);
clear ISP
clear ISP_neg
end
%  %% autocorrelation
%  %Extract a subset of spikes in the middle of the recording
%  mean_time = mean(median(spiketimestamps,1));
%  spike_subset = spiketime_extract(spiketimestamps,mean_time,mean_time+60);
% repeats = ceil(size(spike_subset,2)/Batch_size);
%
% for i = 1:repeats
%     spike_corr = xcorr(Batch(i,Batch_size,spike_subset));
% end
%


%% Get the stimulus time
%In this section the trigger signal is read and information about the
%stimulus is saved in different matrices

%  %Load the stimulus file (Channel 1,2)
%  [stimulus_file,stimulus_path] = uigetfile('*.mat','Select stimulus file');
%  cd (stimulus_path);
%  load(stimulus_file);

%This function recognizes the different Epochs that were used
%during the recording, based on the Epoch codes in the trigger channel
%Epochs = Ar
%duino_Epochs(Ch01_02, 200, SamplingFrequency);

[stim_begin, stim_end, filename, Ch ] = findstimtimes;
stim_begin = stim_begin/Sampling;
stim_end = stim_end/Sampling;
nr_stimuli = length(stim_begin);
stim_nr = string((1:1:nr_stimuli));
prompt = stim_nr;
dlgtitle = 'Stim_name';
dims = [1 100];
stim_list = inputdlg(prompt,dlgtitle,dims);


%create logical arrays for most used stimuli
out.FFF = ~cellfun(@isempty,strfind(stim_list,'FFF'));
out.noise = ~cellfun(@isempty,strfind(stim_list,'noise'));
out.gratings = ~cellfun(@isempty,strfind(stim_list,'grating'));
out.chirp = ~cellfun(@isempty,strfind(stim_list,'chirp'));
FFF = out.FFF;
noise = out.noise;
gratings = out.gratings;
chirp = out.chirp

out.spiketimestamps = spiketimestamps;
out.spiketimestamps_empty = spiketimestamps_empty;
all_clusters = length(spiketimestamps(1,:));
all_spikes = numel(spiketimestamps);
out.all_clusters = all_clusters;
out.all_spikes = all_spikes;
out.centres = centres;
out.filenames = char(HdfFile);
out.stim_begin = stim_begin;
out.stim_end = stim_end;
out.stim_list = stim_list;
out.pathname = pathname;
out.Ch = Ch;
out.cell_indices = cell_indices;
out.hz_threshold.min_distance = hz_threshold.min_distance;
out.hz_threshold.hz_threshold = hz_threshold.hz_threshold;
%Save the data into a .mat file in the working directory
out.savename = strrep(HdfFile, '.hdf5', '_analysis.mat');
savename = out.savename;


%% Check if savefile already exists
choice1 = 0;
if exist(savename,'file') == 2
    answer = questdlg('Do you want to override the existing analysis file?', ...
	'Save choice', ...
	'Yes','No','No');
    % Handle response
    switch answer
    case 'Yes'
        disp([answer ', overwriting existing file'])
        choice1 = 1;
    case 'No'
        disp([answer 'Enter filename'])
        choice1 = 2;
    end
end


if choice1 == 2
     [file, path] = uiputfile('*.mat');
     out.savename = strcat(path,file);
end
%% Save
savename = out.savename;

save_var ={'FFF','noise','gratings','chirp','spiketimestamps',...
    'spiketimestamps_empty','all_clusters','all_spikes','centres',...
    'filename','stim_begin','stim_end','stim_list', 'pathname',...
    'savename', 'Ch', 'cell_indices', 'hz_threshold'...
    };
save([out.savename], save_var{:})
savefig(hdf_fig,[out.pathname,'hdf_figure.fig']);
%savefig(freq_fig,[out.pathname,'freq_figure.fig']);

end