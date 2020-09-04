
function out = bin_responses_app (savepath, add_info)
%This function bins the stimuli selected by the user. It bins the stimulus
%from one trigger signal to the next trigger signal with in the desired
%binsize. This means that each trigger signal is aligned to also be the
%begining of a bin. 



stim_idx = add_info.stim_idx;


%% Memory managment
M = matfile(savepath,'Writable',false);
spikesize = numel(M.spiketimestamps);
stx_before = size(M.spiketimestamps,2);


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


%% Essential information

%Get stim name list
S = load(savepath,'stim_list','cell_indices','Ch');
stim_list = S.stim_list;
cell_indices = S.cell_indices;
Ch = S.Ch;
SamplingFrequency = Ch.SamplingFrequency;


stim_begin = add_info.stim_begin;
stim_end = add_info.stim_end;

%Old stuff
% %Get the right spiketimestamps
% spiketimestamps = load_spiketimestamps_app (savepath, add_info);

% stx = size(spiketimestamps,2);
% sty = size(spiketimestamps,1);
% stz = size(spiketimestamps,3);

for ii = 1:length(stim_begin) %Loop over stimuli
    %Ask for the binsize (can be different for each stimulus)
    
    prompt = {'Enter the binsize you want to use'};
    dlgtitle = 'Binsize';
    dims = [1 35];
    definput = {'0.05'};
    binsize = str2double(inputdlg(prompt,dlgtitle,dims,definput));
    
        
    %Name of the stimulus
    stim_name = stim_list{stim_idx(ii)};
    stim_duration = stim_end(ii) - stim_begin(ii);
    stim_begin_temp = stim_begin(ii);
    stim_end_temp = stim_end(ii);
   
    
    %trigger channel
    stim_begin_frame = stim_begin_temp*SamplingFrequency;
    stim_end_frame = stim_end_temp*SamplingFrequency;
    trigger_ch = Ch.Ch01_02(stim_begin_frame:stim_end_frame);
    
    %identify trigger signals
    noise_trigger_norm = trigger_ch(1,:) > 2500;
    %find the peaks in the differences (which correspond to the beginning of a
    %trigger event
    [~,locs] = findpeaks(gather(double(noise_trigger_norm)),'MinPeakProminence',1,'MinPeakDistance',178);
    locs_s = zeros(1,length(locs)+1);
    locs_s(1) = 0;
    locs_s(2:end) = locs/SamplingFrequency;
    
          
    locs_diff = diff(locs_s);
    locs_diff_norm = ceil05(locs_diff,binsize);
    locs_s(end+1) = locs_s(end)+mean(locs_diff_norm);
    locs_s_real = locs_s+stim_begin_temp;
    nr_ep = length(locs_diff);
    
    
    for bb = 1:batch_nr %Loop over batches of spikes
        
    spiketimestamps = load_spiketimestamps_app (savepath, add_info,...
    'cell_subset',[batch_begins(bb),batch_ends(bb)]);
    cell_indices_temp = cell_indices(batch_begins(bb):batch_ends(bb));
    stx = size(spiketimestamps,2);
    % sty = size(spiketimestamps,1);
    % stz = size(spiketimestamps,3);
    
    spiketimestamps(isnan(spiketimestamps)) = 0;
   
    
    %Here we calculate the last actual value in the matrix for every cell
    %so that we can skip the NaN values (saves time)
    indexL = zeros(1,stx);
    for ll = 1:stx
        try
        indexL(1,ll) = find(spiketimestamps(:,ll),1,'last');
        catch
            continue
        end
    end
    
       
    for bi = 1:nr_ep
        
        Bined_epochs(bi).nr_bins = locs_diff_norm(bi)/binsize;
        Bined_epochs(bi).bin_start = (locs_s_real(bi):binsize:locs_s_real(bi)+...
            locs_diff_norm(bi)-binsize);
        Bined_epochs(bi).bin_end = Bined_epochs(bi).bin_start+binsize;
        Bined_epochs(bi).hist_bins = Bined_epochs(bi).bin_start;
        Bined_epochs(bi).hist_bins(end+1) = Bined_epochs(bi).bin_end(end);
         %This is stupid but in the histocount function the first edge has to be the
        %left edge of the first bin but the last edge has to be the right edge of
        %the last bin. So we have to change the last entry in the bin_start_1
        %matrix accordingly
                         
    end
    Bined_epochs(1).binsize = binsize;
    max_bins = int64(max([Bined_epochs.nr_bins]));
    
    %% Preallocate structures
    for kk = 1:stx
        
        Bined_spikes(kk+batch_ends*(bb-1)).histcounts = zeros;
        Bined_spikes(kk+batch_ends*(bb-1)).cell_idx = cell_indices_temp(kk);
                      
    end
    %%

    for kk = 1:stx
        if indexL(1,kk) == 0 %This checks if there are any spikes at all for that cell
            continue
        else
            spikes_temp = spiketimestamps(1:indexL(1,kk),kk); %Get spikes from one cell
            for bi = 1:nr_ep
                nr_bins = int64(Bined_epochs(bi).nr_bins); %Results are wrong if double precision is used
                if spikes_temp(indexL(1,kk)) < Bined_epochs(bi).hist_bins(1) %Check if last spike is within bin time area
                    continue
                else
                    [Bined_spikes(kk+batch_ends*(bb-1)).histcounts(1:nr_bins,bi), ~]...
                        = histcounts(spikes_temp,Bined_epochs(bi).hist_bins);
                                    

                end
            end
        end
        
        
    
    end
    
    
    
    %% Sort spikes for the raster plot
    
    %preallocate, by calculating sum of spikes in bins
    for kk = 1:stx
        nr_spikes = max(sum(Bined_spikes(kk).histcounts,1));
        Raster_spikes(kk+batch_ends*(bb-1)).spikes = NaN(nr_spikes,nr_ep);
        Raster_spikes(kk+batch_ends*(bb-1)).cell_idx = cell_indices_temp(kk);
    end
    
    
    
    for kk = 1:stx
        
        spikes_temp = spiketimestamps(1:indexL(1,kk),kk); %Get spikes from one cell
        for bi = 1:nr_ep
            %Collect spikes and sort them accodring to the episodes
            spikes_ep_log = logical((spikes_temp > Bined_epochs(bi).hist_bins(1))...
                .*(spikes_temp < Bined_epochs(bi).hist_bins(end)));
            if nnz(spikes_ep_log) == 0
                continue

            else

            spikes_ep_temp = spikes_temp(spikes_ep_log)-Bined_epochs(bi).hist_bins(1);
            nr_spikes = length(spikes_ep_temp);

            Raster_spikes(kk+batch_ends*(bb-1)).spikes(1:nr_spikes,bi) = spikes_ep_temp;

            end
        end
    
    
    
    end

    
    Bined_spikes(1).bins_info = Bined_epochs;
   
   %%
   [~] = sf_organizer(stim_idx(ii),savepath,'variable_name','Bined_spikes',...
       'variable',Bined_spikes);
   [~] = sf_organizer(stim_idx(ii),savepath,'variable_name','Raster_spikes',...
       'variable',Raster_spikes);
   
    
end


out = 1;    
end