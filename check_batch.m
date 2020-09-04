function out = check_batch (savepath)

%This function splits spiketimestamps into chunks of batches which can fit
%into memory at a time. It returns the indeces of the chunks relative to
%the total spiketimestamps, number of total batches and the total amount of
%spikes in the spiketimestmaps (to create final arrays to store all
%information)

%Check how many clusters we have
M = matfile(savepath,'Writable',false);
spikesize = numel(M.spiketimestamps);
stx_before = size(M.spiketimestamps,2);

%One double entry is 8 bytes, we dont want to use more than one fourth of the
%available memory
memory_available = memory;
memory_available = memory_available.MaxPossibleArrayBytes/8;
array_size_available = memory_available/8;

%If loading all spiketimestamps requires more than half the memory
%available, batch mode is activated.
if spikesize>array_size_available
    %Calculate batch size
    batch_nr = ceil(spikesize/array_size_available);
    batch_stx = ceil(stx_before/batch_nr);
    batch_begins = (1:batch_stx:stx_before);
    batch_ends = batch_begins(2:end)-1;
    batch_ends(end+1) = stx_before;
    batch_stx = batch_ends-(batch_begins-1);
    
    out.batch_nr = batch_nr;
    out.batch_stx = batch_stx;
    out.batch_begins = batch_begins;
    out.batch_ends = batch_ends;
    out.stx_before = stx_before;
else
    out.batch_nr = 1;
    out.batch_stx = stx_before;
    out.batch_begins = 1;
    out.batch_ends = stx_before;
    out.stx_before = stx_before;
    
end



    

