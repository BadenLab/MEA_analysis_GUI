function out = plot_all_cells (savepath,add_info)
%% Plot all cells function

%This function plots a raster plot for all cells, or if more than 1500
%cells are in the recordings a random subset of 1500 cells
delete(add_info.settings.whole_recording.panel.Children)
%% Load spikes
M = matfile(savepath,'Writable',false);
spikesize = numel(M.spiketimestamps);

%One double entry is 8 bytes, we dont want to use more than one fourth of the
%available memory
memory_available = memory;
memory_available = memory_available.MaxPossibleArrayBytes/4;
array_size_available = memory_available/8;

%If loading all spiketimestamps requires more than half the memory
%available, batch mode is activated.

    %Calculate batch size
batch_nr = ceil(spikesize/array_size_available);

Ch = M.Ch;
trig_channel = Ch.Ch01_02;
max_time = length(trig_channel)/Ch.SamplingFrequency;
trig_channel_ds = downsample(trig_channel,100);
xvalue = (0:1/(Ch.SamplingFrequency/100):max_time);

xvalue_l = length(xvalue);
trig_channel_ds_l = length(trig_channel_ds);



if ~(trig_channel_ds_l == xvalue_l)
    
    length_diff = trig_channel_ds_l - xvalue_l;
    
    if length_diff > 0
        xvalue(xvalue_l:trig_channel_ds_l) = 0;
    else
        xvalue(trig_channel_ds_l+1:xvalue_l) = [];
    end
        
    
    
end
    


if batch_nr == 1
    
    spiketimestamps = M.spiketimestamps;
    spiketimestamps = full(spiketimestamps);
    spiketimestamps = single(spiketimestamps(:,1:500));
    
    [N,edges] = histcounts(spiketimestamps,'BinWidth',1);
    N = N(2:end);
    edges = edges(2:end-1);
        
    ax(1) = subplot(3,1,1,'Parent',add_info.settings.whole_recording.panel);
    plot_raster(spiketimestamps,'ax',ax);
    ax(2) = subplot(3,1,2,'Parent',add_info.settings.whole_recording.panel);
    plot(ax(2),edges,N,'k');
    ax(3) = subplot(3,1,3,'Parent',add_info.settings.whole_recording.panel);
    plot(ax(3),xvalue,trig_channel_ds);
   linkaxes(ax,'x')
    out = ax(1,1).XLim;
    
end
    
    
    
    
    
    
    


end