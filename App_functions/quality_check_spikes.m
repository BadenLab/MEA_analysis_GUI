function out = quality_check_spikes (spiketimestamps,threshold)

    spiketimestamps_diff = diff(spiketimestamps,1);
    
      
    median_spikes = nanmedian(spiketimestamps_diff,1);
    mean_spikes = nanmean(spiketimestamps_diff,1);
    
    to_plot(1,:) = median_spikes;
    to_plot(2,:) = mean_spikes;
    
%     figure
%     bar(to_plot,'FaceColor','b');
    
    
    spike_percentile = prctile(spiketimestamps_diff,[25 75],1);
    
    
    out = find(spike_percentile(2,:) < threshold);
    
    disp([num2str(length(out)),' Cells have past the quality check']);
    
    
end