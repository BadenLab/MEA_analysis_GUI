function spiketimes = spiketime_extract (spiketimestamps, min, max)
%% This function returns spikes in a given time interval

spiketimestamps_size = size(spiketimestamps);
spiketimes = NaN(spiketimestamps_size(1),spiketimestamps_size(2));
parfor i = 1:spiketimestamps_size(2)
    %Test which spikes are in the interval for a given cluster
    spiketimestamps_temp = full(spiketimestamps(:,i));
    spikes_log = logical((spiketimestamps_temp < max) .* (spiketimestamps_temp > min))
    spikes = spiketimestamps_temp(spikes_log);
    spikes_size = size(spikes);
    
    if isempty(spikes)
        continue
    else
    spikes(spikes_size+1:size(spiketimestamps_temp,1)) = NaN;
    spiketimes(:,i) = spikes;
    end
end



    


end