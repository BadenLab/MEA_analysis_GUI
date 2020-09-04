function out = most_active(Kernels, Channels_idx)
%This function finds the most active Channel from a 3D matrix

%Load the channels in question
nr_channels = numel(Channels_idx);
kernel_size = size(Kernels);
Channels = zeros(nr_channels,kernel_size(1),kernel_size(3));
max_Channels = zeros(nr_channels,kernel_size(3));

for ii = 1:nr_channels
    Channels(ii,:,:) = Kernels(:,Channels_idx(ii),:);
    Channels(ii,:,:) = normalize(Channels(ii,:,:),'center');
    for kk = 1:kernel_size(3)
        max_Channels(ii,kk) = squeeze(max(abs(Channels(ii,:,kk)),[],2));
    end
end

%find the maximal active channel
[~,active_channel] = max(max_Channels,[],'all','linear');
active_channel(2) = ceil(active_channel/kernel_size(3));
active_channel_idx = Channels_idx(active_channel(2));

out = active_channel_idx;





end