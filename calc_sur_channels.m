function channel_matrix = calc_sur_channels (channel, neighbours, boxes)


%create a matrix with all channel indices
channel_x = 1+2*neighbours;
channel_matrix = NaN(channel_x);
a = 1;
for i = -neighbours:neighbours
 channel_matrix(:,a) = ((channel+boxes*i-neighbours):1:(channel+boxes*i+neighbours))';
a = a+1;
end

%change negative channel indices to NaNs
channel_matrix(channel_matrix <1) = NaN;


end