function pixel = channel_id_to_pixel (channel_id, number_channels)

%This scripts calculates the pixel number considering a rectangle
%arrangement of pixels from a channel number considering a linear
%arrangements of pixels
pixel_per_side = sqrt(number_channels);

pixel(1,1) = ceil(channel_id/pixel_per_side);

a = (pixel(1,1)-1) * pixel_per_side;

pixel(1,2) = channel_id - a;




end