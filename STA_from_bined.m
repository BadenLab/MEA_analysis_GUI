trace_average(trace_average < 0) = 0;
[peaks, spiketimes] = findpeaks(trace_average,'MinPeakHeight',0.15,'MinPeakProminence',0.15);
spiketimes_frames = spiketimes*binsize;
a = 1;
for ii = 1:length(peaks)
   upsampled = repelem(spiketimes_frames(ii),ceil(peaks(ii)));
   b = a+length(upsampled)-1; 
   spiketime_frames_up(a:b) = upsampled;   
   a = length(spiketime_frames_up)+1;
  
   
end

spiketimes_frames_int = floor(spiketimes_frames);
frames_start = spiketimes_frames_int-6;
frames_start(frames_start < 0) = 1;


hdf5_file = "D:\Data_MEA_Setup\QDSpy\Stimuli\stimuli_chicken_01_12_19\waves\Noise.h5";
a = 1;
for ii = 1:length(spiketimes_frames_int)
    hdf5_start = [1 frames_start(ii)];
    hdf5_count = [400 6];
    Colour_noise_temp(a:a+5,:,1) = (h5read(hdf5_file, '/UV_Noise', hdf5_start, hdf5_count))';
    Colour_noise_temp(a:a+5,:,2) = (h5read(hdf5_file, '/Blue_Noise', hdf5_start, hdf5_count))';
    Colour_noise_temp(a:a+5,:,3) = (h5read(hdf5_file, '/Green_Noise', hdf5_start, hdf5_count))';
    Colour_noise_temp(a:a+5,:,4) = (h5read(hdf5_file, '/Red_Noise', hdf5_start, hdf5_count))';
    a = a+6;
end

STA = squeeze(sum(Colour_noise_temp,1));
STA_boxes = reshape(STA,[20,20,4]);
UV = squeeze(STA_boxes(:,:,1));
UV_gray = mat2gray(UV);
STA_SD = squeeze(std(double(Colour_noise_temp),[],1));
SD_UV = reshape(STA_SD(:,1),[20,20]);
SD_UV_gray = mat2gray(SD_UV);