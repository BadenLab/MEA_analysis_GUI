function out = kernel_quality_app (savepath, add_info)


width_threshold = add_info.Kernel_info.width_thr; %50
peak_threshold_temp = add_info.Kernel_info.peak_thr; %100
S = load(savepath,'-mat',"Kernel_location");
Kernel_location = S.Kernel_location;
nr_Kernels = size(Kernel_location,2);
%load example trace to get the number of boxes
K = load(Kernel_location(1,1));
Kernel = K.Kernels;
nr_pixel = sqrt(size(Kernel,2));
%bins = size(Kernel,1);



h = waitbar(0,'Please wait...');
for ff = 1:nr_Kernels

waitbar(ff/1000,h,[num2str(ff),' of ', num2str(nr_Kernels)]);
Kernel_location_temp = Kernel_location(1,ff);
K = load(Kernel_location_temp);
Kernels_temp = K.Kernels;
for C_t = 1:4   

%% Old version of finding the kernels
%implay(I)
% %check for the max in the kernel traces
% max_kernels = max(kernels_reshape,[],3);
% %reshape the max kernels to make further operations easier
% max_kernels = reshape(max_kernels,[1,1600]);
% kernel_threshold = std(max_kernels)*thres_fac;
% %calculate the mean of the kernels
% mean_kernel = mean(max_kernels);
% %detect kernels which are higher than the threshold 
% true_kernels = max_kernels > kernel_threshold+mean_kernel;
% %check which channels pass the threshold
% true_kernel_idx = find(true_kernels);
% 
% %The same for negative kernels
% min_kernels = min(kernels_reshape,[],3);
% min_kernels = reshape(min_kernels,[1,1600]);
% kernel_threshold1 = std(min_kernels)*thres_fac;
% mean_kernel_min = mean(min_kernels);
% true_kernels_min = min_kernels < mean_kernel_min - kernel_threshold1;
% true_kernel_idx_temp = find(true_kernels_min);
% true_kernel_idx = [true_kernel_idx true_kernel_idx_temp];



%% New way of kernels
true_kernel_pixel = false(4,nr_pixel*nr_pixel);


parfor i = 1:nr_pixel*nr_pixel
Channel = C_t;
width_thres = width_threshold;
peak_threshold = peak_threshold_temp;
Kernels = Kernels_temp;

kernel_trace = squeeze(Kernels(:,i,Channel));
kernel_trace = smooth(kernel_trace,50);
%find the peaks with the parameters set
kernel_trace_norm = normalize(kernel_trace,'center','mean');
gate = max(kernel_trace_norm);
pixel_test = 0;
if gate>peak_threshold
    pixel_test = findpeaks(kernel_trace_norm,'MinPeakWidth',width_thres,'MinPeakHeight',peak_threshold,'MaxPeakWidth',150);
end

if pixel_test > 0
    true_kernel_pixel(1,i) = true;
end
%The same for negativ peaks (signal fliped up down)
pixel_test = 0;
gate = max(-kernel_trace_norm);
if gate> peak_threshold
    pixel_test = findpeaks(-kernel_trace_norm,'MinPeakWidth',width_thres,'MinPeakHeight',peak_threshold,'MaxPeakWidth',150);
end
pixel_test = abs(pixel_test);
if pixel_test > 0
    true_kernel_pixel(1,i) = true;
    
    
    
end

end
%%Save data into mat file
idx = numel(find(true_kernel_pixel));
try
Kernel_info(ff).true_kernel_idx(1:idx,C_t) = find(true_kernel_pixel(C_t,:));
Kernel_info(ff).nr_true_kernel(1,C_t) = nnz(Kernel_info(ff).true_kernel_idx(:,C_t));

catch
    continue
end

end
%Logical index of cells which pass the test
Kernel_info(ff).true_kernel_log = logical(Kernel_info(ff).nr_true_kernel);

end
close(h) %Close waitbar


Sobj = matfile(savepath, 'Writable',true);
Sobj.Kernel_info = Kernel_info;
out = 1;
end