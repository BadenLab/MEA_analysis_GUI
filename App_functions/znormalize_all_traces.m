function Kernel_normalized = znormalize_all_traces (in)

%This function znormalizes all Kernel traces based on the mean of the last
%third of each trace (trace after the spike happend.
Kernel_size = size(in);

%Check if its an array
if nnz(Kernel_size) == 0
    error("Array is empty")
    
end

%Check if array is 3d

if length(Kernel_size) < 3
    error("Requires 3d input, check array dimensions")
end

ac_length = Kernel_size(1);
ac_future_idx = ceil(ac_length*(2/3));

Kernel_normalized = NaN(Kernel_size(1),Kernel_size(2),Kernel_size(3));
parfor ii = 1:Kernel_size(2)
   
    Kernel = squeeze(in(:,ii,:));
    Kernel_future = Kernel(ac_future_idx:end,:);
    Kernel_future_mean = mean(Kernel_future,'all');
    ac_std = std(Kernel_future,[],'all');
    Kernel_normalized(:,ii,:) = znormalise(Kernel,Kernel_future_mean,ac_std); 
        
    
end






end