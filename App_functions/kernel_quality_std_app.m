function out = kernel_quality_std_app (savepath, add_info)


std_threshold = add_info.Kernel_info.std_threshold; %50
stim_idx = add_info.stim_idx;
S = load(findfile_app(stim_idx,savepath,'Kernel_location'));
Kernel_location = S.Kernel_location;
nr_Kernels = size(Kernel_location,1);
%load example trace to get the number of boxes
%K = load(Kernel_location(1,1));
%Kernel = K.Kernels;
%nr_pixel = sqrt(size(Kernel,2));
max_detect = 5;
%bins = size(Kernel,1);

true_kernel_pixel = zeros(1,max_detect);
h = waitbar(0,'Please wait...');
for kk = 1:nr_Kernels
    waitbar(kk/nr_Kernels,h,[num2str(kk),' of ', num2str(nr_Kernels)]);
    try
    S = load(Kernel_location(kk,1));
    Kernel = S.Kernels;
    for ff = 1:size(Kernel,3)
        Kernel_temp = Kernel(:,:,ff);
        Kernel_std = std(Kernel_temp);
        Kernel_std_all = std(Kernel_temp,0,'all');
        Active_Channels = find(Kernel_std>Kernel_std_all*std_threshold);
        nr_active_Channels = nnz(Active_Channels);
        
        if nr_active_Channels > max_detect
            nr_active_Channels = max_detect;
            AC_values = Kernel_std(Active_Channels);
            [~,I] = sort(AC_values,'descend');
            Active_Channels = Active_Channels(I(1:max_detect));
        end
        
        true_kernel_pixel(1:nr_active_Channels) = Active_Channels;
        
        try
            Kernel_info(kk).true_kernel_idx(:,ff) = true_kernel_pixel';
            Kernel_info(kk).nr_true_kernel(1,ff) = nnz(Kernel_info(kk).true_kernel_idx(:,ff));
    
        catch
        continue
        end
        true_kernel_pixel = zeros(1,max_detect);
    end
   Kernel_info(kk).true_kernel_log = logical(Kernel_info(kk).nr_true_kernel);
   Kernel_info(kk).cell_idx = Kernel_location(kk,2);
    catch
        continue
    end
end
close(h)

out = sf_organizer(stim_idx,savepath,'variable_name','Kernel_info','variable',Kernel_info);


% Sobj = matfile(savepath, 'Writable',true);
% Sobj.Kernel_info = Kernel_info;
out = 1;

end