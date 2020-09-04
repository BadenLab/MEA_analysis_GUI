function out = Get_spatial_kernel_app(savepath, add_info)

%Extract Kernel information
Kernel_info = add_info.Kernel_info;
Kernel_nr = Kernel_info.Cell;
stim_idx = add_info.stim_idx;
%Check if Kernels have been calculated
try
    S = load(findfile_app(stim_idx,savepath,'Kernel_location'));
    kernel_location = S.Kernel_location;
%     S = load(savepath,'-mat','Stimulus_info');
%     kernel_location = S.Stimulus_info(Stimulus_nr).Kernel_location;
catch ME
     message = sprintf('Kernels have not been calculated, run calculate Kernel script first \n%s', ME.message);
     uiwait(warndlg(message));
     return
end
%Load the direction for the Kernel information

Index = find(str2double(kernel_location(:,2)) == Kernel_nr);

folder_name = string(kernel_location{Index,1});

S = load(folder_name);
kernels = S.Kernels;
reshape_factor = sqrt(size(kernels,2));
bins = size(kernels,1);
%Reshape Kernels to the checkerboard size
uv_kernel = reshape(squeeze(kernels(:,:,1)'),[reshape_factor,reshape_factor,bins]);
blue_kernel = reshape(squeeze(kernels(:,:,2)'),[reshape_factor,reshape_factor,bins]);
green_kernel = reshape(squeeze(kernels(:,:,3)'),[reshape_factor,reshape_factor,bins]);
red_kernel = reshape(squeeze(kernels(:,:,4)'),[reshape_factor,reshape_factor,bins]);

%Make greyscale matrix out of it
out.uv_kernel_grey = mat2gray(uv_kernel);
out.blue_kernel_grey = mat2gray(blue_kernel);
out.green_kernel_grey = mat2gray(green_kernel);
out.red_kernel_grey = mat2gray(red_kernel);
out.bins = bins;

end