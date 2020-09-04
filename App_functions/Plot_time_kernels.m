function out = Plot_time_kernels (savepath, add_info)
nr_boxes = add_info.Kernel_info.neighbours;
Kernel_nr = add_info.Kernel_info.Cell; 
Kernel_nr
%Stores the information about which cells is shown
%load kernel info
S = load(findfile_app(add_info.stim_idx,savepath,'Kernel_info'));
Kernel_info = S.Kernel_info;
idx = find(str2double([Kernel_info.cell_idx]) == Kernel_nr);
Kernel_info = Kernel_info(idx);



%load Kernel data
try
    S = load(findfile_app(add_info.stim_idx,savepath,'Kernel_location'));
    kernel_location = S.Kernel_location;
catch ME
     message = sprintf('Kernels have not been calculated, run calculate Kernel script first \n%s', ME.message);
     uiwait(warndlg(message));
     return
end
folder_name = kernel_location(idx,1);
S = load(folder_name);
Kernels = S.Kernels;
boxes = sqrt(size(Kernels,2));
unique_channels = unique(Kernel_info.true_kernel_idx);
unique_channels(unique_channels == 0) = [];
nr_channels = numel(unique_channels);
nr_colours = size(Kernels,3);
if nr_channels > 10
    nr_channels = 10;
end
% colour array
Colour={'#7E2F8E', '#0072BD', '#77AC30', '#A2142F', 'k', 'y', 'c'};



%%
for kk = 1:nr_channels
    
    %Identify the channels
    all_channels = calc_sur_channels(unique_channels(kk),nr_boxes, boxes);
    x_y = size(all_channels,1);
    all_channels_re = reshape(all_channels',[1,x_y^2]);
    figure
    hold on
    for jj = 1:numel(all_channels_re)
        subplot(x_y,x_y,jj)
        for ii = 1:nr_colours
            try
                trace = gpuArray(Kernels(:,all_channels_re(jj),ii));
                plot(trace, 'color', Colour{ii}, 'LineWidth', 1.5)
                title(num2str(all_channels_re(jj)))
                hold on
            catch
                fail = fail+1;
                continue
                
            end
        
        end
        xline(200)
    end



end
out = 1;





end
