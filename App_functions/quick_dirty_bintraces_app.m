function out = quick_dirty_bintraces_app (savepath, add_info)
stim_idx = add_info.stim_idx;
S = load(savepath,'pathname');
pathname = S.pathname;
stim_folder = strcat(pathname, "Stimulus_", num2str(stim_idx),"\");
Bin_name = ("Bins.mat");
cd(stim_folder);
S = load(Bin_name);
binsize = S.binsize;
Binned_sparse = S.Binned_sparse;

S = load(savepath,'Ch');
Ch = S.Ch;
stim_ch = Ch.Ch01_02;
SF = Ch.SamplingFrequency;
stim_ch = stim_ch(add_info.stim_begin*SF:add_info.stim_end*SF);
stim_ch = downsample(stim_ch,100);




prompt = {'Enter the cell you want to use'};
dlgtitle = 'Cell Number';
dims = [1 35];
definput = {'1'};
Cell_number = str2double(inputdlg(prompt,dlgtitle,dims,definput));

Cell = full(Binned_sparse(:,Cell_number));
x_values = (0:binsize:size(Cell,1)*binsize-binsize);
x_values_stim = (0:1/Ch.SamplingFrequency*100:size(Cell,1)*binsize);

length_diff = length(stim_ch)-length(x_values_stim);
stim_ch = stim_ch(1:end-length_diff);
length(stim_ch)
length(x_values_stim)
x_values = x_values/60;
x_values_stim = x_values_stim/60;

%x_values_stim = x_values_stim/60;

ax1 = subplot(2,1,1);
plot(x_values,Cell,'k');
ax2 = subplot(2,1,2);


plot(x_values_stim,stim_ch, 'k');
out = 1;
linkaxes([ax1 ax2], 'x');


end