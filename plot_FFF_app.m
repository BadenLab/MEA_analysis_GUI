function out = plot_FFF_app (savepath, add_info)
FFF_select = add_info.FFF_select;
FFF_select
stim_idx = add_info.stim_idx;
FFF_panel = add_info.FFF_panel;
%% General variables
S = matfile(findfile_app(stim_idx,savepath,'FFF_average.mat'),'Writable',false);

S = S.FFF_average;
idx = [S.cell_idx] == FFF_select;
S = S(idx);
FFF_data = S.traces;

raster_spikes = S.spikes;
clear S

M = matfile(findfile_app(stim_idx,savepath,'Bined_spikes'),'Writable',false);
M = M.Bined_spikes(1,1);
binsize = M.bins_info(1).binsize;


nr_colours = size(FFF_data,2);
nr_bins_per_repeat = size(FFF_data,1);

stim_duration = nr_colours*nr_bins_per_repeat*binsize;
x_values = (binsize:binsize:stim_duration);

y_values = reshape(FFF_data,[1,size(x_values,2)]);
y_values = y_values*(1/binsize);

out.x_values = x_values;
out.y_values = y_values;






%% Raster plot
reshape1 = nr_colours*size(raster_spikes,1);
reshape2 = size(raster_spikes,2);
raster_spikes_plot = NaN(reshape1,reshape2);

for rr = 1:reshape2
    raster_spikes_plot(:,rr) = reshape(raster_spikes(:,rr,:),[reshape1,1]);
    raster_spikes_plot(raster_spikes_plot==0) = NaN;
end

ax1 = subplot(3,1,1,'parent',FFF_panel);

plot_raster(raster_spikes_plot,1,ax1);
raster_lim = size(raster_spikes_plot,2)*0.25+1;
ylim(ax1,[0,raster_lim]);
%ylabel(ax1,'Repeats')


%% Trace plot


% ax1 = 
ax2 = subplot(3,1,2,'parent',FFF_panel);
plot(ax2,x_values,y_values,'k')
%ylabel(ax2,'Spikes in [Hz]');


%% Stimulus plot


%Create colour array for plotting
stim_fig = ones(1, size(x_values,2));
RGB_values = zeros(size(x_values,2),3);
RGB = [1 0 0; 0 1 0; 0 0 1; 1 0 1];

for ii = 1:nr_colours
    idx_start = (ii-1)*nr_bins_per_repeat+1;
    idx_black_start = idx_start + nr_bins_per_repeat/2;

    RGB_values(idx_start:idx_black_start-1,:) = repmat(RGB(ii,:),floor(nr_bins_per_repeat/2),1);
end

%Plot stimulus
ax3 = subplot(3,1,3,'parent',FFF_panel);
stimulus_plot = bar(ax3,x_values,stim_fig,2,'FaceColor','flat','EdgeColor','flat');
%ylabel(ax3,'Stimulus');
%xlabel(ax3,'Time in [s]');

for kk = 1:length(RGB_values)
    
    stimulus_plot.CData(kk,:) = RGB_values(kk,:);
end
linkprop([ax1,ax2,ax3],{'XLim'});
%linkaxes([ax1,ax2,ax3],'x');


end
