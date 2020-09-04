%Plot group characteristics 
nr_colours = 4;
nr_bins = 20;
colour_string = {'m','b','g','r'};
colour_value = [1 0 1; 0 0 1; 0 1 0; 1 0 0];
colour_string_rep = repmat(colour_value,[height(plot_table_sorted) 1]);
colour_name = {'UV','Blue','Green','Red'};
freq_lim = 20;
count_lim = [0 25];
amplitude_lim = 60;
edges_freq = (0:freq_lim/nr_bins:freq_lim);
edges_peaks = (0:amplitude_lim/nr_bins:amplitude_lim); 
%Plot amplitudes
close all
a = 0;
for ii =  1:10
   Indices = plot_table_sorted.Indices(ii);
   Indices = Indices{:}';
   idx = zeros(length(Indices),1);
   for kk = 1:length(Indices)
        idx(kk) = find([cell_kernel_overview.cell_idx] == Indices(kk));
        cell_data(kk).details = cell_kernel_overview(idx(kk)).detailed_info;
        %Get data for different characteristics 
        for cc = 1:nr_colours
        
            On_peaks(kk,cc+nr_colours*a) = cell_data(kk).details.On_peaks(cc);
            Off_peaks(kk,cc+nr_colours*a) = cell_data(kk).details.Off_peaks(cc);
            On_frequency(kk,cc+nr_colours*a) = cell_data(kk).details.On_frequency(cc);
            Off_frequency(kk,cc+nr_colours*a) = cell_data(kk).details.Off_frequency(cc);
            violin_name{1,cc+nr_colours*a} = ['Group ',num2str(ii),' ',colour_name{cc}];
                       
        
        end
        
        
        
   end
   a = a+1;
   
end

On_peaks(On_peaks == 0) = NaN;
Off_peaks(Off_peaks == 0) = NaN;
On_frequency(On_frequency == 0) = NaN;
Off_frequency(Off_frequency == 0) = NaN;

a = figure
vp_on = violinplot(On_peaks,violin_name);

for ii = 1:length(vp_on)
    vp_on(ii).ViolinColor = colour_string_rep(ii,:);
end

hold on

vp_off = violinplot(-Off_peaks,violin_name);

for ii = 1:length(vp_off)
    vp_off(ii).ViolinColor = colour_string_rep(ii,:);
end
hline(0,'k')
title('On and Off peak amplitudes')
set(0, 'DefaultFigureRenderer', 'painters');
saveas(a,'On_Off_peaks','svg');



b = figure
vp_on_f = violinplot(On_frequency,violin_name);

for ii = 1:length(vp_on_f)
    vp_on_f(ii).ViolinColor = colour_string_rep(ii,:);
end

hold on

vp_off_f = violinplot(-Off_frequency,violin_name);

for ii = 1:length(vp_off_f)
    vp_off_f(ii).ViolinColor = colour_string_rep(ii,:);
end
hline(0,'k')
title('On and Off response frequencies')
ylabel ('Frequency in Hz')
set(0, 'DefaultFigureRenderer', 'painters');
saveas(b,'On_Off_freq','svg');




c = figure
vp_diff = violinplot(diff_peaks,violin_name);

for ii = 1:length(vp_diff)
    vp_diff(ii).ViolinColor = colour_string_rep(ii,:);
end































   
   %% Plot the peaks
   a = figure;
   hold on
   for cc = 1:nr_colours
       ax(cc) = subplot(nr_colours,1,cc);
       histogram(On_peaks(:,cc),edges_peaks,'FaceColor',colour_string{cc});
      mean1 = nanmean(On_peaks(:,cc));
       std1 = nanstd(On_peaks(:,cc));
       try
       xline(mean1,'k')
       xline(mean1+std1,'k')
       xline(mean1-std1,'k')
       end
       if cc ~= nr_colours
           set(gca,'xtick',[])
           set(gca,'xticklabel',[])
           set(gca,'ytick',[])
           set(gca,'yticklabel',[])
           
       end
       set(gca,'box','off')
       ylim(count_lim)
       hold off
   end
   linkaxes(ax,'xy')
   sgtitle(['On peaks group: ', num2str(ii)])
   xlabel('Peak amplitude, a.u.')
   ylabel('Counts')
   saveas(a,['On peaks group', num2str(ii)],'svg')
   
   b = figure;
   hold on
   for cc = 1:nr_colours
       ax(cc) = subplot(nr_colours,1,cc);
       histogram(Off_peaks(:,cc),edges_peaks,'FaceColor',colour_string{cc});
       mean1 = nanmean(Off_peaks(:,cc));
       std1 = nanstd(Off_peaks(:,cc));
       try
       xline(mean1,'k')
       xline(mean1+std1,'k')
       xline(mean1-std1,'k')
       end
         if cc ~= nr_colours
           set(gca,'xtick',[])
           set(gca,'xticklabel',[])
           set(gca,'ytick',[])
           set(gca,'yticklabel',[])
         end
       set(gca,'box','off')
        ylim(count_lim)
       hold off

   end
   linkaxes(ax,'xy')
   sgtitle(['Off peaks group: ', num2str(ii)])
   xlabel('Response frequency, a.u.')
   ylabel('Counts')
   saveas(b,['Off peaks group ', num2str(ii)],'svg')
   %% Plot the frequencies
   c = figure
   hold on
   for cc = 1:nr_colours
       ax(cc) = subplot(nr_colours,1,cc);
       histogram(On_frequency(:,cc),edges_freq,'FaceColor',colour_string{cc});
      mean1 = nanmean(On_frequency(:,cc));
       std1 = nanstd(On_frequency(:,cc));
       try
       xline(mean1,'k')
       xline(mean1+std1,'k')
       xline(mean1-std1,'k')
       end
        if cc ~= nr_colours
           set(gca,'xtick',[])
           set(gca,'xticklabel',[])
           set(gca,'ytick',[])
           set(gca,'yticklabel',[])
        end
       set(gca,'box','off')
        ylim(count_lim)
       
       
       hold off
   end
   linkaxes(ax,'xy')
   sgtitle(['On frequencies group: ', num2str(ii)])
   xlabel('Peak amplitude, a.u.')
   saveas(c,['On frequencies group ', num2str(ii)],'svg')
   
   
   d = figure
   hold on
   for cc = 1:nr_colours
       ax(cc) = subplot(nr_colours,1,cc);
       histogram(Off_frequency(:,cc),edges_freq,'FaceColor',colour_string{cc});
       mean1 = nanmean(Off_frequency(:,cc));
       std1 = nanstd(Off_frequency(:,cc));
       try
       xline(mean1,'k')
       xline(mean1+std1,'k')
       xline(mean1-std1,'k')
       end
       if cc ~= nr_colours
           set(gca,'xtick',[])
           set(gca,'xticklabel',[])
           set(gca,'ytick',[])
           set(gca,'yticklabel',[])
       end
       set(gca,'box','off')
        ylim(count_lim)
       hold off

   end
   linkaxes(ax,'xy')
   sgtitle(['Off frequencies group: ', num2str(ii)])
   xlabel('Response frequency, a.u.')
   ylabel('Counts')
    saveas(d,['Off frequencies group ', num2str(ii)],'svg')
   
   
  
       
   
   
   
    
    
