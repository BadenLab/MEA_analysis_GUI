function out = gratings_qc (savepath, add_info)


%% Load Data
stim_idx = add_info.stim_idx;
try
    S = load(findfile_app(stim_idx,savepath,'dir'));
    dir = S.dir;
    
catch
    error('Gratings analysis file could not be found, analyse gratings first');
end
S = load(savepath,'cell_indices');

%% Global variables
cell_indices = S.cell_indices;
nr_cells = size(dir(1).spike_nr_repeats(1,:,1),2);
plot_data = false;
nr_directions = size(dir,2);
stimulus_frequency = 0.75;
stimulus_frequencies = (stimulus_frequency:stimulus_frequency:3*stimulus_frequency);
tolr = 0.25;

%% Preallocate output matrices
cell_spikes = squeeze(dir(1).spike_nr_repeats(:,1,:));
cell_spikes(cell_spikes == 0) = NaN;
nr_repeats = size(cell_spikes,2);
maxspikes_test = max(cell_spikes,[],'all'); 
Q_frequency = zeros(nr_cells,size(dir,2));
Q1 = zeros(nr_cells,size(dir,2));
Q_total = zeros(nr_cells,size(dir,2));

%% Loop
parfor cc = 1:nr_cells
    cc
for dd = 1:nr_directions

%% Load data for current direction
cell_spikes = squeeze(dir(dd).spike_nr_repeats(:,cc,:));
cell_spikes(cell_spikes == 0) = NaN;

%% Bin Data
%Transform data to a series of 0 and 1 at a given bin size
binsize = 0.001; %To create a matrix with 0 and 1 (spike or no spike)
maxspikes = nanmax(cell_spikes,[],'all');
if isnan(maxspikes)
    continue
end
nr_bins = ceil(maxspikes/binsize);

max_bin_time = binsize*nr_bins; %Last edge of the histcount

edges = (binsize:binsize:max_bin_time);

cell_histcount = zeros(nr_bins-1,size(cell_spikes,2));

%Bin all traces
for ii = 1:nr_repeats
    try
    cell_histcount(:,ii) = histcounts(cell_spikes(:,ii),edges);
    catch
        continue
    end
end



%% convolve with gaussian kernel
%each spike falls into a certain time-probability window determined by a
%gaussian kernel

w = gausswin(100);
spikes_smoothed = zeros(size(edges,2)-1,nr_repeats);
for ii = 1:nr_repeats
    
    spikes_smoothed(:,ii) = conv(cell_histcount(:,ii),w,'same');
   
end

%% Plot overview
if plot_data
    ax = [];
    ax(1) = subplot(nr_repeats+1,1,1);
    plot_raster(cell_spikes,1);

for ii = 1:size(cell_spikes,2)
   ax(ii+1) = subplot(nr_repeats+1,1,ii+1);
   plot_data
   plot(edges(1:end-1),spikes_smoothed(:,ii),'k');
   
end
linkaxes(ax,'x')
end


 %% Calculate average trace (old stuff)
% % point wise multiplication
% spikes_smoothed1 = spikes_smoothed+1;
% trace_average = [];
% trace_average(:,1) = spikes_smoothed1(:,1);
% for ii =1:size(spikes_smoothed1,2)-1
%     trace_average = trace_average(:).*spikes_smoothed1(:,ii+1);
% end
% trace_average = trace_average-1;
% %Normalize average trace
% %trace_average = trace_average/max(trace_average);
% if plot_data
%     figure
%     bx = [];
%     bx(1) = subplot(2,1,1);
%     plot_raster(cell_spikes,1);
% 
%     bx(2) = subplot(2,1,2);
%     plot(edges(1:end-1),trace_average,'k')
%     linkaxes(bx,'x')
% end

%% Compare to classical bin approach (old stuff)

binsize_test = 0.1;
%maxspikes_test = max(cell_spikes,[],'all'); 
nr_bins = ceil(maxspikes_test/binsize_test);

max_bin_time_test = binsize_test*nr_bins; %Last edge of the histcount

edges_test = (binsize_test:binsize_test:max_bin_time_test);

cell_histcount_test = zeros(nr_bins-1,size(spikes_smoothed,2));
% cell_gaussian = zeros(nr_bins-1,size(cell_spikes,2));
% figure
% hold on
for ii = 1:nr_repeats
    try
    cell_histcount_test(:,ii) = histcounts(cell_spikes(:,ii),edges_test);
    catch
        continue
    end
end

cell_test_smooth = smoothdata(cell_histcount_test,1,'gaussian',5);
cell_test_smooth_avr = mean(cell_test_smooth,2);
if plot_data
figure
tx = [];
t1x = [];
tx(1) = subplot(nr_repeats+1,1,1);
plot_raster(cell_spikes,1);
for ii = 1:nr_repeats
    tx(ii+1) = subplot(nr_repeats+1,1,ii+1);
    plot(edges_test(1:end-1),cell_test_smooth(:,ii),'k');
       
end
linkaxes(tx,'x')
linkaxes(tx(2:end),'y')
ylim([0,max(cell_test_smooth,[],'all')])

figure
t1x(1) = subplot(2,1,1);
plot_raster(cell_spikes,1);
t1x(2) = subplot(2,1,2);
plot(edges_test(1:end-1),cell_test_smooth_avr,'k');
linkaxes(t1x,'x')
end



%% Compare the power spectra of individial traces
if plot_data
    figure
end
P1 = [];
f1 = [];
P1_cut = [];
f1_cut = [];
sx = [];
for ii = 1:nr_repeats
Fs = 1/binsize;
try
[P1(:,ii),f1(:,ii)] = periodogram(spikes_smoothed(:,ii),[],[],Fs,'power');
%[P1(:,ii),f1(:,ii)] = periodogram(cell_test_smooth(:,ii),[],[],Fs,'power');
catch
    P1(:,ii) = zeros;
    f1(:,ii) = zeros;
    %max_P1(cc,ii,dd) = 0;
    continue
end
f_max = find((f1(:,ii) > 15),1,'first')-1;
f_min = find((f1(:,ii) > 0.5),1,'first')-1;
P1_cut(:,ii) = P1(f_min:f_max,ii);
f1_cut(:,ii) = f1(f_min:f_max,ii);
%[~,max_P1_temp] = max(P1_cut(:,ii));
%max_P1(cc,ii,dd) = f1_cut(max_P1_temp);
% [~,max_P1(2,ii)] = second_max(P1_cut(:,ii));
if plot_data
    
    sx(ii) = subplot(nr_repeats,1,ii);
    plot(f1_cut(:,ii),P1_cut(:,ii),'k')
    grid
    ylabel('P_1')
    title('Power Spectrum')
    linkaxes(sx,'xy')
    ylim([0,max(P1_cut,[],'all')])
end
% sx(2) = subplot(nr_repeats,1,i);
% plot(f2,P2,'r')
% grid
% ylabel('P_2')
% xlabel('Frequency (Hz)')


end

%% Compare mean power spectrum
%Compare the power spectra of these mean traces
try

%Take mean power spectrum    
P2_cut = mean(P1_cut,2);
P2_cut = P2_cut/max(P2_cut);

%find peak in that spectrum
[peaks,locs] = findpeaks(P2_cut);
[max_P2, max_idx] = max(peaks);
%max_P22 = second_max(peaks);

%peak_ratio = max_P2/max_P22;
Q_frequency(cc,dd) = f1_cut(locs(max_idx));

catch
    Q_frequency(cc,dd) = NaN;
    %peak_ratio = 0;
end

if plot_data
   figure
   plot(f1_cut,P2_cut);
    
end


%% Compare to sine wave
cross_coeff = zeros(nr_repeats,1);
for ii = 1:nr_repeats
stim_frequency = Q_frequency(cc,dd);

%Depending on the frequency the expected overlap is higher or lower.

max_peak_width = round(1/stim_frequency/binsize_test);
%Create two sine waves, one with the original frequency, and one with
%double the original frequency


fs_sine = 1/binsize_test; % Sampling frequency (samples per second)
dt = 1/fs_sine; % seconds per sample.
StopTime = max_bin_time_test; % seconds.
t = (0:dt:StopTime-dt)'; % seconds.
F = stim_frequency; % Sine wave frequency (hertz)
sine_org = sin(2*pi*F*t);
sine_org = sine_org(2:end);
sine_org = normalize(sine_org,'range');

% F1 = 2*stim_frequency;
% sine_double = sin(2*pi*F1*t);
% sine_double = sine_double(2:end);

%For debugging
% figure
% plot(sine_org)
% hold on
% plot(sine_double)

%Cross correlate the sine wave with the responses

% cross_coeff_double = zeros(nr_repeats,1);
try

   cross_coeff_trace = xcorr(sine_org,cell_test_smooth(:,ii),'coeff');
   [~,~,~,prominence] = findpeaks(cross_coeff_trace,'MinPeakProminence',0.1,'MaxPeakWidth',...
       max_peak_width,'Annotate','extents');
   
   %cross_coeff_double(ii) = max(xcorr(sine_double,cell_test_smooth(:,ii),'coeff'));
      
    %Compare to sine autocorrelation
    sine_coeff_trace = xcorr(sine_org,'coeff');
   [~,~,~,prominence_sine] = findpeaks(sine_coeff_trace,'MinPeakProminence',0.1,'MaxPeakWidth',...
   max_peak_width,'Annotate','extents');
   

    if isempty(prominence)
        cross_coeff(ii) = 0;
    else
        cross_coeff(ii) = max(prominence)/max(prominence_sine);
    end
    
%     cross_coeff_double(isnan(cross_coeff_double)) = 0;
   
    
catch
    continue
end
end
cross_coeff(isnan(cross_coeff)) = 0;
Q_total(cc,dd) = max(nanmean(cross_coeff));






 %% Compare power spectrum of mean trace (old_version)
% Compare the power spectra of these mean traces
% try
% Fs = 1/binsize;
% Fs_test = 1/binsize_test;
% %[P1,f1] = periodogram(trace_average(:,1),[],[],Fs,'power');
% [P2,f2] = periodogram(cell_test_smooth_avr(:,1),[],[],Fs_test,'power');
% % [P2,f2] = periodogram(trace_average(:,1),[],[],Fs_test,'power');
% P2_cut(:,1) = P2(f_min:f_max,1);
% f2_cut(:,1) = f2(f_min:f_max,1);
% [peaks,locs] = findpeaks(P2_cut);
% [max_P2, max_idx] = max(peaks);
% max_P22 = second_max(peaks);
% 
% peak_ratio = max_P2/max_P22;
% Q_frequency(cc,dd) = f2_cut(locs(max_idx));
% 
% catch
%     Q_frequency(cc,dd) = NaN;
%     peak_ratio = 0;
% end
% 
% 
% 
% 
% % sx(1) = subplot(2,1,1);
% % plot(f1,P1,'k')
% % grid
% % ylabel('P_1')
% % title('Power Spectrum')
% % 
% if plot_data
%     figure
%     plot(f2_cut,P2_cut,'r')
%     grid
%     ylabel('P_2')
%     xlabel('Frequency (Hz)')
% end
% % 
% % linkaxes(sx,'x')

% %% Calculate quality matrix
% 
% unique_peaks = unique(max_P1(cc,:,dd));
% count_similar = zeros(1,numel(unique_peaks));
% %Count how many peaks are observed at the same frequency
% for ii = 1:numel(unique_peaks)
%     if unique_peaks(ii) < 2*tolerance_idx
%         count_similar(ii) = 0;
%         continue
%     else
%     count_similar(ii) = nnz((max_P1(cc,:,dd) > unique_peaks(ii)-tolerance_idx)...
%         .*(max_P1(cc,:,dd) < unique_peaks(ii)+tolerance_idx));
%     end
% end
% Q1(cc,dd) = max(count_similar./8);
% Q_total(cc,dd) = Q1(cc,dd) * peak_ratio;


end
end


%% Calculate the final QC

[QC, QCI] = max(Q_total,[],2);
QC_frequency = zeros(nr_cells,1);
for ii = 1:nr_cells
    QC_frequency(ii,1) = Q_frequency(ii,QCI(ii));
end

QC_pass = logical(QC > 0.5);

for kk = 1:length(stimulus_frequencies)
    freq_pass(:,kk) = (QC_frequency > stimulus_frequencies(kk)-tolr).*...
        (QC_frequency < stimulus_frequencies(kk)+tolr);
end

freq_pass = sum(freq_pass,2);
QC_pass = logical(QC_pass.*freq_pass);


%% Create output file
for ii = 1:nr_cells
    Grating_QC(ii).QC = QC(ii,1);
    Grating_QC(ii).Q_total = Q_total(ii,:);
    Grating_QC(ii).QC_frequency = QC_frequency(ii,1);
    Grating_QC(ii).cell_idx = cell_indices(ii);
    Grating_QC(ii).QC_pass = QC_pass(ii); 
    Grating_QC(ii).freq_pass = freq_pass(ii);
end


[~] = sf_organizer(stim_idx,savepath,'variable_name','Grating_QC','variable',Grating_QC);

%  out = sf_organizer(stim_idx,savepath,'variable_name','dir','variable',dir);
%% Plot the histogram of QI
grating_qc_figure = figure;
bar(cell_indices,QC,'k')
title('Quality index for all cells')
xlabel('Cell ID')
ylabel("Quality index total")
ylim([0 1])
hline(0.5,'-',"Threshold")
%Save as matlab figure
sf_organizer(add_info.stim_idx,savepath,'variable_name','Grating_qc_figure',...
    'variable',grating_qc_figure);

grating_qc_hs = figure;
histogram(QC,'FaceColor','k')
title('QC histogram for all cells');
xlim([0 1])
vline(0.5,'-',"Threshold")
ylabel("Counts")
xlabel("Quality index, all cells")
%Save as matlab figure
sf_organizer(add_info.stim_idx,savepath,'variable_name','Grating_qc_hs',...
    'variable',grating_qc_hs);


figure
histogram(QC_frequency)
title('Frequency histogram')
xlabel('Frequency')
xlim([0 5])


%Plot comparison histogram
edges = (0:0.02:1);

idx_QC = [Grating_QC.QC_pass] == 1;
idx_freq = [Grating_QC.freq_pass] == 1;

grating_qc_summary = figure;
hold on
histogram([Grating_QC.QC],edges,'FaceColor','k')
histogram([Grating_QC(idx_QC).QC],edges)
histogram([Grating_QC(idx_freq).QC],edges)
legend("All cells","Passed QC, and frequency test","Passed frequency test");
hold off
ylabel('Cell count')
xlabel('Quality index')
title("Quality index for all cells summarized")
sf_organizer(add_info.stim_idx,savepath,'variable_name','grating_qc_summary',...
    'variable',grating_qc_summary);


idx05 = [Grating_QC.QC]>0.5;
frequency_qc_summary = figure;
edges1 = (0:0.1:5);
hold on
histogram([Grating_QC.QC_frequency],edges1,'FaceColor','k')
histogram([Grating_QC(idx05).QC_frequency],edges1)
histogram([Grating_QC(idx_QC).QC_frequency],edges1)
title("Frequency test summary")
ylabel("Cell counts")
xlabel("Frequency in [Hz]");
hold off
sf_organizer(add_info.stim_idx,savepath,'variable_name','frequency_qc_summary',...
    'variable',frequency_qc_summary);

% figure
% histogram([Grating_QC.],edges,'FaceColor','k')



out = 1;



end
















%% Experimental stuff (old)

% %% Test for cross correlation
% figure
% [c,lags] = xcorr(spikes_smoothed,'coeff');
% stem(lags,c,'color','k')
% 
% lag_window = 1; %In seconds
% lag_window_bins = lag_window/binsize;
% 
% %find 0 time lag
% 
% lag_0 = find(lags == 0);
% lag_cut(1) = lag_0-0.5*lag_window_bins+1;
% lag_cut(2) = lag_0+0.5*lag_window_bins;
% 
% c_new = c(lag_cut(1):lag_cut(2),:);
% 
% 
% %Find the autocorrelation indices
% auto_ones = true(size(c,2),1);
% auto_idx = (1:nr_repeats:size(c,2));
% add_idx = (0:1:7);
% auto_idx = auto_idx + add_idx;
% auto_ones(auto_idx) = false;
% 
% 
% 
% c_new_cross = c_new(:,auto_ones);
% figure
% stem(lags(lag_cut(1):lag_cut(2)),c_new_cross,'color','k')
% c_new_max = max(c_new_cross,[],1);
% figure
% boxplot(c_new_max);
% mean(c_new_max)
% 
% 
% 
% %%
% binsize = 0.001;
% maxspikes = max(cell_spikes,[],'all');
% nr_bins = ceil(maxspikes/binsize);
% 
% max_bin_time = binsize*nr_bins;
% 
% edges = (binsize:binsize:max_bin_time);
% 
% w = gausswin(50);
% 
% cell_histcount = zeros(nr_bins-1,size(spikes_smoothed,2));
% % cell_gaussian = zeros(nr_bins-1,size(cell_spikes,2));
% % figure
% % hold on
% for ii = 1:size(cell_spikes,2)
%     
%     cell_histcount(:,ii) = histcounts(cell_spikes(:,ii),edges);
% %     cell_gaussian(:,ii) = conv(cell_histcount(:,ii),w,'same');
%     plot(edges(2:end),cell_histcount(:,ii)+10*ii,'k')
% 
% end
% 
% 
% 
% % point wise multiplication
% trace_average(:,1) = spikes_smoothed(:,1);
% for ii =1:size(spikes_smoothed,2)-1
%     trace_average = trace_average(:).*spikes_smoothed(:,ii+1);
%     
% end
% 
% 
% 
% 
% combinations = combnk(1:size(cell_spikes,2),2);
% schreiber_test_result = zeros(length(combinations),1);
% 
% for ii = 1:length(combinations)
%     c1 = combinations(ii,1);
%     c2 = combinations(ii,2);
%     test_inner = dot(spikes_smoothed(:,c1),spikes_smoothed(:,c2));
% 
%     test_norm = norm(spikes_smoothed(:,c1));
%     test_norm2 = norm(spikes_smoothed(:,c2));
% 
% 
%     schreiber_test_result(ii) = test_inner / (test_norm*test_norm2);
% end
% figure
% boxplot(schreiber_test_result);
% 
% 
% 
% Fs = 1/binsize;
% T = binsize;
% L = length(cell_gaussian);
% t = (0:L-1)*T; 
% 
% Y = fft(cell_gaussian,[],1);
% fft_counts = zeros(floor(L/2+1),size(cell_spikes,2));
% fft_counts_var = zeros(1,size(cell_spikes,2));
% figure
% for ii = 1:size(cell_spikes,2)
%     P2 = abs(Y(:,ii)/L);
%     P1 = P2(ceil(1:L/2+1));
%     P1(2:end-1) = 2*P1(2:end-1);
%     fft_counts(:,ii) = P1;
%     fft_counts_var(ii) = var(fft_counts(1:95,ii));
% 
% 
%     f = Fs*(0:(L/2))/L;
%     
%     ax(ii) = subplot(size(cell_spikes,2),1,ii);
%     plot(f,P1) 
% end
% linkaxes(ax,'x');
% xlim([0,20])
% 
% 
% 
% fft_mean = mean(fft_counts(1:95,:),2);
% 
% fft_var_mean = var(fft_mean);
% fft_mean_var = mean(fft_counts_var);
% 
% QI = fft_var_mean/fft_mean_var;
% 
% test_weight = test_mean.*test_std;
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% %Experimental stuff
% 
% 
% [P1,f1] = periodogram(trace_average(:,1),[],[],Fs,'power');
% [P2,f2] = periodogram(cell_gaussian(:,3),[],[],Fs,'power');
% 
% ax(1) = subplot(2,1,1)
% plot(f1,P1,'k')
% grid
% ylabel('P_1')
% title('Power Spectrum')
% 
% ax(2) = subplot(2,1,2)
% plot(f2,P2,'r')
% grid
% ylabel('P_2')
% xlabel('Frequency (Hz)')
% 
% linkaxes(ax,'x');
% 
% [pk1,lc1] = findpeaks(P1,'SortStr','descend');
% P1peakFreqs = f1(lc1)
% 
% hold on
% 
% [pk2,lc2] = findpeaks(P2,'SortStr','descend');
% P2peakFreqs = f2(lc2)
% 
% 
% [Cxy,f] = mscohere(cell_gaussian(:,2),cell_gaussian(:,3),[],[],[],Fs);
% 
% 
% thresh = 0.75;
% [pks,locs] = findpeaks(Cxy,'MinPeakHeight',thresh);
% MatchingFreqs = f(locs)
% 
% 
% figure
% plot(f,Cxy)
% ax = gca;
% grid
% xlabel('Frequency (Hz)')
% title('Coherence Estimate')
% ax.XTick = MatchingFreqs;
% ax.YTick = thresh;
% axis([0 200 0 1])
% 
% 
% %More stuff
% 
% 
% level = 7;
% wpt = wpdec(cell_gaussian(:,2),level,'coif2');
% 
% figure;
% [S,T,F] = wpspectrum(wpt,Fs,'plot');
%     
% spectrogram(cell_histcount(:,2),100,90,10,'yaxis');
% wt = cwt(cell_histcount(:,2),Fs);
% 
% figure
% 
% ax(1) = subplot(2,1,1);
% hold on
% plot(edges(2:end),cell_histcount(:,2)+10*ii,'k')
% ax(2) = subplot(2,1,2);
% 
% plot(edges(2:end),cell_gaussian(:,2)+10*ii,'k')
% linkaxes(ax,'x');
% 








