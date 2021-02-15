function out = analyse_gratings_app (savepath, add_info)

nr_dir = 8; %Number of directions 
batch_info = check_batch(savepath);
directions = [0,180, 45,225, 90,270, 135,315];
%directions_or = directions(1:length(directions)/2);

[dir_sorted, dir_order] = sort(directions);
%Extract the important information from input
dir_radians = deg2rad(dir_sorted);
dir_radians_or = dir_radians(1:2:end);

diff_dir = diff(dir_radians(1:2));
diff_dir_or = diff(dir_radians_or(1:2));
p_value = 0.05;
cut_factor = 1/5;


ch = load(savepath,'-mat','Ch');
ch = ch.Ch;
Ch01_02 = ch.Ch01_02;

cell_indices = load(savepath,'cell_indices');

loadpath = load(savepath,'pathname'); 
pathname = loadpath.pathname;

SamplingFrequency = ch.SamplingFrequency;

%Load cell indices
S = load(savepath,'-mat','cell_indices');
cell_indices = S.cell_indices;

%Check how many stimulus repeas are in the file
stim_repeats = numel(add_info.stim_begin);

for nn = 1:stim_repeats
    %Get information for the current stimulus
    stim_begin = add_info.stim_begin(nn);
    stim_end = add_info.stim_end(nn);
    stim_idx = add_info.stim_idx(nn);
    
    %Get the right trigger trace
    start_trigger = stim_begin*SamplingFrequency;
    end_trigger = stim_end*SamplingFrequency;
    
    trigger_ch = Ch01_02(int64(start_trigger-100):int64(end_trigger+100));
    diff_trigger_ch = diff(trigger_ch(1,:));
    trigger_ch_norm = diff_trigger_ch(1,:) > 500;
    %Get trigger times
    [~,locs] = findpeaks(gather(double(trigger_ch_norm)),'MinPeakProminence',1,'MinPeakDistance',178);
    nr_locs = length(locs);
    
  
%     if locs(1) ~= 0
%         locs = locs - locs(1); %In case the trigger channel started before the first trigger signal
%        
%     end
    
    locs_s = locs/SamplingFrequency;
    
    nr_repeats = ceil(nr_locs/nr_dir);
    repeat_idx = (1:1:nr_dir);
    repeat_idx = repmat(repeat_idx,[1,nr_repeats]); %This gives each trigger an index 
    
    
    %Next we need to load the spikes and loop over nr of batches (if
    %necessary)
    
    %Create output structure
        for cc = 1:batch_info.stx_before
            Data_circular(cc).spikes_deg = []; %Preallocate structure
            Data_circular(cc).circ_mean = [];
            Data_circular(cc).circ_var = [];
            Data_circular(cc).circ_std = [];
            Data_circular(cc).circ_confmean = [];
            Data_circular(cc).circ_skewness = [];
            Data_circular(cc).circ_kurtosis = [];
            Data_circular(cc).circ_rtest = [];
            Data_circular(cc).circ_otest = [];
            Data_circular(cc).circ_raotest = [];
            Data_circular(cc).circ_vmpdf = [];
            Data_circular(cc).circ_rtest_sig = 0;
            Data_circular(cc).circ_otest_sig = 0;
            Data_circular(cc).circ_r = [];
            Data_circular(cc).cell_idx = cell_indices(cc);
            Data_circular(cc).spikes_deg_w = [];
            
            % For orientation selectivity
            
            Data_circular(cc).spikes_deg_or = []; %Preallocate structure
            Data_circular(cc).circ_mean_or = [];
            Data_circular(cc).circ_var_or = [];
            Data_circular(cc).circ_std_or = [];
            Data_circular(cc).circ_confmean_or = [];
            Data_circular(cc).circ_skewness_or = [];
            Data_circular(cc).circ_kurtosis_or = [];
            Data_circular(cc).circ_rtest_or = [];
            Data_circular(cc).circ_otest_or = [];
            Data_circular(cc).circ_raotest_or = [];
            Data_circular(cc).circ_vmpdf_or = [];
            Data_circular(cc).circ_rtest_sig_or = 0;
            Data_circular(cc).circ_otest_sig_or = 0;
            Data_circular(cc).circ_r_or = [];
            Data_circular(cc).cell_idx_or = cell_indices(cc);
            Data_circular(cc).spikes_deg_w_or = [];
            
    
        end
        Data_circular(1).dir_radians = dir_radians;
    %%
    for bb = 1:batch_info.batch_nr
        
        
        %%

        new_info.stim_idx = stim_idx;
        new_info.stim_begin = stim_begin;
        new_info.stim_end = stim_end;
        spiketimestamps = load_spiketimestamps_app(savepath,new_info,...
            [batch_info.batch_begins(bb),batch_info.batch_ends(bb)]);
        %normalize spiketimestamps to stimulus begin
        spiketimestamps = spiketimestamps - stim_begin;
        
        %Preallocate
        for ss = 1:nr_dir %Preallocate strcture to save information for each direction
            dir(ss).spike_nr_repeats = [];
        end
        nr_repeats_real = zeros(1,nr_dir); %Saves how many repeats there are 
        %for real (stimulus may have end before all directions were shown
        %%
        
        for ii = 1:nr_dir
            %check which direction 
            trigger_idx = repeat_idx(ii);
            %the end of this direction is the beginning of the next
            end_idx = repeat_idx(ii+1);
            %find location of similar events throughout the stimulus
            event_idx = repeat_idx == trigger_idx;
            event_end = repeat_idx == end_idx;
            %Get the actual time of the begin and end for each trigger
            %signal
            trigger_start = locs(event_idx(1:length(locs)))/SamplingFrequency;
            trigger_end = locs(event_end(1:length(locs)))/SamplingFrequency;
            
            %I guess the last trigger end based on the mean intervalt, but I think thats ok here...
            if ii == nr_dir
                trigger_end = trigger_end(2:end);
             end
            
            %Check if trigger end and start are the same length
            if length(trigger_start) - length(trigger_end) == -1
                trigger_end = trigger_end(1:length(trigger_start));
            elseif length(trigger_start) - length(trigger_end) == +1
                trigger_start = trigger_start(1:length(trigger_end));
            end
            
            trigger_dur = mean(trigger_end-trigger_start);%Mean dir duration
            cut = trigger_dur*cut_factor;
            
            
            trigger_start_w = trigger_start + cut;%This creates an extra time window, so the onset 
            trigger_end_w = trigger_end - cut; %response and offset response is not counted
            
            nr_repeats_real(ii) = length(trigger_start); %The stimulus may have
            %been aborted early which will result in unenven numbers of
            %repeats for different directions, we will correct for that
            %later when taking the mean over number of repeats
            
            
            %Here we estimate the maximal number of spikes that could fall
            %into one repeat to create an 3D array of an appropriate size.
            mean_fire_rate = nanmean(diff(spiketimestamps,1),1);
            max_spikes_repeat = ceil(nanmean((trigger_dur./mean_fire_rate))*10);
            spike_nr_repeats = zeros(max_spikes_repeat,batch_info.batch_stx(bb),nr_repeats);
            spikes_nr_window = zeros(max_spikes_repeat,batch_info.batch_stx(bb),nr_repeats);
            for kk = 1:length(trigger_start)
                
               spiketime_log = logical((spiketimestamps > trigger_start(kk)).*...
                    (spiketimestamps<=trigger_end(kk)));
               spiketime_log_window = logical((spiketimestamps > trigger_start_w(kk)).*...
                    (spiketimestamps<=trigger_end_w(kk)));
                for k1 = 1:batch_info.batch_stx(bb)
                    spiketime_end = sum(spiketime_log(:,k1) ~=0,1);
                    spiketime_end_window = sum(spiketime_log_window(:,k1) ~= 0,1);
                    spike_nr_repeats(1:spiketime_end,k1,kk) =...
                        spiketimestamps(spiketime_log(:,k1),k1)-trigger_start(kk);
                    spikes_nr_window(1:spiketime_end_window,k1,kk) =...
                        spiketimestamps(spiketime_log_window(:,k1),k1)-trigger_start(kk);
                end
                
                    
                
            end
            
            dir(ii).spike_nr_repeats = spike_nr_repeats;
            dir(ii).spike_nr_repeats_window = spikes_nr_window;
        end
            %%   
       % repeat_dur = mean(trigger_end-trigger_start);

dir = dir(dir_order);
        
       
        
          
spikes_nr_dir = zeros(batch_info.batch_stx(bb),nr_dir);
spikes_nr_dir_w = zeros(batch_info.batch_stx(bb),nr_dir);


std_dir = zeros(batch_info.batch_stx(bb),nr_dir);
for cc = 1:batch_info.batch_stx(bb)
    for kk = 1:nr_dir
    


    spikes = squeeze(dir(kk).spike_nr_repeats(:,cc,:));
    spikes_w = squeeze(dir(kk).spike_nr_repeats_window(:,cc,:));


    %Here we count how many non NaN values can be found and divide
    %by number of repeats which is equal to the spikes per
    %direction
%     spikes_nr_dir(cc,kk) = nnz(spikes);
    spikes_nr_dir(cc,kk) = nnz(spikes)/nr_repeats_real(kk);
    spikes_nr_dir_w(cc,kk) = nnz(spikes_w)/nr_repeats_real(kk);
    
    
    std_dir(cc,kk) = std(sum(spikes),1);
    
    
    
    
    
    

    
    

        
   
   
    end
        
        
        
end
    
%spikes_nr_dir = spikes_nr_dir(:,dir_order); 

%Organize the data so that same orientations are coombined into same bins

spikes_nr_dir_or = zeros(batch_info.batch_stx(bb),nr_dir/2);
spikes_nr_dir_w_or = zeros(batch_info.batch_stx(bb),nr_dir/2);
for cc = 1:batch_info.batch_stx(bb)
for dd = 1:nr_dir/2
    spikes_nr_dir_or(cc,dd) = spikes_nr_dir(cc,dd) + spikes_nr_dir(cc,dd+nr_dir/2);
    spikes_nr_dir_w_or(cc,dd) = spikes_nr_dir_w(cc,dd) + spikes_nr_dir_w(dd+nr_dir/2);
    
end
end
    
    
    


%%
for cc = batch_info.batch_begins(bb):batch_info.batch_ends(bb)
    try
    %% First direction selectivity
    Data_circular(cc).spikes_deg = spikes_nr_dir(cc,:)'; %Transform to fit the toolbox
    %Calculate different statistics 
    %Functions witch one output
    Data_circular(cc).circ_mean = circ_mean(dir_radians',Data_circular(cc).spikes_deg,1);
    Data_circular(cc).circ_var = circ_var(dir_radians',Data_circular(cc).spikes_deg,diff_dir);
    Data_circular(cc).circ_rtest = circ_rtest(dir_radians',Data_circular(cc).spikes_deg,diff_dir);
    %Data_circular(cc).circ_otest = circ_otest(dir_radians',[],Data_circular(cc).spikes_deg);
    Data_circular(cc).circ_r = circ_r(dir_radians',Data_circular(cc).spikes_deg,diff_dir);
   
    %Create logical indication if test returned significant result or not
    if Data_circular(cc).circ_rtest < p_value
        Data_circular(cc).circ_rtest_sig = 1;
    end
    if Data_circular(cc).circ_otest < p_value
        Data_circular(cc).circ_otest_sig = 1;
    end
       
    
    %Functions with two outputs
    [a,b] = circ_std(dir_radians',Data_circular(cc).spikes_deg,diff_dir);
    Data_circular(cc).circ_std = [a,b];
    [a,b] = circ_skewness(dir_radians',Data_circular(cc).spikes_deg);
    Data_circular(cc).circ_skewness = [a,b];
    [a,b] = circ_kurtosis(dir_radians',Data_circular(cc).spikes_deg);
    Data_circular(cc).circ_kurtosis = [a,b];
    %%
    Data_circular(cc).confmean = circ_confmean(dir_radians',[],Data_circular(cc).spikes_deg,diff_dir');
    

%     Data_circular(cc).circ_raotest = circ_raotest(Data_circular(cc).spikes_deg);
%     Data_circular(cc).circ_vmpdf = circ_vmpdf(Data_circular(cc).spikes_deg);
    catch    
        continue
    end
    %% Next orientation selectivity
    try
     Data_circular(cc).spikes_deg_or = spikes_nr_dir_or(cc,:)'; %Transform to fit the toolbox
    %Calculate different statistics 
    %Functions witch one output
    Data_circular(cc).mean_or = circ_mean(dir_radians_or',Data_circular(cc).spikes_deg_or,1);
    Data_circular(cc).var_or = circ_var(dir_radians_or',Data_circular(cc).spikes_deg_or,diff_dir_or);
    Data_circular(cc).rtest_or = circ_rtest(dir_radians_or',Data_circular(cc).spikes_deg_or,diff_dir_or);
    %Data_circular(cc).circ_otest = circ_otest(dir_radians',[],Data_circular(cc).spikes_deg);
    Data_circular(cc).circ_r_or = circ_r(dir_radians_or',Data_circular(cc).spikes_deg_or,diff_dir_or);
   
    %Create logical indication if test returned significant result or not
    if Data_circular(cc).rtest_or < p_value
        Data_circular(cc).circ_rtest_sig_or = 1;
    end
%     if Data_circular(cc).otest_or < p_value
%         Data_circular(cc).otest_sig_or = 1;
%     end
%        
    
    %Functions with two outputs
    [a,b] = circ_std(dir_radians_or',Data_circular(cc).spikes_deg_or,diff_dir_or);
    Data_circular(cc).or_std = [a,b];
    [a,b] = circ_skewness(dir_radians_or',Data_circular(cc).spikes_deg_or);
    Data_circular(cc).or_skewness = [a,b];
    [a,b] = circ_kurtosis(dir_radians_or',Data_circular(cc).spikes_deg_or);
    Data_circular(cc).or_kurtosis = [a,b];
    %%
    Data_circular(cc).confmean = circ_confmean(dir_radians_or',[],Data_circular(cc).spikes_deg_or,diff_dir_or');
    

%     Data_circular(cc).circ_raotest = circ_raotest(Data_circular(cc).spikes_deg);
%     Data_circular(cc).circ_vmpdf = circ_vmpdf(Data_circular(cc).spikes_deg);
    catch    
        continue
    end
end


%The same for the window data

for cc = batch_info.batch_begins(bb):batch_info.batch_ends(bb)
    try
    Data_circular(cc).spikes_deg_w = spikes_nr_dir_w(cc,:)'; %Transform to fit the toolbox
    %Calculate different statistics 
    %Functions witch one output
    Data_circular(cc).circ_mean_w = circ_mean(dir_radians',Data_circular(cc).spikes_deg_w,1);
    Data_circular(cc).circ_var_w = circ_var(dir_radians',Data_circular(cc).spikes_deg_w,diff_dir);
    Data_circular(cc).circ_rtest_w = circ_rtest(dir_radians',Data_circular(cc).spikes_deg_w,diff_dir);
    %Data_circular(cc).circ_otest = circ_otest(dir_radians',[],Data_circular(cc).spikes_deg);
    Data_circular(cc).circ_r_w = circ_r(dir_radians',Data_circular(cc).spikes_deg_w,diff_dir);
   
    %Create logical indication if test returned significant result or not
    if Data_circular(cc).circ_rtest_w < p_value
        Data_circular(cc).circ_rtest_sig_w = 1;
    end
    if Data_circular(cc).circ_otest_w < p_value
        Data_circular(cc).circ_otest_sig_w = 1;
    end
       
    
    %Functions with two outputs
    [a,b] = circ_std(dir_radians',Data_circular(cc).spikes_deg_w,diff_dir);
    Data_circular(cc).circ_std_w = [a,b];
    [a,b] = circ_skewness(dir_radians',Data_circular(cc).spikes_deg_w);
    Data_circular(cc).circ_skewness_w = [a,b];
    [a,b] = circ_kurtosis(dir_radians',Data_circular(cc).spikes_deg_w);
    Data_circular(cc).circ_kurtosis_w = [a,b];
    %%
    Data_circular(cc).confmean_w = circ_confmean(dir_radians',[],Data_circular(cc).spikes_deg_w,diff_dir');
    

%     Data_circular(cc).circ_raotest = circ_raotest(Data_circular(cc).spikes_deg);
%     Data_circular(cc).circ_vmpdf = circ_vmpdf(Data_circular(cc).spikes_deg);
    catch    
        continue
    end
    
%     %% Orientation selectivity
%      try
%      Data_circular(cc).spikes_deg_w_or = spikes_nr_dir_w_or(cc,:)'; %Transform to fit the toolbox
%     %Calculate different statistics 
%     %Functions witch one output
%     Data_circular(cc).mean_w_or = circ_mean(dir_radians_or',Data_circular(cc).spikes_deg_w_or,1);
%     Data_circular(cc).var_or = circ_var(dir_radians_or',Data_circular(cc).spikes_deg_w_or,diff_dir_or);
%     Data_circular(cc).rtest_or = circ_rtest(dir_radians_or',Data_circular(cc).spikes_deg_w_or,diff_dir_or);
%     %Data_circular(cc).circ_otest = circ_otest(dir_radians',[],Data_circular(cc).spikes_deg);
%     Data_circular(cc).circ_r_or = circ_r(dir_radians_or',Data_circular(cc).spikes_deg_w_or,diff_dir_or);
%    
%     %Create logical indication if test returned significant result or not
%     if Data_circular(cc).rtest_or < p_value
%         Data_circular(cc).circ_rtest_sig_or = 1;
%     end
% %     if Data_circular(cc).otest_or < p_value
% %         Data_circular(cc).otest_sig_or = 1;
% %     end
% %        
%     
%     %Functions with two outputs
%     [a,b] = circ_std(dir_radians_or',Data_circular(cc).spikes_deg_w_or,diff_dir_or);
%     Data_circular(cc).or_std = [a,b];
%     [a,b] = circ_skewness(dir_radians_or',Data_circular(cc).spikes_deg_w_or);
%     Data_circular(cc).or_skewness = [a,b];
%     [a,b] = circ_kurtosis(dir_radians_or',Data_circular(cc).spikes_deg_w_or);
%     Data_circular(cc).or_kurtosis = [a,b];
%     %%
%     Data_circular(cc).confmean = circ_confmean(dir_radians_or',[],Data_circular(cc).spikes_deg_w_or,diff_dir_or');
%     
% 
% %     Data_circular(cc).circ_raotest = circ_raotest(Data_circular(cc).spikes_deg);
% %     Data_circular(cc).circ_vmpdf = circ_vmpdf(Data_circular(cc).spikes_deg);
%     catch    
%         continue
%     end
    
end



end
% test = Data_circular(273).spikes_deg;
% circ_plot(test,'hist',[],20,true,true,'linewidth',2,'color','r')
% 
% 
% plot the activity of the three neurons
% figure
% sub_nr = numSubplots(length(test_p));
% for j = 1:length(test_p)
%     try
%   idx = test_p(j);
%   subplot(sub_nr(1),sub_nr(2),j)
%   
%   % compute and plot mean resultant vector length and direction
%   mw = max(Data_circular(idx).spikes_deg);
%   r = circ_r(dir_radians',Data_circular(idx).spikes_deg,diff_dir) * mw;
%   r
%   phi = circ_mean(dir_radians',Data_circular(idx).spikes_deg);
%   
%   
%   hold on;
%   zm = r*exp(1i*phi');
%   plot([0 real(zm)], [0, imag(zm)],'r','linewidth',1.5)
%   % plot the tuning function of the three neurons 
%   polar([dir_radians dir_radians(1)], [Data_circular(idx).spikes_deg' Data_circular(idx).spikes_deg(1)],'k')
%  
%   
%   % draw a unit circle
%   zz = exp(1i*linspace(0, 2*pi, 101)) * mw;
%   plot(real(zz),imag(zz),'k:')
%   plot([-mw mw], [0 0], 'k:', [0 0], [-mw mw], 'k:')
%   title(num2str(idx))
% 
%   formatSubplot(gca,'ax','square','box','off','lim',[-mw mw -mw mw])
%   set(gca,'xtick',[])
%   set(gca,'ytick',[])
%     catch
%         continue
%     end
% 
% end

    

%     
%     
%     
%     
%     
%     
%     
% end
% 
% 
% 
% 
% 
%    
%     
%     %order spikes according to direction
%     spikes_nr_dir = spikes_nr_dir(:,dir_order);
%     std_dir = std_dir(:,dir_order);
%     wrong_cells = sum(spikes_nr_dir,2)<nr_dir;
% 
%     spikes_nr_dir_q = spikes_nr_dir(~wrong_cells,:);
%     std_dir_q = std_dir(~wrong_cells,:);
% 
%     spikes_nr_max = nanmax(spikes_nr_dir_q,[],2); 
%     spikes_nr_dir_norm = spikes_nr_dir_q./spikes_nr_max;
%     spikes_nr_min = nanmin(spikes_nr_dir_q,[],2);
%     ds_index = spikes_nr_min./spikes_nr_max;
%     for ss = 1:nr_dir
%         idx = dir_order(ss);
%         dir(idx).spike_nr_repeats_q = dir(ss).spike_nr_repeats(:,~wrong_cells,:);
%         
%     end
% 
% 
%     ds_filtered = find(ds_index<0.1);
% 
%     sub_nr = numSubplots(length(ds_filtered));
%     figure
%     for kk = 1:length(ds_filtered)
%         idx = ds_filtered(kk);
%         dx(kk) = subplot(sub_nr(1),sub_nr(2),kk);
%         polarplot(dir_radians,spikes_nr_dir_q(idx,:))
%     end
%     
%     
% 
% 
% 
% %% asdaw
% 
% 
% 
% for ii = 1:20
% idx = ii;
% figure
% sub_nr = numSubplots(nr_dir);
% for kk = 1:nr_dir
%     
%    %Collect all spikes
%    ax(kk) = subplot(sub_nr(1),sub_nr(2),kk);
%    for cc = 1
%    
%    spikes = squeeze(dir(kk).spike_nr_repeats(:,idx,:));
%    spikes(spikes == 0) = NaN;
%    yvalue = ones(1,length(spikes));
%    
%    for rr = 1:nr_repeats
%        last_idx = lastNaN(spikes(:,rr),1);
%        scatter(spikes(1:last_idx,rr),(yvalue(1:last_idx)*0.1*rr)+(cc-1)*yvalue(1:last_idx),'k','.');
%        hold on
%    end
%    end
%     ylim([0,cc+0.1*nr_repeats])
% end
%   
%    linkaxes(ax,'x');
% end
%     

    out = sf_organizer(stim_idx,savepath,'variable_name','Data_circular','variable',Data_circular);
    out = sf_organizer(stim_idx,savepath,'variable_name','dir','variable',dir);
end
   

end


                
                
                

            
            
            
        
        
        
    
    
    
    
    
    
    
    
    
    
    
    


    










% for ii = 1:nr_dir
%            %check which index is the right
%             trigger_idx = repeat_idx(ii);
%             end_idx = repeat_idx(ii+1);
%             %find location of similar events
%             event_idx = repeat_idx == trigger_idx;
%             event_end = repeat_idx == end_idx;
%             trigger_start = locs(event_idx(1:length(locs)))/SamplingFrequency;
%             trigger_end = locs(event_end(1:length(locs)))/SamplingFrequency;
%             %I guess the last trigger end based on the mean intervalt, but I think thats ok here...
%             if length(trigger_start)-length(trigger_end) == 1
%                 trigger_end(end+1) = trigger_end(end)+ mean(diff(trigger_end));
%             end
%             
%             
%             %find spikes within that range
%             spike_nr_repeats = NaN(length(trigger_times),batch_info.batch_stx);
%             for kk = 1:length(trigger_start)
%                 
%                 spike_nr_repeats(kk,:) = sum((spiketimestamps > trigger_start(kk)).*...
%                     (spiketimestamps<trigger_end(kk))~=0,1);
%             end
%             spike_nr_dir(ii,:) = nanmean(spike_nr_repeats,1);
%             std_dir(ii,:) = nanstd(spike_nr_repeats,1);
%         end