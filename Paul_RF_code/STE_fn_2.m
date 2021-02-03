%function [STE_index_unique_arr,STE_rep_vec,Num_unique_STEs] = STE_fn(trig_times_vec,spike_times_vec,p)
function STE = STE_fn_2(trig_times_vec,spike_times_vec,p)
% Calculate the Spike-Triggered Ensemble (STE)
% Modified from STE_fn on 22,09,2020 to deal with repeated chuncks of
% frozen noise w/o gaps periodicity and problems with spike after final
% stim frame

if p.Time_Choice == 1 % stimulus frames
    
    STE_index_vec = NaN(p.length_spike_times,1);
    
    
    if p.Gap_Ind == 1 %  w/o gaps PAR Mod 17,09,2020 (whole if statement, replacing first for loop below)
        
        for i = 1:p.length_spike_times
            STE_index_vec(i) = sum((trig_times_vec-spike_times_vec(i))<0);
        end
        
    else % p.Gap_Ind == 2 with gaps
        
        trig_times_vec_mod = [trig_times_vec,trig_times_vec(end)+p.stim_int]; % (trig_times_vec is a row vec)
        
        for i = 1:p.length_spike_times
            STE_index_vec(i) = sum((trig_times_vec_mod-spike_times_vec(i))<0);
        end
        
    end
    
    %%% Store only unique spike triggering stimuli and a record of the number of
    %%% repeats of each
    % 'stable' - keeps rows in same order
    [STE.STE_index_unique_arr,~,ic] = unique(STE_index_vec,'stable');
    STE.STE_rep_vec             = hist(ic,unique(ic))';
    STE.Num_unique_STEs         = length(STE.STE_rep_vec);
    
elseif p.Time_Choice == 2 % define own time grid
    
    STE_index_mat = NaN(p.length_spike_times,p.Num_STE_bins);
    
    if p.Gap_Ind == 1 %  w/o gaps PAR Mod 17,09,2020 (whole if statement)
        
        for i = 1:p.length_spike_times
            
            loop_spike_time = spike_times_vec(i);
            
            loop_sample_times_vec = loop_spike_time + p.stim_timesample_vec;
            
            loop_sample_times_vec_mod = mod(loop_sample_times_vec,(trig_times_vec(end)+p.stim_int)); % PAR Mod 17,09,2020
            
            for j = 1:p.Num_STE_bins
                
                % Find which stimulis frame was active at each of the sample times
                % before the spike. Indices are recorded (use sum to find index as
                % number of triggers fired by the j'th time).
                STE_index_mat(i,j) = sum((trig_times_vec-loop_sample_times_vec_mod(j))<0); % PAR Mod 17,09,2020 (loop_sample_times_vec --> loop_sample_times_vec_mod)
                
            end
            
        end
        
    else % p.Gap_Ind == 2 with gaps
        
        trig_times_vec_mod = [trig_times_vec,trig_times_vec(end)+p.stim_int]; % (trig_times_vec is a row vec)
        
        for i = 1:p.length_spike_times
            
            loop_spike_time = spike_times_vec(i);
            
            loop_sample_times_vec = loop_spike_time + p.stim_timesample_vec;
            
            for j = 1:p.Num_STE_bins
                
                % Find which stimulis frame was active at each of the sample times
                % before the spike. Indices are recorded (use sum to find index as
                % number of triggers fired by the j'th time).
                STE_index_mat(i,j) = sum((trig_times_vec_mod -loop_sample_times_vec(j))<0);
                
            end
            
        end
        
        % PAR Mod 03,02,2021 (occasionally the sample comes from a time slightly before the first
        % trigger time due to hard-to-avoid minor inaccuracies in reading the trigger channel,
        % this error crops up in the spike mapping stage, but such spikes should be in the first stim window
        % given the spike removal stage, so this corrects the error)
        STE_index_mat(STE_index_mat==0)=1;
        
    end
    
    %%% Store only unique spike triggering stimuli and a record of the number of
    %%% repeats of each
    % 'stable' - keeps rows in same order
    % 'rows'   - keeps full rows in output
    [STE.STE_index_unique_arr,~,ic] = unique(STE_index_mat,'stable','rows');
    STE.STE_rep_vec             = hist(ic,unique(ic))';
    STE.Num_unique_STEs         = length(STE.STE_rep_vec);
    
end

