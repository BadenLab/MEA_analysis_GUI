%function [STE_index_unique_arr,STE_rep_vec,Num_unique_STEs] = STE_fn(trig_times_vec,spike_times_vec,p)
function STE = STE_fn(trig_times_vec,spike_times_vec,p)
% Calculate the Spike-Triggered Ensemble (STE)

if p.Time_Choice == 1 % stimulus frames
    
    STE_index_vec = NaN(p.length_spike_times,1);
    
    for i = 1:p.length_spike_times
        
        STE_index_vec(i) = sum((trig_times_vec-spike_times_vec(i))<0);
        
    end
    
    %%% Store only unique spike triggering stimuli and a record of the number of
    %%% repeats of each
    % 'stable' - keeps rows in same order
    [STE.STE_index_unique_arr,~,ic] = unique(STE_index_vec,'stable');
    STE.STE_rep_vec             = hist(ic,unique(ic))';
    STE.Num_unique_STEs         = length(STE.STE_rep_vec);
    
elseif p.Time_Choice == 2 % define own time grid
    
    STE_index_mat = NaN(p.length_spike_times,p.Num_STE_bins);
    
    for i = 1:p.length_spike_times
        
        loop_spike_time = spike_times_vec(i);
        
        loop_sample_times_vec = loop_spike_time + p.stim_timesample_vec;
        
        for j = 1:p.Num_STE_bins
            
            % Find which stimulis frame was active at each of the sample times
            % before the spike. Indices are recorded (use sum to find index as
            % number of triggers fired by the j'th time).
            STE_index_mat(i,j) = sum((trig_times_vec-loop_sample_times_vec(j))<0);
            
        end
        
    end
    
    %%% Store only unique spike triggering stimuli and a record of the number of
    %%% repeats of each
    % 'stable' - keeps rows in same order
    % 'rows'   - keeps full rows in output
    [STE.STE_index_unique_arr,~,ic] = unique(STE_index_mat,'stable','rows');
    STE.STE_rep_vec             = hist(ic,unique(ic))';
    STE.Num_unique_STEs         = length(STE.STE_rep_vec);
    
end

