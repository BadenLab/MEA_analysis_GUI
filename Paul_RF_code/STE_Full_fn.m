function STE_Full = STE_Full_fn(stimulus_arr,trig_times_vec,spike_times_vec,p)
% Calculate the Full Spike-Triggered Ensemble (STE) array

STE_Full = NaN(p.stim_rows,p.stim_columns,p.Num_STE_bins,p.length_spike_times);
%Paul    = squeeze(STE_Full(:,:,:,1));

if p.Time_Choice == 1 % stimulus frames
    
    for i = 1:p.length_spike_times
        
        STE_index_loop    = sum((trig_times_vec-spike_times_vec(i))<0);
        
        STE_Full(:,:,:,i) = stimulus_arr(:,:,STE_index_loop-p.Num_STE_bins:STE_index_loop-1);
        
    end
    
    %%% Store only unique spike triggering stimuli and a record of the number of
    %%% repeats of each
    % 'stable' - keeps rows in same order
    %     [STE.STE_index_unique_arr,~,ic] = unique(STE_index_vec,'stable');
    %     STE.STE_rep_vec             = hist(ic,unique(ic))';
    %     STE.Num_unique_STEs         = length(STE.STE_rep_vec);
    
elseif p.Time_Choice == 2 % define own time grid
    
    STE_index_vec = NaN(p.Num_STE_bins,1);
    
    for i = 1:p.length_spike_times
        
        loop_sample_times_vec = spike_times_vec(i) + p.stim_timesample_vec;
        
        for j = 1:p.Num_STE_bins
            
            % Find which stimulis frame was active at each of the sample times
            % before the spike. Indices are recorded (use sum to find index as
            % number of triggers fired by the j'th time).
            STE_index_vec(j) = sum((trig_times_vec-loop_sample_times_vec(j))<0);
            
        end
        
        STE_Full(:,:,:,i) = stimulus_arr(:,:,STE_index_vec);
        
    end
    
    %%% Store only unique spike triggering stimuli and a record of the number of
    %%% repeats of each
    % 'stable' - keeps rows in same order
    % 'rows'   - keeps full rows in output
%     [STE.STE_index_unique_arr,~,ic] = unique(STE_index_mat,'stable','rows');
%     STE.STE_rep_vec             = hist(ic,unique(ic))';
%     STE.Num_unique_STEs         = length(STE.STE_rep_vec);
    
end

