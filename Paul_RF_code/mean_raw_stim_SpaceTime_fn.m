function mean_raw_stim_arr = mean_raw_stim_SpaceTime_fn(stimulus_arr,trig_times_vec,p)
% Calculate the mean of the raw stimuli

mean_raw_stim_arr_loop = zeros(p.stim_rows,p.stim_columns,p.Num_STE_bins);

if p.Time_Choice == 1 || p.Mean_Stim_Choice  == 1 % stimulus frames
    
    for i = 1:p.Num_Raw_Stim
        mean_raw_stim_arr_loop = mean_raw_stim_arr_loop + stimulus_arr(:,:,i:p.Num_STE_bins+i-1);
    end
    
elseif p.Time_Choice == 2 && p.Mean_Stim_Choice  == 2 % define own time grid
    
    for i = 1:p.Num_Raw_Stim
        
        loop_sample_times_vec = p.Num_Raw_Stim_t_vec(i) + p.stim_timesample_vec;
        loop_index_vec        = zeros(p.Num_STE_bins,1);
        
        % Find which stimulis frame was active at each of the sample times
        % Indices are recorded (use sum to find index as
        % number of triggers fired by the j'th time).
        for j = 1:p.Num_STE_bins
            loop_index_vec(j) = sum((trig_times_vec-loop_sample_times_vec(j))<=1e-10);
            % was <0, make <= 0 to allow for equality case which arrises here.
            % make 1e-10 since can get values like 4.4409*1e-16 for
            % difference beteween first entries of 'trig_times_vec' and
            % 'loop_sample_times_vec' and hence loop_index_vec(1) = 0.
        end
        
        mean_raw_stim_arr_loop = mean_raw_stim_arr_loop + stimulus_arr(:,:,loop_index_vec);
        
    end
    
end

mean_raw_stim_arr = mean_raw_stim_arr_loop/p.Num_Raw_Stim;