function STA = STA_fn(stimulus_arr,p,STE,mean_raw_stim_arr)
% Calculate the STA

STA_sum_arr  = zeros(p.stim_rows,p.stim_columns,p.Num_STE_bins);

if p.Time_Choice == 1 % stimulus frames
        
    for i = 1:STE.Num_unique_STEs
        
        % Choose Num_STE_bins stim frames prior to frame in which spike occurs
        STA_loop_arr = stimulus_arr(:,:,STE.STE_index_unique_arr(i)-p.Num_STE_bins:STE.STE_index_unique_arr(i)-1);
        
        if     p.STA_Choice  == 1 % don't subtract average stim;
            STA_sum_arr = STA_sum_arr + STE.STE_rep_vec(i)*STA_loop_arr;
        elseif p.STA_Choice  == 2 % subtract average stim;
            STA_sum_arr = STA_sum_arr + STE.STE_rep_vec(i)*(STA_loop_arr - mean_raw_stim_arr);
        end
        
    end
    
elseif p.Time_Choice == 2 % define own time grid
    
    for i = 1:STE.Num_unique_STEs
        
        STA_loop_arr = stimulus_arr(:,:,STE.STE_index_unique_arr(i,:));
        
        if p.STA_Choice  == 1     % don't subtract average stim;
            STA_sum_arr  = STA_sum_arr + STE.STE_rep_vec(i)*STA_loop_arr;
        elseif p.STA_Choice  == 2 % subtract average stim;
            STA_sum_arr  = STA_sum_arr + STE.STE_rep_vec(i)*(STA_loop_arr - mean_raw_stim_arr);
        end
        
    end
    
end

STA = STA_sum_arr/p.length_spike_times;