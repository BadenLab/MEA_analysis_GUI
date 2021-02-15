function STA = meanSTS_fn_2(stimulus_arr,p,STE)
% Calculate the mean spike triggered stimulus (w/o removing the mean raw stim)

STA_sum_arr  = zeros(p.stim_rows,p.stim_columns,p.Num_STE_bins);

if p.Time_Choice == 1 % stimulus frames
    
    
    if p.Gap_Ind == 1 %  w/o gaps PAR Mod 17,09,2020 (whole if statement)
        
        for i = 1:STE.Num_unique_STEs
            
            loop_sample_indices_vec = mod(STE.STE_index_unique_arr(i)-p.Num_STE_bins:STE.STE_index_unique_arr(i)-1,p.noise_length);
            loop_sample_indices_vec(loop_sample_indices_vec==0) = p.noise_length;
            % Choose Num_STE_bins stim frames prior to frame in which spike occurs
            STA_loop_arr = stimulus_arr(:,:,loop_sample_indices_vec);
            
            STA_sum_arr = STA_sum_arr + STE.STE_rep_vec(i)*STA_loop_arr;
            
        end
        
    else % p.Gap_Ind == 2 with gaps
        
        for i = 1:STE.Num_unique_STEs
            
            % Choose Num_STE_bins stim frames prior to frame in which spike occurs
            STA_loop_arr = stimulus_arr(:,:,STE.STE_index_unique_arr(i)-p.Num_STE_bins:STE.STE_index_unique_arr(i)-1);
            
            STA_sum_arr = STA_sum_arr + STE.STE_rep_vec(i)*STA_loop_arr;
            
        end
        
    end
    
    
elseif p.Time_Choice == 2 % define own time grid
    
    
    if p.Gap_Ind == 1 %  w/o gaps PAR Mod 17,09,2020 (whole if statement)
        
        for i = 1:STE.Num_unique_STEs
            
            STA_loop_arr = stimulus_arr(:,:,STE.STE_index_unique_arr(i,:));
            
            STA_sum_arr  = STA_sum_arr + STE.STE_rep_vec(i)*STA_loop_arr;
            
        end
        
    else % p.Gap_Ind == 2 with gaps
        
        for i = 1:STE.Num_unique_STEs
            
            STA_loop_arr = NaN(p.stim_rows,p.stim_columns,p.Num_STE_bins); % Preallocate
            for j = 1:p.Num_STE_bins
                if STE.STE_index_unique_arr(i,j) <= p.noise_length
                    STA_loop_arr(:,:,j) = stimulus_arr(:,:,STE.STE_index_unique_arr(i,j));
                else % STE.STE_index_unique_arr(i,j) == p.noise_length + 1;
                    STA_loop_arr(:,:,j) = zeros(p.stim_rows,p.stim_columns);
                end
            end
            
            STA_sum_arr  = STA_sum_arr + STE.STE_rep_vec(i)*STA_loop_arr;
            
        end
        
    end
    
end

STA = STA_sum_arr/p.length_spike_times;