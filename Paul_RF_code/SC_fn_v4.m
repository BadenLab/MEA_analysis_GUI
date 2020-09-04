function [Stixel_covar,Pixel_covar] = SC_fn_v4(stimulus_arr,STE,meanSpikeTrigStim,p)
% Calculate the SC array

%STE.STE_index_unique_arr
%STE.STE_rep_vec 
%STE.Num_unique_STEs

Stixel_covar = zeros(p.stim_rows,p.stim_columns,p.Num_STE_bins);%NaN(p.stim_rows,p.stim_columns,p.Num_STE_bins)
%Pixel_covar = NaN(p.stim_rows,p.stim_columns);

for i = 1:p.Num_STE_bins
    
    for j = 1:p.stim_columns
        
        for k = 1:p.stim_rows
            
            %Stixel_vec_loop = NaN(1,STE.Num_unique_STEs);
            if p.Time_Choice == 1 % stimulus frames
                 
                for l = 1:STE.Num_unique_STEs
                    %Stixel_vec_loop(l) = stimulus_arr(k,j,STE.STE_index_unique_arr(l)-i);
                    Stixel_covar(k,j,i) = Stixel_covar(k,j,i) + STE.STE_rep_vec(l)*(stimulus_arr(k,j,STE.STE_index_unique_arr(l)-p.Num_STE_bins+i-1) - meanSpikeTrigStim(k,j,i))^2;
                end
            elseif p.Time_Choice == 2 % define own time grid
                
                
                
                for l = 1:STE.Num_unique_STEs
                    %Stixel_vec_loop(l) = stimulus_arr(k,j,STE.STE_index_unique_arr(l,i));
                    Stixel_covar(k,j,i) = Stixel_covar(k,j,i) + STE.STE_rep_vec(l)*(stimulus_arr(k,j,STE.STE_index_unique_arr(l,i)) - meanSpikeTrigStim(k,j,i))^2;
                end
            end
            
            %Stixel_covar(k,j,i) = var(Stixel_vec_loop); % or cov
            
        end
        
    end
    
end

Stixel_covar = Stixel_covar/(p.length_spike_times - 1);

Pixel_covar = min(Stixel_covar,[],3); % min(Stixel_covar,[],3) mean(Stixel_covar,3) max(Stixel_covar,[],3)