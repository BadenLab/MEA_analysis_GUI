function [Stixel_covar,Pixel_covar] = SC_fn_v5(stimulus_arr,STE,meanSpikeTrigStim,p)
% Calculate the SC array
% Modified from SC_fn_v4 on 25,09,2020 to deal with repeated chunks of
% frozen noise w/o gaps periodicity and problems with spike after final
% stim frame

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
                 
                
                if p.Gap_Ind == 1 %  w/o gaps PAR Mod 17,09,2020 (whole if statement)
                    
                    for l = 1:STE.Num_unique_STEs
                        loop_sample_index = mod(STE.STE_index_unique_arr(l)-p.Num_STE_bins+i-1,p.noise_length);
                        loop_sample_index(loop_sample_index==0) = p.noise_length;
                        Stixel_covar(k,j,i) = Stixel_covar(k,j,i) + STE.STE_rep_vec(l)*(stimulus_arr(k,j,loop_sample_index) - meanSpikeTrigStim(k,j,i))^2;
                    end
                    
                else % p.Gap_Ind == 2 with gaps
                    
                    for l = 1:STE.Num_unique_STEs
                        Stixel_covar(k,j,i) = Stixel_covar(k,j,i) + STE.STE_rep_vec(l)*(stimulus_arr(k,j,STE.STE_index_unique_arr(l)-p.Num_STE_bins+i-1) - meanSpikeTrigStim(k,j,i))^2;
                    end
                    
                end
                
                
            elseif p.Time_Choice == 2 % define own time grid
                
                
                if p.Gap_Ind == 1 %  w/o gaps PAR Mod 17,09,2020 (whole if statement)
                    
                    for l = 1:STE.Num_unique_STEs
                        Stixel_covar(k,j,i) = Stixel_covar(k,j,i) + STE.STE_rep_vec(l)*(stimulus_arr(k,j,STE.STE_index_unique_arr(l,i)) - meanSpikeTrigStim(k,j,i))^2;
                    end
                    
                else % p.Gap_Ind == 2 with gaps
                    
                    for l = 1:STE.Num_unique_STEs
                        
                        if STE.STE_index_unique_arr(l,i) <= p.noise_length
                            Stim_loop = stimulus_arr(k,j,STE.STE_index_unique_arr(l,i));
                        else % STE.STE_index_unique_arr(l,i) == p.noise_length + 1;
                            Stim_loop = 0;
                        end
                        
                        Stixel_covar(k,j,i) = Stixel_covar(k,j,i) + STE.STE_rep_vec(l)*(Stim_loop - meanSpikeTrigStim(k,j,i))^2;
                        
                    end
                    
                end
                
                
            end
            
            %Stixel_covar(k,j,i) = var(Stixel_vec_loop); % or cov
            
        end
        
    end
    
end

Stixel_covar = Stixel_covar/(p.length_spike_times - 1);

Pixel_covar = min(Stixel_covar,[],3); % min(Stixel_covar,[],3) mean(Stixel_covar,3) max(Stixel_covar,[],3)