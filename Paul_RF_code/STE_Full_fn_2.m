function STE_Full = STE_Full_fn_2(stimulus_arr,trig_times_vec,spike_times_vec,p)
% Calculate the Full Spike-Triggered Ensemble (STE) array
% Modified from STE_Full_fn on 23,09,2020 to deal with repeated chuncks of
% frozen noise w/o gaps periodicity and problems with spike after final
% stim frame

STE_Full = NaN(p.stim_rows,p.stim_columns,p.Num_STE_bins,p.length_spike_times);
%Paul    = squeeze(STE_Full(:,:,:,1));

if p.Time_Choice == 1 % stimulus frames
    
    if p.Gap_Ind == 1 %  w/o gaps PAR Mod 17,09,2020 (whole if statement)
        
        for i = 1:p.length_spike_times
            
            STE_index_loop    = sum((trig_times_vec-spike_times_vec(i))<0);
            loop_sample_indices_vec = mod(STE_index_loop-p.Num_STE_bins:STE_index_loop-1,p.noise_length); % PAR Mod 17,09,2020
            loop_sample_indices_vec(loop_sample_indices_vec==0) = p.noise_length;
            STE_Full(:,:,:,i) = stimulus_arr(:,:,loop_sample_indices_vec);
        end
        
    else % p.Gap_Ind == 2 with gaps
        
        trig_times_vec_mod = [trig_times_vec,trig_times_vec(end)+p.stim_int]; % (trig_times_vec is a row vec)
        
        for i = 1:p.length_spike_times
            STE_index_loop    = sum((trig_times_vec_mod-spike_times_vec(i))<0);
            STE_Full(:,:,:,i) = stimulus_arr(:,:,STE_index_loop-p.Num_STE_bins:STE_index_loop-1);
        end
        
    end
    
    
elseif p.Time_Choice == 2 % define own time grid
    
    STE_index_vec = NaN(p.Num_STE_bins,1);
    
    if p.Gap_Ind == 1 %  w/o gaps PAR Mod 17,09,2020 (whole if statement)
        
        for i = 1:p.length_spike_times
            
            loop_sample_times_vec = spike_times_vec(i) + p.stim_timesample_vec;
            
            loop_sample_times_vec_mod = mod(loop_sample_times_vec,(trig_times_vec(end)+p.stim_int)); % PAR Mod 17,09,2020
            
            for j = 1:p.Num_STE_bins
                
                % Find which stimulis frame was active at each of the sample times
                % before the spike. Indices are recorded (use sum to find index as
                % number of triggers fired by the j'th time).
                STE_index_vec(j) = sum((trig_times_vec-loop_sample_times_vec_mod(j))<0); % PAR Mod 17,09,2020 (loop_sample_times_vec --> loop_sample_times_vec_mod)
                
            end
            
            STE_Full(:,:,:,i) = stimulus_arr(:,:,STE_index_vec);
            
        end
        
    else % p.Gap_Ind == 2 with gaps
        
        trig_times_vec_mod = [trig_times_vec,trig_times_vec(end)+p.stim_int]; % (trig_times_vec is a row vec) PAR Mod 17,09,2020
        
        for i = 1:p.length_spike_times
            
            loop_sample_times_vec = spike_times_vec(i) + p.stim_timesample_vec;
            
            for j = 1:p.Num_STE_bins
                
                % Find which stimulis frame was active at each of the sample times
                % before the spike. Indices are recorded (use sum to find index as
                % number of triggers fired by the j'th time).
                STE_index = sum((trig_times_vec_mod-loop_sample_times_vec(j))<0); % was STE_index_vec(j) PAR Mod 17,09,2020
                
                % PAR Mod 03,02,2021 (occasionally the sample comes from a time slightly before the first
                % trigger time due to hard-to-avoid minor inaccuracies in reading the trigger channel,
                % this error crops up in the spike mapping stage, but such spikes should be in the first stim window
                % given the spike removal stage, so this corrects the error)
                if STE_index==0 
                    STE_index=1;
                end
                
                if STE_index == p.noise_length + 1 % PAR Mod 17,09,2020 (To deal with cases where we sample from a time where there was no stimulus)
                    STE_Full(:,:,j,i) = zeros(p.stim_rows,p.stim_columns);
                else
                    STE_Full(:,:,j,i) = stimulus_arr(:,:,STE_index);
                end
                
            end
            
        end
        
    end
    
end

