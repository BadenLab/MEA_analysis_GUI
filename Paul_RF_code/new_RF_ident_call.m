%% Call RF Identification Function
% MEA data: stimulus_array := 40 rows, 40 columns,6000 frames, 4 colours
a = 1;
%RF_Ident = cell(stx,1);

if p.RF_Ident_Meth_vec(1) == 1 % STA-SD method
    if p.STA_Choice  == 2 % subtract average stim
        p.mean_raw_stim_arr = NaN(p.stim_rows,p.stim_columns,p.Num_STE_bins,p.Spectral_Dim);
        for i = 1:p.Spectral_Dim
            p.mean_raw_stim_arr(:,:,:,i) = mean_raw_stim_SpaceTime_fn_v2(stimulus_arr(:,:,:,i),trig_times_vec_trunc,p); % PAR Mod 27,08,2020 was 'mean_raw_stim_SpaceTime_fn' and 'trig_times_vec'
        end
      
    end
end

tic;
if Parpool == 1 % Parpool on %@mars: This doesnt work if a parpool is already active
    
    %parpool(Num_Cores);
    
    parfor i = 1:stx % for/parfor
        i
        settings = add_info;  
        p_par = p;
        trig_times_vec_par = trig_times_vec; %@mars: added this declaration
        %to avoid uneccessary communication between the workers in the parloop.
        spike_times_vec_loop      = spiketimestamps(:,i);
        %Check if all values are NAN
        test_for_NaN = find(all(isnan(spike_times_vec_loop),1));
        if test_for_NaN
            RF_Ident(i).RF_results = NaN;
            RF_Ident(i).cell_idx = cell_indices_temp(i);
             
            RF_Ident(i).add_info = settings;
            
            continue
        end
        
        spike_times_vec_loop(isnan(spike_times_vec_loop)) = []; % spiketimestamps(~isnan(spiketimestamps(:,Cell_Choice)),Cell_Choice);
        spike_times_vec_loop      = spike_times_vec_loop(spike_times_vec_loop>(first_trig_time + min_start_int));
        %spike_times_vec_loop      = spike_times_vec_loop(spike_times_vec_loop<(last_trig_time + min_end_int));
        if (nr_trig_last_repeat >= p.Num_STE_bins + 1) || (nr_trig_last_repeat == 0)  % PAR Mod 17,09,2020 (whole if statement, replaces 'spike_times_vec_loop      = spike_times_vec_loop(spike_times_vec_loop<(last_trig_time + min_end_int));')
            spike_times_vec_loop      = spike_times_vec_loop(spike_times_vec_loop<(last_trig_time + min_end_int));
        else % Num_Trig_Final_NoiseChunk < p.Num_STE_bins + 1
            spike_times_vec_loop      = spike_times_vec_loop(spike_times_vec_loop<(trig_times_vec(p.stim_frames*(p.Num_FNoise_rep_ceil-1)) + min_end_int));
        end
        
        
        
        
        
        if p_par.Num_trigs > p_par.stim_frames % Frozen noise repeated case % PAR Mod 27,08,2020 (whole if statement)
            
            if sum(inter_noise_int_vec <= max_trig_int) == p_par.Num_FNoise_rep_ceil-1 % If the frozen noise repeates w/o a gap
                
                % No need to remove further spikes in this case
                
                for j = 2:p_par.Num_FNoise_rep_ceil % map spikes in repeated stimulus chunks to original chunk
                    
                    if j==p_par.Num_FNoise_rep_ceil
                        frame_end_loop = p_par.Num_trigs - p_par.stim_frames*(p_par.Num_FNoise_rep_ceil-1);
                    else
                        frame_end_loop = p_par.stim_frames;
                    end
                    
                    for k = 1:frame_end_loop
                        if (k < frame_end_loop) || (j < p_par.Num_FNoise_rep_ceil)
                            indices_loop = find((spike_times_vec_loop>trig_times_vec_par((j-1)*p_par.stim_frames + k))&&(spike_times_vec_loop<trig_times_vec_par((j-1)*p_par.stim_frames + k + 1)));
                            a_loop      = trig_times_vec_par(k);
                            b_loop      = trig_times_vec_par(k+1);
                            a_dash_loop = trig_times_vec_par((j-1)*p_par.stim_frames + k);
                            b_dash_loop = trig_times_vec_par((j-1)*p_par.stim_frames + k + 1);
                            %spike_times_vec_loop(indices_loop) = a_loop + (spike_times_vec_loop(indices_loop) - a_dash_loop)*(b_loop - a_loop)/(b_dash_loop - a_dash_loop);
                        else % (k == frame_end_loop) && (j == p.Num_FNoise_rep_ceil) % No end time is given for the last frame so have to specify differently
                            indices_loop = find((spike_times_vec_loop>trig_times_vec_par((j-1)*p_par.stim_frames + k))&&(spike_times_vec_loop<trig_times_vec_par((j-1)*p_par.stim_frames + k) + min_end_int));
                            a_loop      = trig_times_vec_par(k);
                            b_loop      = trig_times_vec_par(k+1);
                            a_dash_loop = trig_times_vec_par((j-1)*p_par.stim_frames + k);
                            b_dash_loop = trig_times_vec_par((j-1)*p_par.stim_frames + k) + min_end_int;
                            %spike_times_vec_loop(indices_loop) = a_loop + (spike_times_vec_loop(indices_loop) - a_dash_loop)*(b_loop - a_loop)/(b_dash_loop - a_dash_loop);
                        end
                        spike_times_vec_loop(indices_loop) = a_loop + (spike_times_vec_loop(indices_loop) - a_dash_loop)*(b_loop - a_loop)/(b_dash_loop - a_dash_loop); % PAR Mod 17,09,2020 (moved from inside if statement above)
                    end
                    
                end
                
                
                 % Map spike in interval after last stimulus window % PAR Mod 17,09,2020 (everything from here to the next 'else' is new)
                indices_loop = find((spike_times_vec_loop>trig_times_vec(p.Num_trigs) + stim_int)...
                                   &(spike_times_vec_loop<trig_times_vec(p.Num_trigs) + min_end_int));
                a_loop      = trig_times_vec(1);
                b_loop      = trig_times_vec(2);
                a_dash_loop = trig_times_vec(p.Num_trigs) + stim_int;
                b_dash_loop = trig_times_vec(p.Num_trigs) + min_end_int;
                spike_times_vec_loop(indices_loop) = a_loop + (spike_times_vec_loop(indices_loop) - a_dash_loop)*(b_loop - a_loop)/(b_dash_loop - a_dash_loop);
                
            else % If there are gaps between (any) frozen noise chuncks
                
                if (nr_trig_last_repeat >= p.Num_STE_bins + 1) || (nr_trig_last_repeat == 0) % PAR Mod 17,09,2020 (if statement is new, but first for loop it contains was here originally)
                    for j = 1:p.Num_FNoise_rep_ceil-1 % remove spikes with stimulus windows that fall between the noise chunks
                        spike_times_vec_loop((spike_times_vec_loop>(trig_times_vec(j*p.stim_frames) + min_end_int))&(spike_times_vec_loop<trig_times_vec(j*p.stim_frames+p.Num_STE_bins+1))) = [];
                    end
                    
                elseif p.Num_FNoise_rep_ceil>2 % still need to remove spike between earlier chunks, priveded there were more than 2 originally
                    for j = 1:p.Num_FNoise_rep_ceil-2 % remove spikes with stimulus windows that fall between the noise chunks
                        spike_times_vec_loop((spike_times_vec_loop>(trig_times_vec(j*p.stim_frames) + min_end_int))&(spike_times_vec_loop<trig_times_vec(j*p.stim_frames+p.Num_STE_bins+1))) = [];
                    end
                end
                
                for j = 2:p_par.Num_FNoise_rep_ceil % map spikes in repeated stimulus chunks to original chunk
                    
                    if j==p_par.Num_FNoise_rep_ceil
                        frame_end_loop = p_par.Num_trigs - p_par.stim_frames*(p_par.Num_FNoise_rep_ceil-1);
                    else
                        frame_end_loop = p_par.stim_frames;
                    end
                    
                    for k = 1:frame_end_loop+1
                        if k<frame_end_loop % frame_end_loop<p.stim_frames
                            
                            indices_loop = find((spike_times_vec_loop>=trig_times_vec((j-1)*p_par.stim_frames + k)).*(spike_times_vec_loop<trig_times_vec((j-1)*p_par.stim_frames + k + 1)));
                            a_loop      = trig_times_vec(k);
                            b_loop      = trig_times_vec(k+1);
                            a_dash_loop = trig_times_vec((j-1)*p_par.stim_frames + k);
                            b_dash_loop = trig_times_vec((j-1)*p_par.stim_frames + k + 1);
                            %spike_times_vec_loop(indices_loop) = a_loop + (spike_times_vec_loop(indices_loop) - a_dash_loop)*(b_loop - a_loop)/(b_dash_loop - a_dash_loop);
                       elseif  k==frame_end_loop % No end time is given for the last frame so have to specify differently % PAR Mod 17,09,2020 (was 'else')
                            indices_loop = find((spike_times_vec_loop>trig_times_vec((j-1)*p.stim_frames + k))&(spike_times_vec_loop<trig_times_vec((j-1)*p.stim_frames + k) + stim_int)); % PAR Mod 17,09,2020 (& not &&)
                            a_loop      = trig_times_vec(k);
                            b_loop      = trig_times_vec(k) + stim_int; % PAR Mod 17,09,2020 (trig_times_vec(k) + stim_int not trig_times_vec(k+1))
                            a_dash_loop = trig_times_vec((j-1)*p_par.stim_frames + k);
                            b_dash_loop = trig_times_vec((j-1)*p_par.stim_frames + k) + min_end_int;
                            %spike_times_vec_loop(indices_loop) = a_loop + (spike_times_vec_loop(indices_loop) - a_dash_loop)*(b_loop - a_loop)/(b_dash_loop - a_dash_loop);
                        
                        else % k==frame_end_loop+1 % PAR Mod 17,09,2020 (everything between this else statement and the end of the if statement is new)
                            indices_loop = find((spike_times_vec_loop>trig_times_vec((j-1)*p.stim_frames + k-1) + stim_int)...
                                               &(spike_times_vec_loop<trig_times_vec((j-1)*p.stim_frames + k-1) + min_end_int));
                            if frame_end_loop == p.stim_frames
                                a_loop      = trig_times_vec(p.stim_frames) + stim_int;
                                b_loop      = trig_times_vec(p.stim_frames) + min_end_int;
                                a_dash_loop = trig_times_vec(j*p.stim_frames) + stim_int;
                                b_dash_loop = trig_times_vec(j*p.stim_frames) + min_end_int;
                            else %frame_end_loop < p.stim_frames
                                a_loop      = trig_times_vec(k);
                                b_loop      = trig_times_vec(k+1);
                                a_dash_loop = trig_times_vec((j-1)*p.stim_frames + k-1) + stim_int;
                                b_dash_loop = trig_times_vec((j-1)*p.stim_frames + k-1) + min_end_int;
                            end
                             spike_times_vec_loop(indices_loop) = a_loop + (spike_times_vec_loop(indices_loop) - a_dash_loop)*(b_loop - a_loop)/(b_dash_loop - a_dash_loop); % PAR Mod 17,09,2020 (moved from inside if statement above)
                        end    
                                        
                        
                    end
                    
                end
                
            end
            
        end
        
        %p.length_spike_times      = length(spike_times_vec_loop);
        length_spike_times_loop   = length(spike_times_vec_loop); % To allow parfor
        
        %RF_Ident{i} = RF_Ident_fn_v5(stimulus_arr,trig_times_vec_trunc,spike_times_vec_loop,length_spike_times_loop,p_par); % PAR Mod 27,08,2020 --> was 'RF_Ident_fn_v4' and 'trig_times_vec'
        RF_Ident(i).RF_results = RF_Ident_fn_v7(stimulus_arr,trig_times_vec_trunc,spike_times_vec_loop,length_spike_times_loop,p); % PAR Mod 17,09,2020 (RF_Ident_fn_v5 -->  RF_Ident_fn_v6) %%% PAR Mod 27,08,2020 --> was 'RF_Ident_fn_v4' and 'trig_times_vec'
        RF_Ident(i).cell_idx = cell_indices_temp(i);
        
        RF_Ident(i).add_info = settings.settings.kernel_new;
        
        
               
        
        disp(i);
        
    end
    
    delete(gcp('nocreate'));
    
else %  Parpool == 2 % Parpool off
    
    for i = 1:True_Num_Cells
        
        spike_times_vec_loop      = spiketimestamps(:,i);
        
        test_for_NaN = find(all(isnan(spike_times_vec_loop),1));
        if test_for_NaN
            RF_Ident(i).RF_results = NaN;
            RF_Ident(i).cell_idx = cell_indices_temp(i);
            if i == 1&& bb == 1 
                RF_Ident(i).add_info = add_info.settings.kernel_new;
            end
            continue
        end
        
        spike_times_vec_loop(isnan(spike_times_vec_loop)) = []; % spiketimestamps(~isnan(spiketimestamps(:,Cell_Choice)),Cell_Choice);
        spike_times_vec_loop      = spike_times_vec_loop(spike_times_vec_loop>(first_trig_time + min_start_int));
        spike_times_vec_loop      = spike_times_vec_loop(spike_times_vec_loop<(last_trig_time + min_end_int));
        
        if p.Num_trigs > p.stim_frames % Frozen noise repeated case % PAR Mod 27,08,2020 (whole if statement)
            
            if sum(inter_noise_int_vec <= max_trig_int) == p.Num_FNoise_rep_ceil-1 % If the frozen noise repeates w/o a gap
                
                % No need to remove further spikes in this case
                
                for j = 2:p.Num_FNoise_rep_ceil % map spikes in repeated stimulus chunks to original chunk
                    
                    if j==p.Num_FNoise_rep_ceil
                        frame_end_loop = p.Num_trigs - p.stim_frames*(p.Num_FNoise_rep_ceil-1);
                    else
                        frame_end_loop = p.stim_frames;
                    end
                    
                    for k = 1:frame_end_loop
                        if (k < frame_end_loop) || (j < p.Num_FNoise_rep_ceil)
                            indices_loop = find((spike_times_vec_loop>trig_times_vec((j-1)*p.stim_frames + k))&(spike_times_vec_loop<trig_times_vec((j-1)*p.stim_frames + k + 1)));
                            a_loop      = trig_times_vec(k);
                            b_loop      = trig_times_vec(k+1);
                            a_dash_loop = trig_times_vec((j-1)*p.stim_frames + k);
                            b_dash_loop = trig_times_vec((j-1)*p.stim_frames + k + 1);
                            spike_times_vec_loop(indices_loop) = a_loop + (spike_times_vec_loop(indices_loop) - a_dash_loop)*(b_loop - a_loop)/(b_dash_loop - a_dash_loop);
                        else % (k == frame_end_loop) && (j == p.Num_FNoise_rep_ceil) % No end time is given for the last frame so have to specify differently
                            indices_loop = find((spike_times_vec_loop>trig_times_vec((j-1)*p.stim_frames + k))&(spike_times_vec_loop<trig_times_vec((j-1)*p.stim_frames + k) + stim_int));
                            a_loop      = trig_times_vec(k);
                            b_loop      = trig_times_vec(k+1);
                            a_dash_loop = trig_times_vec((j-1)*p.stim_frames + k);
                            b_dash_loop = trig_times_vec((j-1)*p.stim_frames + k) + stim_int;
                            spike_times_vec_loop(indices_loop) = a_loop + (spike_times_vec_loop(indices_loop) - a_dash_loop)*(b_loop - a_loop)/(b_dash_loop - a_dash_loop);
                        end
                    end
                    
                end
                
            else % If there are gaps between (any) frozen noise chuncks
                
                for j = 1:p.Num_FNoise_rep_ceil-1 % remove spikes with stimulus windows that fall between the noise chunks
                    spike_times_vec_loop((spike_times_vec_loop>(trig_times_vec(j*p.stim_frames) + min_end_int))&(spike_times_vec_loop<trig_times_vec(j*p.stim_frames+p.Num_STE_bins+1))) = [];
                end
                
                for j = 2:p.Num_FNoise_rep_ceil % map spikes in repeated stimulus chunks to original chunk
                    
                    if j==p.Num_FNoise_rep_ceil
                        frame_end_loop = p.Num_trigs - p.stim_frames*(p.Num_FNoise_rep_ceil-1);
                    else
                        frame_end_loop = p.stim_frames;
                    end
                    
                    for k = 1:frame_end_loop
                        if k<frame_end_loop % frame_end_loop<p.stim_frames
                            indices_loop = find((spike_times_vec_loop>trig_times_vec((j-1)*p.stim_frames + k))&&(spike_times_vec_loop<trig_times_vec((j-1)*p.stim_frames + k + 1)));
                            a_loop      = trig_times_vec(k);
                            b_loop      = trig_times_vec(k+1);
                            a_dash_loop = trig_times_vec((j-1)*p.stim_frames + k);
                            b_dash_loop = trig_times_vec((j-1)*p.stim_frames + k + 1);
                            spike_times_vec_loop(indices_loop) = a_loop + (spike_times_vec_loop(indices_loop) - a_dash_loop)*(b_loop - a_loop)/(b_dash_loop - a_dash_loop);
                        else % k==frame_end_loop % frame_end_loop==p.stim_frames % No end time is given for the last frame so have to specify differently
                            indices_loop = find((spike_times_vec_loop>trig_times_vec((j-1)*p.stim_frames + k))&&(spike_times_vec_loop<trig_times_vec((j-1)*p.stim_frames + k) + stim_int));
                            a_loop      = trig_times_vec(k);
                            b_loop      = trig_times_vec(k+1);
                            a_dash_loop = trig_times_vec((j-1)*p.stim_frames + k);
                            b_dash_loop = trig_times_vec((j-1)*p.stim_frames + k) + stim_int;
                            spike_times_vec_loop(indices_loop) = a_loop + (spike_times_vec_loop(indices_loop) - a_dash_loop)*(b_loop - a_loop)/(b_dash_loop - a_dash_loop);
                        end
                    end
                    
                end
                
            end
            
        end
        
        %p.length_spike_times      = length(spike_times_vec_loop);
        length_spike_times_loop   = length(spike_times_vec_loop); % To allow parfor
        
        RF_Ident(i).RF_results = RF_Ident_fn_v7(stimulus_arr,trig_times_vec_trunc,spike_times_vec_loop,length_spike_times_loop,p); % PAR Mod 17,09,2020 (RF_Ident_fn_v5 -->  RF_Ident_fn_v6) %%% PAR Mod 27,08,2020 --> was 'RF_Ident_fn_v4' and 'trig_times_vec'
        RF_Ident(i).cell_idx = cell_indices_temp(i);
        
        %Add info if first cell
        
        RF_Ident(i).add_info = add_info.settings.kernel_new;
        
           
        
    end
    
end