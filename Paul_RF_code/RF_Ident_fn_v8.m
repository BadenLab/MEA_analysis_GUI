function RF_Ident = RF_Ident_fn_v8(stimulus_arr,trig_times_vec,spike_times_vec,length_spike_times,p)
% RF Indentification Function - ver. 8
% 18,13,2020 onwards
% Corrected problem below line 1310 (calculating RF STA etc) which arises
% if there are no significant pixels


% trig_times_vec: the vector input is 'trig_times_vec_trunc'

% Structures
% p.   --> parameters and variables from CL code
% q.   --> general parameters and variables from RF_Ident_fn
% STE. --> variables from STE_fn

% MEA data: stimulus_array := 40 rows x 40 columns x 6000 frames x 4 colours


% Find the spectral dimension (number of colours) of the stimulus
%Spectral_Dim = size(stimulus_arr,4); % (or just feed in as a structure)

% Add this way to avoid problems with paralel looping of this function
p.length_spike_times = length_spike_times;

%%% Calculate the Spike-Triggered Ensemble (STE) indices
if p.RF_Ident_Meth_vec(1) == 1 % STA-SD method
    STE = STE_fn_2(trig_times_vec,spike_times_vec,p);
end

% Preallocate and loop over spectra
if (p.STE_Full_Choice == 1) && (p.RF_Ident_Meth_vec(2) == 1 || p.RF_Ident_Meth_vec(3) == 1 || p.RF_Ident_Meth_vec(4) == 1) % LC/MI/SC method % PAR Mod 27,08,2020 --> added '(p.STE_Full_Choice == 1) &&'
    STE_Full = NaN(p.stim_rows,p.stim_columns,p.Num_STE_bins,p.length_spike_times,p.Spectral_Dim);
end

if p.RF_Ident_Meth_vec(1) == 1
    
    STA    = NaN(p.stim_rows,p.stim_columns,p.Num_STE_bins,p.Spectral_Dim);
    STA_SD = NaN(p.stim_rows,p.stim_columns,p.Spectral_Dim);
    
    % Now calculated outside function so as not to repeatedly calculate
    % within loop.
    % for each cell.
%     if p.STA_Choice  == 2 % subtract average stim
%         mean_raw_stim_arr = NaN(p.stim_rows,p.stim_columns,p.Num_STE_bins,p.Spectral_Dim);
%     end
    
end

if p.RF_Ident_Meth_vec(2) == 1
    
    LCA_Stixel_covar = NaN(p.stim_rows,p.stim_columns,p.Num_STE_bins,p.Spectral_Dim);
    LCA_Pixel_covar  = NaN(p.stim_rows,p.stim_columns,p.Spectral_Dim);
    
end

if p.RF_Ident_Meth_vec(3) == 1

    
    
end

if p.RF_Ident_Meth_vec(4) == 1
    
    meanSpikeTrigStim = NaN(p.stim_rows,p.stim_columns,p.Num_STE_bins,p.Spectral_Dim); % PAR Mod 27,08,2020
    SCA_Stixel_covar = NaN(p.stim_rows,p.stim_columns,p.Num_STE_bins,p.Spectral_Dim);
    SCA_Pixel_covar  = NaN(p.stim_rows,p.stim_columns,p.Spectral_Dim);
    
end


if p.RF_Ident_Meth_vec(1) == 1
    
    RF_Ident.STA_SD_Sig_Thresh                 = NaN(p.Spectral_Dim,1);
    
    if     p.RF_Type(1) == 1   % Box
        RF_Ident.STASD_Box_RF_coords           = cell(p.Spectral_Dim,1); % Num_RF_pixels x 2
        RF_Ident.STASD_Box_RF_coords_centre    = cell(p.Spectral_Dim,1); % 1 x 2
        RF_Ident.STASD_Box_RF_coords_noncentre = cell(p.Spectral_Dim,1); % Num_RF_pixels-1 x 2
        RF_Ident.STASD_Box_Num_RF_rows         = cell(p.Spectral_Dim,1); % 1
        RF_Ident.STASD_Box_Num_RF_cols         = cell(p.Spectral_Dim,1); % 1
        RF_Ident.STASD_Box_Num_RF_pixels       = cell(p.Spectral_Dim,1); % 1
        
        RF_Ident.STASD_Box_RF_STA     = cell(p.Spectral_Dim,1); % Num_RF_pixels x p.Num_STE_bins % PAR Mod 09,10,2020
        RF_Ident.STASD_Box_RF_STA_SD  = cell(p.Spectral_Dim,1); % Num_RF_pixels                  % PAR Mod 09,10,2020
        RF_Ident.STASD_Box_RF_STA_time  = cell(p.Spectral_Dim,1); % p.Num_STE_bins               % PAR Mod 09,10,2020
        
        if p.Spectral_Dim > 1 % PAR Mod 09,10,2020
            RF_Ident.STASD_Box_FullRF_STA     = cell(p.Spectral_Dim,1); % Num_RF_pixels x p.Num_STE_bins % PAR Mod 09,10,2020
            RF_Ident.STASD_Box_FullRF_STA_SD  = cell(p.Spectral_Dim,1); % Num_RF_pixels                  % PAR Mod 09,10,2020
            RF_Ident.STASD_Box_FullRF_STA_time  = cell(p.Spectral_Dim,1); % p.Num_STE_bins               % PAR Mod 09,10,2020
        end
    end
    if p.RF_Type(2) == 1   % All Significant Pixels
        RF_Ident.STASD_ASP_RF_coords           = cell(p.Spectral_Dim,1); % STA_SD_Num_Sig_Pix x 2
        RF_Ident.STASD_ASP_Num_RF_pixels       = cell(p.Spectral_Dim,1); % 1
        
        RF_Ident.STASD_ASP_RF_STA     = cell(p.Spectral_Dim,1); % Num_RF_pixels x p.Num_STE_bins % PAR Mod 09,10,2020
        RF_Ident.STASD_ASP_RF_STA_SD  = cell(p.Spectral_Dim,1); % Num_RF_pixels                  % PAR Mod 09,10,2020
        RF_Ident.STASD_ASP_RF_STA_time  = cell(p.Spectral_Dim,1); % p.Num_STE_bins               % PAR Mod 09,10,2020
        
        if p.Spectral_Dim > 1 % PAR Mod 09,10,2020
            RF_Ident.STASD_ASP_FullRF_STA     = cell(p.Spectral_Dim,1); % Num_RF_pixels x p.Num_STE_bins % PAR Mod 09,10,2020
            RF_Ident.STASD_ASP_FullRF_STA_SD  = cell(p.Spectral_Dim,1); % Num_RF_pixels                  % PAR Mod 09,10,2020
            RF_Ident.STASD_ASP_FullRF_STA_time  = cell(p.Spectral_Dim,1); % p.Num_STE_bins               % PAR Mod 09,10,2020
        end
    end
    if p.RF_Type(3) == 1   % Gaussian
        RF_Ident.STASD_Gaus_RF_coords          = cell(p.Spectral_Dim,1); % Num_RF_pixels x 2
        RF_Ident.STASD_Gaus_Num_RF_pixels      = cell(p.Spectral_Dim,1); % 1
        RF_Ident.STASD_Gaus_Gaussian_mean      = cell(p.Spectral_Dim,1); % 2
        RF_Ident.STASD_Gaus_Gaussian_covar     = cell(p.Spectral_Dim,1); % 2 x 2
        RF_Ident.STASD_Gaus_eval_descend_vec   = cell(p.Spectral_Dim,1); % 2 x 1
        RF_Ident.STASD_Gaus_evec_descend_mat   = cell(p.Spectral_Dim,1); % 2 x 2
        RF_Ident.STASD_Gaus_Angle_Major_Axis   = NaN(p.Spectral_Dim,1);  % 1
        RF_Ident.STASD_Gaus_Angle_Minor_Axis   = NaN(p.Spectral_Dim,1);  % 1
        RF_Ident.STASD_Gaus_Axis_eval_Ratio    = NaN(p.Spectral_Dim,1);  % 1
        RF_Ident.STASD_Gaus_Rel_Ellipse_Area   = NaN(p.Spectral_Dim,1);  % 1
        RF_Ident.STASD_Gaus_Abs_Ellipse_Area   = NaN(p.Spectral_Dim,1);  % 1
        RF_Ident.STASD_Gaus_gmPDF              = cell(p.Spectral_Dim,1); % 1 (Struct)
        RF_Ident.STASD_Gaus_Thresh_height      = cell(p.Spectral_Dim,1); % 1
        
        RF_Ident.STASD_Gaus_RF_STA    = cell(p.Spectral_Dim,1); % Num_RF_pixels x p.Num_STE_bins % PAR Mod 09,10,2020
        RF_Ident.STASD_Gaus_RF_STA_SD = cell(p.Spectral_Dim,1); % Num_RF_pixels                  % PAR Mod 09,10,2020
        RF_Ident.STASD_Gaus_RF_STA_time = cell(p.Spectral_Dim,1); % p.Num_STE_bins               % PAR Mod 09,10,2020
        
        if p.Spectral_Dim > 1 % PAR Mod 09,10,2020
            RF_Ident.STASD_Gaus_FullRF_STA    = cell(p.Spectral_Dim,1); % Num_RF_pixels x p.Num_STE_bins % PAR Mod 09,10,2020
            RF_Ident.STASD_Gaus_FullRF_STA_SD = cell(p.Spectral_Dim,1); % Num_RF_pixels                  % PAR Mod 09,10,2020
            RF_Ident.STASD_Gaus_FullRF_STA_time = cell(p.Spectral_Dim,1); % p.Num_STE_bins               % PAR Mod 09,10,2020
        end
    end
    
end

if p.RF_Ident_Meth_vec(2) == 1
    
    RF_Ident.LC_Sig_Thresh                  = NaN(p.Spectral_Dim,1);
    
    if     p.RF_Type(1) == 1   % Box
        RF_Ident.LC_Box_RF_coords           = cell(p.Spectral_Dim,1); % Num_RF_pixels x 2
        RF_Ident.LC_Box_RF_coords_centre    = cell(p.Spectral_Dim,1); % 1 x 2
        RF_Ident.LC_Box_RF_coords_noncentre = cell(p.Spectral_Dim,1); % Num_RF_pixels-1 x 2
        RF_Ident.LC_Box_Num_RF_rows         = cell(p.Spectral_Dim,1); % 1
        RF_Ident.LC_Box_Num_RF_cols         = cell(p.Spectral_Dim,1); % 1
        RF_Ident.LC_Box_Num_RF_pixels       = cell(p.Spectral_Dim,1); % 1
    end
    if p.RF_Type(2) == 1   % All Significant Pixels
        RF_Ident.LC_ASP_RF_coords           = cell(p.Spectral_Dim,1); % STA_SD_Num_Sig_Pix x 2
        RF_Ident.LC_ASP_Num_RF_pixels       = cell(p.Spectral_Dim,1); % 1
    end
    if p.RF_Type(3) == 1   % Gaussian
        RF_Ident.LC_Gaus_RF_coords          = cell(p.Spectral_Dim,1); % Num_RF_pixels x 2
        RF_Ident.LC_Gaus_Num_RF_pixels      = cell(p.Spectral_Dim,1); % 1
        RF_Ident.LC_Gaus_Gaussian_mean      = cell(p.Spectral_Dim,1); % 2
        RF_Ident.LC_Gaus_Gaussian_covar     = cell(p.Spectral_Dim,1); % 2 x 2
        RF_Ident.LC_Gaus_eval_descend_vec   = cell(p.Spectral_Dim,1); % 2 x 1
        RF_Ident.LC_Gaus_evec_descend_mat   = cell(p.Spectral_Dim,1); % 2 x 2
        RF_Ident.LC_Gaus_Angle_Major_Axis   = NaN(p.Spectral_Dim,1);  % 1
        RF_Ident.LC_Gaus_Angle_Minor_Axis   = NaN(p.Spectral_Dim,1);  % 1
        RF_Ident.LC_Gaus_Axis_eval_Ratio    = NaN(p.Spectral_Dim,1);  % 1
        RF_Ident.LC_Gaus_Rel_Ellipse_Area   = NaN(p.Spectral_Dim,1);  % 1
        RF_Ident.LC_Gaus_Abs_Ellipse_Area   = NaN(p.Spectral_Dim,1);  % 1
        RF_Ident.LC_Gaus_gmPDF              = cell(p.Spectral_Dim,1); % 1 (Struct)
        RF_Ident.LC_Gaus_Thresh_height      = cell(p.Spectral_Dim,1); % 1
    end
    
end

if p.RF_Ident_Meth_vec(3) == 1
    
    RF_Ident.MI_Sig_Thresh                  = NaN(p.Spectral_Dim,1);
    
    if     p.RF_Type(1) == 1   % Box
        RF_Ident.MI_Box_RF_coords           = cell(p.Spectral_Dim,1); % Num_RF_pixels x 2
        RF_Ident.MI_Box_RF_coords_centre    = cell(p.Spectral_Dim,1); % 1 x 2
        RF_Ident.MI_Box_RF_coords_noncentre = cell(p.Spectral_Dim,1); % Num_RF_pixels-1 x 2
        RF_Ident.MI_Box_Num_RF_rows         = cell(p.Spectral_Dim,1); % 1
        RF_Ident.MI_Box_Num_RF_cols         = cell(p.Spectral_Dim,1); % 1
        RF_Ident.MI_Box_Num_RF_pixels       = cell(p.Spectral_Dim,1); % 1
    end
    if p.RF_Type(2) == 1   % All Significant Pixels
        RF_Ident.MI_ASP_RF_coords           = cell(p.Spectral_Dim,1); % STA_SD_Num_Sig_Pix x 2
        RF_Ident.MI_ASP_Num_RF_pixels       = cell(p.Spectral_Dim,1); % 1
    end
    if p.RF_Type(3) == 1   % Gaussian
        RF_Ident.MI_Gaus_RF_coords          = cell(p.Spectral_Dim,1); % Num_RF_pixels x 2
        RF_Ident.MI_Gaus_Num_RF_pixels      = cell(p.Spectral_Dim,1); % 1
        RF_Ident.MI_Gaus_Gaussian_mean      = cell(p.Spectral_Dim,1); % 2
        RF_Ident.MI_Gaus_Gaussian_covar     = cell(p.Spectral_Dim,1); % 2 x 2
        RF_Ident.MI_Gaus_eval_descend_vec   = cell(p.Spectral_Dim,1); % 2 x 1
        RF_Ident.MI_Gaus_evec_descend_mat   = cell(p.Spectral_Dim,1); % 2 x 2
        RF_Ident.MI_Gaus_Angle_Major_Axis   = NaN(p.Spectral_Dim,1);  % 1
        RF_Ident.MI_Gaus_Angle_Minor_Axis   = NaN(p.Spectral_Dim,1);  % 1
        RF_Ident.MI_Gaus_Axis_eval_Ratio    = NaN(p.Spectral_Dim,1);  % 1
        RF_Ident.MI_Gaus_Rel_Ellipse_Area   = NaN(p.Spectral_Dim,1);  % 1
        RF_Ident.MI_Gaus_Abs_Ellipse_Area   = NaN(p.Spectral_Dim,1);  % 1
        RF_Ident.MI_Gaus_gmPDF              = cell(p.Spectral_Dim,1); % 1 (Struct)
        RF_Ident.MI_Gaus_Thresh_height      = cell(p.Spectral_Dim,1); % 1
    end
    
end

if p.RF_Ident_Meth_vec(4) == 1
    
    RF_Ident.SC_Sig_Thresh                  = NaN(p.Spectral_Dim,1);
    
    if     p.RF_Type(1) == 1   % Box
        RF_Ident.SC_Box_RF_coords           = cell(p.Spectral_Dim,1); % Num_RF_pixels x 2
        RF_Ident.SC_Box_RF_coords_centre    = cell(p.Spectral_Dim,1); % 1 x 2
        RF_Ident.SC_Box_RF_coords_noncentre = cell(p.Spectral_Dim,1); % Num_RF_pixels-1 x 2
        RF_Ident.SC_Box_Num_RF_rows         = cell(p.Spectral_Dim,1); % 1
        RF_Ident.SC_Box_Num_RF_cols         = cell(p.Spectral_Dim,1); % 1
        RF_Ident.SC_Box_Num_RF_pixels       = cell(p.Spectral_Dim,1); % 1
        
        RF_Ident.SC_Box_RF_Stixel_covar  = cell(p.Spectral_Dim,1); % Num_RF_pixels x p.Num_STE_bins % PAR Mod 09,10,2020
        RF_Ident.SC_Box_RF_Pixel_covar   = cell(p.Spectral_Dim,1); % Num_RF_pixels                  % PAR Mod 09,10,2020
        RF_Ident.SC_Box_RF_Stixel_covar_time  = cell(p.Spectral_Dim,1); % p.Num_STE_bins            % PAR Mod 09,10,2020
        
        if p.Spectral_Dim > 1 % PAR Mod 09,10,2020
            RF_Ident.SC_Box_FullRF_Stixel_covar  = cell(p.Spectral_Dim,1); % Num_RF_pixels x p.Num_STE_bins % PAR Mod 09,10,2020
            RF_Ident.SC_Box_FullRF_Pixel_covar   = cell(p.Spectral_Dim,1); % Num_RF_pixels                  % PAR Mod 09,10,2020
            RF_Ident.SC_Box_FullRF_Stixel_covar_time  = cell(p.Spectral_Dim,1); % p.Num_STE_bins            % PAR Mod 09,10,2020
        end
    end
    if p.RF_Type(2) == 1   % All Significant Pixels
        RF_Ident.SC_ASP_RF_coords           = cell(p.Spectral_Dim,1); % STA_SD_Num_Sig_Pix x 2
        RF_Ident.SC_ASP_Num_RF_pixels       = cell(p.Spectral_Dim,1); % 1
        
        RF_Ident.SC_ASP_RF_Stixel_covar  = cell(p.Spectral_Dim,1); % Num_RF_pixels x p.Num_STE_bins % PAR Mod 09,10,2020
        RF_Ident.SC_ASP_RF_Pixel_covar   = cell(p.Spectral_Dim,1); % Num_RF_pixels                  % PAR Mod 09,10,2020
        RF_Ident.SC_ASP_RF_Stixel_covar_time  = cell(p.Spectral_Dim,1); % p.Num_STE_bins            % PAR Mod 09,10,2020
        
        if p.Spectral_Dim > 1 % PAR Mod 09,10,2020
            RF_Ident.SC_ASP_FullRF_Stixel_covar  = cell(p.Spectral_Dim,1); % Num_RF_pixels x p.Num_STE_bins % PAR Mod 09,10,2020
            RF_Ident.SC_ASP_FullRF_Pixel_covar   = cell(p.Spectral_Dim,1); % Num_RF_pixels                  % PAR Mod 09,10,2020
            RF_Ident.SC_ASP_FullRF_Stixel_covar_time  = cell(p.Spectral_Dim,1); % p.Num_STE_bins            % PAR Mod 09,10,2020
        end
    end
    if p.RF_Type(3) == 1   % Gaussian
        RF_Ident.SC_Gaus_RF_coords          = cell(p.Spectral_Dim,1); % Num_RF_pixels x 2
        RF_Ident.SC_Gaus_Num_RF_pixels      = cell(p.Spectral_Dim,1); % 1
        RF_Ident.SC_Gaus_Gaussian_mean      = cell(p.Spectral_Dim,1); % 2
        RF_Ident.SC_Gaus_Gaussian_covar     = cell(p.Spectral_Dim,1); % 2 x 2
        RF_Ident.SC_Gaus_eval_descend_vec   = cell(p.Spectral_Dim,1); % 2 x 1
        RF_Ident.SC_Gaus_evec_descend_mat   = cell(p.Spectral_Dim,1); % 2 x 2
        RF_Ident.SC_Gaus_Angle_Major_Axis   = NaN(p.Spectral_Dim,1);  % 1
        RF_Ident.SC_Gaus_Angle_Minor_Axis   = NaN(p.Spectral_Dim,1);  % 1
        RF_Ident.SC_Gaus_Axis_eval_Ratio    = NaN(p.Spectral_Dim,1);  % 1
        RF_Ident.SC_Gaus_Rel_Ellipse_Area   = NaN(p.Spectral_Dim,1);  % 1
        RF_Ident.SC_Gaus_Abs_Ellipse_Area   = NaN(p.Spectral_Dim,1);  % 1
        RF_Ident.SC_Gaus_gmPDF              = cell(p.Spectral_Dim,1); % 1 (Struct)
        RF_Ident.SC_Gaus_Thresh_height      = cell(p.Spectral_Dim,1); % 1
        
        RF_Ident.SC_Gaus_RF_Stixel_covar = cell(p.Spectral_Dim,1); % Num_RF_pixels x p.Num_STE_bins % PAR Mod 09,10,2020
        RF_Ident.SC_Gaus_RF_Pixel_covar  = cell(p.Spectral_Dim,1); % Num_RF_pixels                  % PAR Mod 09,10,2020
        RF_Ident.SC_Gaus_RF_Stixel_covar_time = cell(p.Spectral_Dim,1); % p.Num_STE_bins            % PAR Mod 09,10,2020
        
        if p.Spectral_Dim > 1 % PAR Mod 09,10,2020
            RF_Ident.SC_Gaus_FullRF_Stixel_covar = cell(p.Spectral_Dim,1); % Num_RF_pixels x p.Num_STE_bins % PAR Mod 09,10,2020
            RF_Ident.SC_Gaus_FullRF_Pixel_covar  = cell(p.Spectral_Dim,1); % Num_RF_pixels                  % PAR Mod 09,10,2020
            RF_Ident.SC_Gaus_FullRF_Stixel_covar_time = cell(p.Spectral_Dim,1); % p.Num_STE_bins            % PAR Mod 09,10,2020
        end
    end
    
end

for i = 1:p.Spectral_Dim
    
    %%% Calculate the Full Spike-Triggered Ensemble (STE) array
    % p.stim_rows  x  p.stim_columns  x  p.Num_STE_bins  x  p.length_spike_times
    if p.STE_Full_Choice == 1 % PAR Mod 27,08,2020 --> added if statement
        STE_Full(:,:,:,:,i) = STE_Full_fn_2(stimulus_arr(:,:,:,i),trig_times_vec,spike_times_vec,p);
    end
    
    %%% Calculate the Spike-Triggered Average (STA)
    if p.RF_Ident_Meth_vec(1) == 1
        
        if p.STA_Choice  == 2 % subtract average stim
            % Calculate the mean of the raw stimuli
            % Performed outside function to avoid repeated calc
            % mean_raw_stim_arr(:,:,:,i) = mean_raw_stim_SpaceTime_fn(stimulus_arr(:,:,:,i),trig_times_vec,p);
            if p.STE_Full_Choice == 1 % STE_Full % PAR Mod 27,08,2020 (whole if statement)
                STA(:,:,:,i)               = STA_fn_v2(STE_Full(:,:,:,:,i),p.mean_raw_stim_arr(:,:,:,i),p);
            else % p.STE_Full_Choice == 1 % No STE_Full
                STA(:,:,:,i)               = STA_fn_2(stimulus_arr(:,:,:,i),p,STE,p.mean_raw_stim_arr(:,:,:,i));
            end
            
            if (p.STE_Full_Choice == 2) && (p.RF_Ident_Meth_vec(4) == 1) % No STE_Full and % PAR Mod 27,08,2020 (whole if statement)
               meanSpikeTrigStim(:,:,:,i) = meanSTS_fn_2(stimulus_arr(:,:,:,i),p,STE); 
            end
        else
            STA(:,:,:,i)               = STA_fn_2(stimulus_arr(:,:,:,i),p,STE,[]);
        end
        
        %%% Calculate the STA-SD
        STA_SD(:,:,i) = std(STA(:,:,:,i),[],3);
        
    end
    
    %%% Calculate the Local Covariance Array   (LCA)
    if p.RF_Ident_Meth_vec(2) == 1
        % function LCA = LCA_fn(stimulus_arr,p,STE)
        [LCA_Stixel_covar(:,:,:,i),LCA_Pixel_covar(:,:,i)] = LC_fn_v2(STE_Full(:,:,:,:,i),p);
    end
    
    %%% Calculate the Mutual Information Array (MIA)
    if p.RF_Ident_Meth_vec(3) == 1
        % function MIA = MIA_fn(stimulus_arr,p,STE)
    end
    
    %%% Calculate the Self Covariance Array (SCA)
    if p.RF_Ident_Meth_vec(4) == 1
        if p.STE_Full_Choice == 1 % PAR Mod 27,08,2020 --> added if statement
            [SCA_Stixel_covar(:,:,:,i),SCA_Pixel_covar(:,:,i)] = SC_fn(STE_Full(:,:,:,:,i),p);
            %[SCA_Stixel_covar(:,:,:,i),SCA_Pixel_covar(:,:,i)] = SC_fn_v3(STE_Full(:,:,:,:,i),STA(:,:,:,i),p); % PAR Mod % 27,08,2020 (Not used, keep commented)
        else % p.STE_Full_Choice == 2
            %[SCA_Stixel_covar(:,:,:,i),SCA_Pixel_covar(:,:,i)] = SC_fn_v2(stimulus_arr(:,:,:,i),STE,STA(:,:,:,i),p); % PAR Mod 27,08,2020 --> changed from 'SC_fn(STE_Full(:,:,:,:,i),p)' to 'SC_fn_v2()'
            [SCA_Stixel_covar(:,:,:,i),SCA_Pixel_covar(:,:,i)] = SC_fn_v5(stimulus_arr(:,:,:,i),STE,meanSpikeTrigStim(:,:,:,i),p); % PAR Mod 27,08,2020
        end
    end
    
    %%% Find the Indices and Coordinates of Pixels Meeting the QC
    
    % STA SD RF
    if p.RF_Ident_Meth_vec(1) == 1
        
        if p.STA_SD_QC_Type == 1
            
            Mean_STA_SD                    = mean(STA_SD(:,:,i),'all');
            SD_STA_SD                      = std(STA_SD(:,:,i),[],'all');
            RF_Ident.STA_SD_Sig_Thresh(i)  = Mean_STA_SD + p.STA_SD_Thresh*SD_STA_SD;
            STA_SD_Num_Sig_Pix             = sum(STA_SD(:,:,i)>RF_Ident.STA_SD_Sig_Thresh(i),'all');
            
        else % if p.STA_SD_QC_Type == 2
            
            STA_SD_prob_dist               = fitdist(reshape(STA_SD(:,:,i),[p.stim_pixels,1]),'Normal');
            STA_SD_Upr_bd                  = 1 - 0.01*p.STA_SD_CI_Per;
            RF_Ident.STA_SD_Sig_Thresh(i)  = icdf(STA_SD_prob_dist,STA_SD_Upr_bd);
            STA_SD_Num_Sig_Pix             = sum(STA_SD(:,:,i)>RF_Ident.STA_SD_Sig_Thresh(i),'all');
            
        end
        
        if STA_SD_Num_Sig_Pix > 0
            
            if     p.RF_Type(1) == 1   % Box
                
                [RF_Ident.STASD_Box_RF_coords{i},...
                    RF_Ident.STASD_Box_RF_coords_centre{i},...
                    RF_Ident.STASD_Box_RF_coords_noncentre{i},...
                    RF_Ident.STASD_Box_Num_RF_rows{i},...
                    RF_Ident.STASD_Box_Num_RF_cols{i},...
                    RF_Ident.STASD_Box_Num_RF_pixels{i}] = Box_RF_fn_v2(STA_SD(:,:,i),p,1);
                
            end
            
            if p.RF_Type(2) == 1   % All Significant Pixels
                
                [RF_Ident.STASD_ASP_RF_coords{i},RF_Ident.STASD_ASP_Num_RF_pixels{i}] = AllSigPix_RF_fn_v3(STA_SD(:,:,i),RF_Ident.STA_SD_Sig_Thresh(i),1,STA_SD_Num_Sig_Pix);
                
            end
            
            if p.RF_Type(3) == 1 && p.RF_Type(2) == 1   % Gaussian
                
                [RF_Ident.STASD_Gaus_RF_coords{i},...
                    RF_Ident.STASD_Gaus_Num_RF_pixels{i},...
                    RF_Ident.STASD_Gaus_Gaussian_mean{i},...
                    RF_Ident.STASD_Gaus_Gaussian_covar{i},...
                    RF_Ident.STASD_Gaus_eval_descend_vec{i},...
                    RF_Ident.STASD_Gaus_evec_descend_mat{i},...
                    RF_Ident.STASD_Gaus_Angle_Major_Axis(i),...
                    RF_Ident.STASD_Gaus_Angle_Minor_Axis(i),...
                    RF_Ident.STASD_Gaus_Axis_eval_Ratio(i),...
                    RF_Ident.STASD_Gaus_Rel_Ellipse_Area(i),...
                    RF_Ident.STASD_Gaus_Abs_Ellipse_Area(i),...
                    RF_Ident.STASD_Gaus_gmPDF{i},...
                    RF_Ident.STASD_Gaus_Thresh_height{i}] = Gaussian_RF_fn_v3(RF_Ident.STASD_ASP_RF_coords{i},RF_Ident.STASD_ASP_Num_RF_pixels{i},p);
                
            elseif p.RF_Type(3) == 1 && p.RF_Type(2) == 0
                
                [STASD_ASP_RF_coords_temp,STASD_ASP_Num_RF_pixels_temp] = AllSigPix_RF_fn_v3(STA_SD(:,:,i),RF_Ident.STA_SD_Sig_Thresh(i),1,STA_SD_Num_Sig_Pix);
                
                [RF_Ident.STASD_Gaus_RF_coords{i},...
                    RF_Ident.STASD_Gaus_Num_RF_pixels{i},...
                    RF_Ident.STASD_Gaus_Gaussian_mean{i},...
                    RF_Ident.STASD_Gaus_Gaussian_covar{i},...
                    RF_Ident.STASD_Gaus_eval_descend_vec{i},...
                    RF_Ident.STASD_Gaus_evec_descend_mat{i},...
                    RF_Ident.STASD_Gaus_Angle_Major_Axis(i),...
                    RF_Ident.STASD_Gaus_Angle_Minor_Axis(i),...
                    RF_Ident.STASD_Gaus_Axis_eval_Ratio(i),...
                    RF_Ident.STASD_Gaus_Rel_Ellipse_Area(i),...
                    RF_Ident.STASD_Gaus_Abs_Ellipse_Area(i),...
                    RF_Ident.STASD_Gaus_gmPDF{i},...
                    RF_Ident.STASD_Gaus_Thresh_height{i}] = Gaussian_RF_fn_v3(STASD_ASP_RF_coords_temp,STASD_ASP_Num_RF_pixels_temp,p);
                
            end
            
        else
            
            %disp(sprintf('There are no significant pixels as judged by the STA-SD criterion for colour %i',i));
            
        end
        
    end
    
    % LCA RF
    if p.RF_Ident_Meth_vec(2) == 1
        
        if p.LC_QC_Type == 1
            
            Mean_LC                    = mean(LCA_Pixel_covar(:,:,i),'all');
            SD_LC                      = std(LCA_Pixel_covar(:,:,i),[],'all');
            RF_Ident.LC_Sig_Thresh(i)  = Mean_LC + p.LC_Thresh*SD_LC;
            LC_Num_Sig_Pix             = sum(LCA_Pixel_covar(:,:,i)>RF_Ident.LC_Sig_Thresh(i),'all');
            
        else % if p.LC_QC_Type == 2
            
            LC_prob_dist               = fitdist(reshape(LCA_Pixel_covar(:,:,i),[p.stim_pixels,1]),'Normal');
            LC_Upr_bd                  = 1 - 0.01*p.LC_CI_Per;
            RF_Ident.LC_Sig_Thresh(i)  = icdf(LC_prob_dist,LC_Upr_bd);
            LC_Num_Sig_Pix             = sum(LCA_Pixel_covar(:,:,i)>RF_Ident.LC_Sig_Thresh(i),'all');
            
        end
        
        if LC_Num_Sig_Pix > 0
            
            if     p.RF_Type(1) == 1   % Box
                
                [RF_Ident.LC_Box_RF_coords{i},...
                    RF_Ident.LC_Box_RF_coords_centre{i},...
                    RF_Ident.LC_Box_RF_coords_noncentre{i},...
                    RF_Ident.LC_Box_Num_RF_rows{i},...
                    RF_Ident.LC_Box_Num_RF_cols{i},...
                    RF_Ident.LC_Box_Num_RF_pixels{i}] = Box_RF_fn_v2(LCA_Pixel_covar(:,:,i),p,2);
                
            end
            
            if p.RF_Type(2) == 1   % All Significant Pixels
                
                [RF_Ident.LC_ASP_RF_coords{i},RF_Ident.LC_ASP_Num_RF_pixels{i}] = AllSigPix_RF_fn_v3(LCA_Pixel_covar(:,:,i),RF_Ident.LC_Sig_Thresh(i),2,LC_Num_Sig_Pix);
                
            end
            
            if p.RF_Type(3) == 1 && p.RF_Type(2) == 1  % Gaussian
                
                [RF_Ident.LC_Gaus_RF_coords{i},...
                    RF_Ident.LC_Gaus_Num_RF_pixels{i},...
                    RF_Ident.LC_Gaus_Gaussian_mean{i},...
                    RF_Ident.LC_Gaus_Gaussian_covar{i},...
                    RF_Ident.LC_Gaus_eval_descend_vec{i},...
                    RF_Ident.LC_Gaus_evec_descend_mat{i},...
                    RF_Ident.LC_Gaus_Angle_Major_Axis(i),...
                    RF_Ident.LC_Gaus_Angle_Minor_Axis(i),...
                    RF_Ident.LC_Gaus_Axis_eval_Ratio(i),...
                    RF_Ident.LC_Gaus_Rel_Ellipse_Area(i),...
                    RF_Ident.LC_Gaus_Abs_Ellipse_Area(i),...
                    RF_Ident.LC_Gaus_gmPDF{i},...
                    RF_Ident.LC_Gaus_Thresh_height{i}] = Gaussian_RF_fn_v3(RF_Ident.LC_ASP_RF_coords{i},RF_Ident.LC_ASP_Num_RF_pixels{i},p);
                
            elseif p.RF_Type(3) == 1 && p.RF_Type(2) == 0
                
                [LC_ASP_RF_coords_temp,LC_ASP_Num_RF_pixels_temp] = AllSigPix_RF_fn_v3(LCA_Pixel_covar(:,:,i),RF_Ident.LC_Sig_Thresh(i),2,LC_Num_Sig_Pix);
                
                [RF_Ident.LC_Gaus_RF_coords{i},...
                    RF_Ident.LC_Gaus_Num_RF_pixels{i},...
                    RF_Ident.LC_Gaus_Gaussian_mean{i},...
                    RF_Ident.LC_Gaus_Gaussian_covar{i},...
                    RF_Ident.LC_Gaus_eval_descend_vec{i},...
                    RF_Ident.LC_Gaus_evec_descend_mat{i},...
                    RF_Ident.LC_Gaus_Angle_Major_Axis(i),...
                    RF_Ident.LC_Gaus_Angle_Minor_Axis(i),...
                    RF_Ident.LC_Gaus_Axis_eval_Ratio(i),...
                    RF_Ident.LC_Gaus_Rel_Ellipse_Area(i),...
                    RF_Ident.LC_Gaus_Abs_Ellipse_Area(i),...
                    RF_Ident.LC_Gaus_gmPDF{i},...
                    RF_Ident.LC_Gaus_Thresh_height{i}] = Gaussian_RF_fn_v3(LC_ASP_RF_coords_temp,LC_ASP_Num_RF_pixels_temp,p);
                
            end
            
        else
            
            %disp(sprintf('There are no significant pixels as judged by the LC criterion for colour %i',i));
            %LC_RF = [];
            
        end
        
    end
    
    % MIA RF
    if p.RF_Ident_Meth_vec(3) == 1
        
        if p.MI_QC_Type == 1
            
            Mean_MI                    = mean(MIA(:,:,i),'all');
            SD_MI                      = std(MIA(:,:,i),[],'all');
            RF_Ident.MI_Sig_Thresh(i)  = Mean_MI + p.MI_Thresh*SD_MI;
            MI_Num_Sig_Pix             = sum(MIA(:,:,i)>RF_Ident.MI_Sig_Thresh(i),'all');
            
        else % if p.MI_QC_Type == 2
            
            MI_prob_dist               = fitdist(reshape(MIA(:,:,i),[p.stim_pixels,1]),'Normal');
            MI_Upr_bd                  = 1 - 0.01*p.MI_CI_Per;
            RF_Ident.MI_Sig_Thresh(i)  = icdf(MI_prob_dist,MI_Upr_bd);
            MI_Num_Sig_Pix             = sum(MIA(:,:,i)>RF_Ident.MI_Sig_Thresh(i),'all');
            
        end
        
        if MI_Num_Sig_Pix > 0
            
            if     p.RF_Type(1) == 1   % Box
                
                [RF_Ident.MI_Box_RF_coords{i},...
                    RF_Ident.MI_Box_RF_coords_centre{i},...
                    RF_Ident.MI_Box_RF_coords_noncentre{i},...
                    RF_Ident.MI_Box_Num_RF_rows{i},...
                    RF_Ident.MI_Box_Num_RF_cols{i},...
                    RF_Ident.MI_Box_Num_RF_pixels{i}] = Box_RF_fn_v2(MIA(:,:,i),p,3);
                
            end
            
            if p.RF_Type(2) == 1   % All Significant Pixels
                
                [RF_Ident.MI_ASP_RF_coords{i},RF_Ident.MI_ASP_Num_RF_pixels{i}] = AllSigPix_RF_fn_v3(MIA(:,:,i),RF_Ident.MI_Sig_Thresh(i),3,MI_Num_Sig_Pix);
                
            end
            
            
            if p.RF_Type(3) == 1 && p.RF_Type(2) == 1  % Gaussian
                
                [RF_Ident.MI_Gaus_RF_coords{i},...
                    RF_Ident.MI_Gaus_Num_RF_pixels{i},...
                    RF_Ident.MI_Gaus_Gaussian_mean{i},...
                    RF_Ident.MI_Gaus_Gaussian_covar{i},...
                    RF_Ident.MI_Gaus_eval_descend_vec{i},...
                    RF_Ident.MI_Gaus_evec_descend_mat{i},...
                    RF_Ident.MI_Gaus_Angle_Major_Axis(i),...
                    RF_Ident.MI_Gaus_Angle_Minor_Axis(i),...
                    RF_Ident.MI_Gaus_Axis_eval_Ratio(i),...
                    RF_Ident.MI_Gaus_Rel_Ellipse_Area(i),...
                    RF_Ident.MI_Gaus_Abs_Ellipse_Area(i),...
                    RF_Ident.MI_Gaus_gmPDF{i},...
                    RF_Ident.MI_Gaus_Thresh_height{i}] = Gaussian_RF_fn_v3(RF_Ident.MI_ASP_RF_coords{i},RF_Ident.MI_ASP_Num_RF_pixels{i},p);
                
            elseif p.RF_Type(3) == 1 && p.RF_Type(2) == 0
                
                [MI_ASP_RF_coords_temp,MI_ASP_Num_RF_pixels_temp] = AllSigPix_RF_fn_v3(MIA(:,:,i),RF_Ident.MI_Sig_Thresh(i),3,MI_Num_Sig_Pix);
                
                 [RF_Ident.MI_Gaus_RF_coords{i},...
                    RF_Ident.MI_Gaus_Num_RF_pixels{i},...
                    RF_Ident.MI_Gaus_Gaussian_mean{i},...
                    RF_Ident.MI_Gaus_Gaussian_covar{i},...
                    RF_Ident.MI_Gaus_eval_descend_vec{i},...
                    RF_Ident.MI_Gaus_evec_descend_mat{i},...
                    RF_Ident.MI_Gaus_Angle_Major_Axis(i),...
                    RF_Ident.MI_Gaus_Angle_Minor_Axis(i),...
                    RF_Ident.MI_Gaus_Axis_eval_Ratio(i),...
                    RF_Ident.MI_Gaus_Rel_Ellipse_Area(i),...
                    RF_Ident.MI_Gaus_Abs_Ellipse_Area(i),...
                    RF_Ident.MI_Gaus_gmPDF{i},...
                    RF_Ident.MI_Gaus_Thresh_height{i}] = Gaussian_RF_fn_v3(MI_ASP_RF_coords_temp,MI_ASP_Num_RF_pixels_temp,p);
                
            end
            
        else
            
            %disp(sprintf('There are no significant pixels as judged by the MI criterion for colour %i',i));
            
        end
        
    end
    
    % SCA RF
    if p.RF_Ident_Meth_vec(4) == 1
        
        if p.SC_QC_Type == 1
            
            Mean_SC                    = mean(SCA_Pixel_covar(:,:,i),'all');
            SD_SC                      = std(SCA_Pixel_covar(:,:,i),[],'all');
            RF_Ident.SC_Sig_Thresh(i)  = Mean_SC - p.SC_Thresh*SD_SC; % minus here as small var is better
            SC_Num_Sig_Pix             = sum(SCA_Pixel_covar(:,:,i)<RF_Ident.SC_Sig_Thresh(i),'all'); % < here as small var is better
            
        else % if p.SC_QC_Type == 2
            
            SC_prob_dist               = fitdist(reshape(SCA_Pixel_covar(:,:,i),[p.stim_pixels,1]),'Normal');
            SC_Lwr_bd                  = 0.01*p.SC_CI_Per; % rather than 1 - 0.01*... as small var is better
            RF_Ident.SC_Sig_Thresh(i)  = icdf(SC_prob_dist,SC_Lwr_bd);
            SC_Num_Sig_Pix             = sum(SCA_Pixel_covar(:,:,i)<RF_Ident.SC_Sig_Thresh(i),'all'); % < here as small var is better
            
        end
        
        if SC_Num_Sig_Pix > 0
            
            if     p.RF_Type(1) == 1   % Box
                
                [RF_Ident.SC_Box_RF_coords{i},...
                    RF_Ident.SC_Box_RF_coords_centre{i},...
                    RF_Ident.SC_Box_RF_coords_noncentre{i},...
                    RF_Ident.SC_Box_Num_RF_rows{i},...
                    RF_Ident.SC_Box_Num_RF_cols{i},...
                    RF_Ident.SC_Box_Num_RF_pixels{i}] = Box_RF_fn_v2(SCA_Pixel_covar(:,:,i),p,4);
                
            end
            
            if p.RF_Type(2) == 1   % All Significant Pixels
                
                [RF_Ident.SC_ASP_RF_coords{i},RF_Ident.SC_ASP_Num_RF_pixels{i}] = AllSigPix_RF_fn_v3(SCA_Pixel_covar(:,:,i),RF_Ident.SC_Sig_Thresh(i),4,SC_Num_Sig_Pix);
                
            end
            
            if p.RF_Type(3) == 1 && p.RF_Type(2) == 1   % Gaussian
                
                [RF_Ident.SC_Gaus_RF_coords{i},...
                    RF_Ident.SC_Gaus_Num_RF_pixels{i},...
                    RF_Ident.SC_Gaus_Gaussian_mean{i},...
                    RF_Ident.SC_Gaus_Gaussian_covar{i},...
                    RF_Ident.SC_Gaus_eval_descend_vec{i},...
                    RF_Ident.SC_Gaus_evec_descend_mat{i},...
                    RF_Ident.SC_Gaus_Angle_Major_Axis(i),...
                    RF_Ident.SC_Gaus_Angle_Minor_Axis(i),...
                    RF_Ident.SC_Gaus_Axis_eval_Ratio(i),...
                    RF_Ident.SC_Gaus_Rel_Ellipse_Area(i),...
                    RF_Ident.SC_Gaus_Abs_Ellipse_Area(i),...
                    RF_Ident.SC_Gaus_gmPDF{i},...
                    RF_Ident.SC_Gaus_Thresh_height{i}] = Gaussian_RF_fn_v3(RF_Ident.SC_ASP_RF_coords{i},RF_Ident.SC_ASP_Num_RF_pixels{i},p);
                
            elseif p.RF_Type(3) == 1 && p.RF_Type(2) == 0
                
                [SC_ASP_RF_coords_temp,SC_ASP_Num_RF_pixels_temp] = AllSigPix_RF_fn_v3(SCA_Pixel_covar(:,:,i),RF_Ident.SC_Sig_Thresh(i),4,SC_Num_Sig_Pix);
                
                [RF_Ident.SC_Gaus_RF_coords{i},...
                    RF_Ident.SC_Gaus_Num_RF_pixels{i},...
                    RF_Ident.SC_Gaus_Gaussian_mean{i},...
                    RF_Ident.SC_Gaus_Gaussian_covar{i},...
                    RF_Ident.SC_Gaus_eval_descend_vec{i},...
                    RF_Ident.SC_Gaus_evec_descend_mat{i},...
                    RF_Ident.SC_Gaus_Angle_Major_Axis(i),...
                    RF_Ident.SC_Gaus_Angle_Minor_Axis(i),...
                    RF_Ident.SC_Gaus_Axis_eval_Ratio(i),...
                    RF_Ident.SC_Gaus_Rel_Ellipse_Area(i),...
                    RF_Ident.SC_Gaus_Abs_Ellipse_Area(i),...
                    RF_Ident.SC_Gaus_gmPDF{i},...
                    RF_Ident.SC_Gaus_Thresh_height{i}] = Gaussian_RF_fn_v3(SC_ASP_RF_coords_temp,SC_ASP_Num_RF_pixels_temp,p);
                
            end
            
        else
            
            %disp(sprintf('There are no significant pixels as judged by the SC criterion for colour %i',i));
            
        end
        
    end
    
    
end


%% Include additional variables in RF_Indent

if     p.RF_Ident_Meth_vec(1) == 1
    
    RF_Ident.STA_SD = STA_SD;
    
end

if p.RF_Ident_Meth_vec(2) == 1
    
    RF_Ident.LCA_Stixel_covar = LCA_Stixel_covar;
    RF_Ident.LCA_Pixel_covar  = LCA_Pixel_covar;
    
end

if p.RF_Ident_Meth_vec(3) == 1
    
    
    
end

if p.RF_Ident_Meth_vec(4) == 1
    
    RF_Ident.SCA_Stixel_covar = SCA_Stixel_covar;
    RF_Ident.SCA_Pixel_covar  = SCA_Pixel_covar;
    
end


%% Find the Full RF (over all colours)

if p.Spectral_Dim > 1
    
    if     p.RF_Ident_Meth_vec(1) == 1 % STA-SD
        
        % Check if there are any significant pixels across all of the colours
        Sig_Pix_Indicator_STASD = 0;
        Colour_with_Sig_Pix_STASD_vec = NaN(p.Spectral_Dim,1);
        for i = 1:p.Spectral_Dim
            if ~isempty(RF_Ident.STASD_ASP_Num_RF_pixels{i})
                Colour_with_Sig_Pix_STASD_vec(i) = i;
            end
            Sig_Pix_Indicator_STASD = Sig_Pix_Indicator_STASD + ~isempty(RF_Ident.STASD_ASP_Num_RF_pixels{i});
        end
        Colour_with_Sig_Pix_STASD_vec(isnan(Colour_with_Sig_Pix_STASD_vec)) = [];
        
        if Sig_Pix_Indicator_STASD == 0
            
                if p.RF_Type(1) == 1 % Box
                    RF_Ident.STASD_Box_FullRF_Num_pixels = 0;
                end
                if p.RF_Type(2) == 1 % All Significant Pixels
                    RF_Ident.STASD_ASP_FullRF_Num_pixels = 0;
                end
                if p.RF_Type(3) == 1 % Gaussian
                    RF_Ident.STASD_Gaus_FullRF_Num_pixels = 0;
                end
             
        elseif Sig_Pix_Indicator_STASD == 1
            
                if p.RF_Type(1) == 1 % Box
                    RF_Ident.STASD_Box_FullRF_coords     = RF_Ident.STASD_Box_RF_coords{Colour_with_Sig_Pix_STASD_vec};
                    RF_Ident.STASD_Box_FullRF_Num_rows   = RF_Ident.STASD_Box_Num_RF_rows{Colour_with_Sig_Pix_STASD_vec};
                    RF_Ident.STASD_Box_FullRF_Num_cols   = RF_Ident.STASD_Box_Num_RF_cols{Colour_with_Sig_Pix_STASD_vec};
                    RF_Ident.STASD_Box_FullRF_Num_pixels = RF_Ident.STASD_Box_Num_RF_pixels{Colour_with_Sig_Pix_STASD_vec};
                end
                if p.RF_Type(2) == 1 % All Significant Pixels
                    RF_Ident.STASD_ASP_FullRF_coords     = RF_Ident.STASD_ASP_RF_coords{Colour_with_Sig_Pix_STASD_vec};
                    RF_Ident.STASD_ASP_FullRF_Num_pixels = RF_Ident.STASD_ASP_Num_RF_pixels{Colour_with_Sig_Pix_STASD_vec};
                end
                if p.RF_Type(3) == 1 % Gaussian
                    RF_Ident.STASD_Gaus_FullRF_coords           = RF_Ident.STASD_Gaus_RF_coords{Colour_with_Sig_Pix_STASD_vec};
                    RF_Ident.STASD_Gaus_FullRF_Num_pixels       = RF_Ident.STASD_Gaus_Num_RF_pixels{Colour_with_Sig_Pix_STASD_vec};
                    RF_Ident.STASD_Gaus_FullRF_Gaussian_mean    = RF_Ident.STASD_Gaus_Gaussian_mean{Colour_with_Sig_Pix_STASD_vec};
                    RF_Ident.STASD_Gaus_FullRF_Gaussian_covar   = RF_Ident.STASD_Gaus_Gaussian_covar{Colour_with_Sig_Pix_STASD_vec};
                    RF_Ident.STASD_Gaus_FullRF_eval_descend_vec = RF_Ident.STASD_Gaus_eval_descend_vec{Colour_with_Sig_Pix_STASD_vec};
                    RF_Ident.STASD_Gaus_FullRF_evec_descend_mat = RF_Ident.STASD_Gaus_evec_descend_mat{Colour_with_Sig_Pix_STASD_vec};
                    RF_Ident.STASD_Gaus_FullRF_Angle_Major_Axis = RF_Ident.STASD_Gaus_Angle_Major_Axis(Colour_with_Sig_Pix_STASD_vec);
                    RF_Ident.STASD_Gaus_FullRF_Angle_Minor_Axis = RF_Ident.STASD_Gaus_Angle_Minor_Axis(Colour_with_Sig_Pix_STASD_vec);
                    RF_Ident.STASD_Gaus_FullRF_Axis_eval_Ratio  = RF_Ident.STASD_Gaus_Axis_eval_Ratio(Colour_with_Sig_Pix_STASD_vec);
                    RF_Ident.STASD_Gaus_FullRF_Rel_Ellipse_Area = RF_Ident.STASD_Gaus_Rel_Ellipse_Area(Colour_with_Sig_Pix_STASD_vec);
                    RF_Ident.STASD_Gaus_FullRF_Abs_Ellipse_Area = RF_Ident.STASD_Gaus_Abs_Ellipse_Area(Colour_with_Sig_Pix_STASD_vec);
                    RF_Ident.STASD_Gaus_FullRF_gmPDF            = RF_Ident.STASD_Gaus_gmPDF{Colour_with_Sig_Pix_STASD_vec};
                    RF_Ident.STASD_Gaus_FullRF_Thresh_height    = RF_Ident.STASD_Gaus_Thresh_height{Colour_with_Sig_Pix_STASD_vec};
                end
            
        else % Sig_Pix_Indicator_STASD > 1
            
            if     p.RF_Type(1) == 1 % Box
                
                STASD_Box_RF_coords_All_centre = [];
                for i = 1:p.Spectral_Dim
                    STASD_Box_RF_coords_All_centre    = [STASD_Box_RF_coords_All_centre;RF_Ident.STASD_Box_RF_coords_centre{i}];
                end
                STASD_Box_RF_coords_All_centre        = unique(STASD_Box_RF_coords_All_centre,'rows');
                % If draw box around all sig pix.
                % [RF_Ident.STASD_Box_FullRF_coords,...
                %  RF_Ident.STASD_Box_FullRF_Num_rows,...
                %  RF_Ident.STASD_Box_FullRF_Num_cols,...
                %  RF_Ident.STASD_Box_FullRF_Num_pixels] = Box_RF_fn_alt(RF_Ident.STASD_ASP_FullRF_coords,p);
                % If draw box around best pix for each colour (do this as more consistent with single colour boxes).
                [RF_Ident.STASD_Box_FullRF_coords,...
                    RF_Ident.STASD_Box_FullRF_Num_rows,...
                    RF_Ident.STASD_Box_FullRF_Num_cols,...
                    RF_Ident.STASD_Box_FullRF_Num_pixels] = Box_RF_fn_alt_v3(STASD_Box_RF_coords_All_centre,p); % Mod 18,12,2020 (was Box_RF_fn_alt_v2)
                
            end
            
            if p.RF_Type(2) == 1     % All Significant Pixels
            
                % I began to write this to deal with the case where Gauss is
                % selected but ASP isn't. I will leave it for now, perhaps
                % writing the GUI such that ASP has to be chosen if Gauss
                % is or just removing the ASP variables from RF_ident at
                % the end of this code.
%                 STASD_ASP_RF_coords     = RF_Ident.STASD_ASP_RF_coords;
%                 STASD_ASP_Num_RF_pixels = RF_Ident.STASD_ASP_Num_RF_pixels;
%                 
%             elseif p.RF_Type(2) == 0 && p.RF_Type(3) == 1
%                 
%                 
%                 
%             end
%          
%             if (p.RF_Type(2) == 1) || (p.RF_Type(2) == 0 && p.RF_Type(3) == 1)
            
                % Old code from v3
%                 RF_Ident.STASD_ASP_FullRF_coords  = [];
%                 for i = 1:p.Spectral_Dim
%                     RF_Ident.STASD_ASP_FullRF_coords  = [RF_Ident.STASD_ASP_FullRF_coords;RF_Ident.STASD_ASP_RF_coords{i}];
%                 end
%                 RF_Ident.STASD_ASP_FullRF_coords      = unique(RF_Ident.STASD_ASP_FullRF_coords,'rows');
%                 RF_Ident.STASD_ASP_FullRF_Num_pixels  = size(RF_Ident.STASD_ASP_FullRF_coords,1);
                

                % Find max STA-SD for each colour
                max_STA_SD_val_vec = NaN(p.Spectral_Dim,1);
                for i = 1:Sig_Pix_Indicator_STASD
                    max_loop = max(STA_SD(RF_Ident.STASD_ASP_RF_coords{Colour_with_Sig_Pix_STASD_vec(i)}(:,1),RF_Ident.STASD_ASP_RF_coords{Colour_with_Sig_Pix_STASD_vec(i)}(:,2),Colour_with_Sig_Pix_STASD_vec(i)),[],'all'); 
                    max_loop = max_loop(1);
                    max_STA_SD_val_vec(Colour_with_Sig_Pix_STASD_vec(i)) = max_loop;
                end
                
                % Find the colour with the best STA-SD (if multiple choose first)
                Best_STASD_Colour  = find(max_STA_SD_val_vec == max(max_STA_SD_val_vec));
                if length(Best_STASD_Colour) > 1
                    Best_STASD_Colour = Best_STASD_Colour(1);
                end
                
                % Add pixels from best colour RF
                RF_Ident.STASD_ASP_FullRF_coords = RF_Ident.STASD_ASP_RF_coords{Best_STASD_Colour};
                
                % Find coordinates of ASP RF pixels and their neighbours
                STASD_ASP_Pixel_Neighbour_Coords = NaN(8*RF_Ident.STASD_ASP_Num_RF_pixels{Best_STASD_Colour},2);
                index_count = 1;
                for i = 1:RF_Ident.STASD_ASP_Num_RF_pixels{Best_STASD_Colour}
                    [Neigh_coords_loop,Num_Neigh_pixels_loop] = Neighbour_fn(RF_Ident.STASD_ASP_RF_coords{Best_STASD_Colour}(i,1),RF_Ident.STASD_ASP_RF_coords{Best_STASD_Colour}(i,2),p);
                    STASD_ASP_Pixel_Neighbour_Coords(index_count:index_count+Num_Neigh_pixels_loop-1,:) = Neigh_coords_loop;
                    index_count = index_count + Num_Neigh_pixels_loop;
                end
                STASD_ASP_Pixel_Neighbour_Coords(isnan(STASD_ASP_Pixel_Neighbour_Coords(:,1)),:) = []; % remove NaN
                STASD_ASP_Pixel_and_Neighbour_Coords = [STASD_ASP_Pixel_Neighbour_Coords;RF_Ident.STASD_ASP_RF_coords{Best_STASD_Colour}]; % append ASP RF coords
                STASD_ASP_Pixel_and_Neighbour_Coords = unique(STASD_ASP_Pixel_and_Neighbour_Coords,'rows'); % remove repeated coords
                
                % Find which RFs of other colours overlap with strongest RF or its neighbouring pixels
                remaining_colours_vec = Colour_with_Sig_Pix_STASD_vec;
                remaining_colours_vec(remaining_colours_vec==Best_STASD_Colour) = [];
                length_remaining_colours_vec = length(remaining_colours_vec);
                for i = 1:length_remaining_colours_vec
                    intersection_loop = intersect(RF_Ident.STASD_ASP_RF_coords{remaining_colours_vec(i)},STASD_ASP_Pixel_and_Neighbour_Coords,'rows');
                    if ~isempty(intersection_loop)
                        RF_Ident.STASD_ASP_FullRF_coords = [RF_Ident.STASD_ASP_FullRF_coords;RF_Ident.STASD_ASP_RF_coords{remaining_colours_vec(i)}];
                    end
                end
                RF_Ident.STASD_ASP_FullRF_coords      = unique(RF_Ident.STASD_ASP_FullRF_coords,'rows'); % remove repeated coords
                RF_Ident.STASD_ASP_FullRF_Num_pixels  = size(RF_Ident.STASD_ASP_FullRF_coords,1);
                %NB will need to add condition/rewrite in Gaussian section
                %below where use diff variables in ASP (i.e. when calc Gaussian but not ASP so ASP output has different name)
                
            end
            
            if p.RF_Type(3) == 1     % Gaussian
                
                [RF_Ident.STASD_Gaus_FullRF_coords,...
                    RF_Ident.STASD_Gaus_FullRF_Num_pixels,...
                    RF_Ident.STASD_Gaus_FullRF_Gaussian_mean,...
                    RF_Ident.STASD_Gaus_FullRF_Gaussian_covar,...
                    RF_Ident.STASD_Gaus_FullRF_eval_descend_vec,...
                    RF_Ident.STASD_Gaus_FullRF_evec_descend_mat,...
                    RF_Ident.STASD_Gaus_FullRF_Angle_Major_Axis,...
                    RF_Ident.STASD_Gaus_FullRF_Angle_Minor_Axis,...
                    RF_Ident.STASD_Gaus_FullRF_Axis_eval_Ratio,...
                    RF_Ident.STASD_Gaus_FullRF_Rel_Ellipse_Area,...
                    RF_Ident.STASD_Gaus_FullRF_Abs_Ellipse_Area,...
                    RF_Ident.STASD_Gaus_FullRF_gmPDF,...
                    RF_Ident.STASD_Gaus_FullRF_Thresh_height] = Gaussian_RF_fn_v3(RF_Ident.STASD_ASP_FullRF_coords,RF_Ident.STASD_ASP_FullRF_Num_pixels,p);
                
            end
            
        end
        
    end
    
    
    if     p.RF_Ident_Meth_vec(2) == 1 % LC
        
        % Check if there are any significant pixels across all of the colours
        Sig_Pix_Indicator_LC    = 0;
        Colour_with_Sig_Pix_LC_vec = NaN(p.Spectral_Dim,1);
        for i = 1:p.Spectral_Dim
            if ~isempty(RF_Ident.LC_ASP_Num_RF_pixels{i})
                Colour_with_Sig_Pix_LC_vec(i) = i;
            end
            Sig_Pix_Indicator_LC    = Sig_Pix_Indicator_LC    + ~isempty(RF_Ident.LC_ASP_Num_RF_pixels{i});
        end
        Colour_with_Sig_Pix_LC_vec(isnan(Colour_with_Sig_Pix_LC_vec)) = [];
        
        if Sig_Pix_Indicator_LC == 0
            
            if p.RF_Type(1) == 1 % Box
                RF_Ident.LC_Box_FullRF_Num_pixels = 0;
            end
            if p.RF_Type(2) == 1 % All Significant Pixels
                RF_Ident.LC_ASP_FullRF_Num_pixels = 0;
            end
            if p.RF_Type(3) == 1 % Gaussian
                RF_Ident.LC_Gaus_FullRF_Num_pixels = 0;
            end
            
        elseif Sig_Pix_Indicator_LC == 1
            
            if p.RF_Type(1) == 1 % Box
                RF_Ident.LC_Box_FullRF_coords     = RF_Ident.LC_Box_RF_coords{Colour_with_Sig_Pix_LC_vec};
                RF_Ident.LC_Box_FullRF_Num_rows   = RF_Ident.LC_Box_Num_RF_rows{Colour_with_Sig_Pix_LC_vec};
                RF_Ident.LC_Box_FullRF_Num_cols   = RF_Ident.LC_Box_Num_RF_cols{Colour_with_Sig_Pix_LC_vec};
                RF_Ident.LC_Box_FullRF_Num_pixels = RF_Ident.LC_Box_Num_RF_pixels{Colour_with_Sig_Pix_LC_vec};
            end
            if p.RF_Type(2) == 1 % All Significant Pixels
                RF_Ident.LC_ASP_FullRF_coords     = RF_Ident.LC_ASP_RF_coords{Colour_with_Sig_Pix_LC_vec};
                RF_Ident.LC_ASP_FullRF_Num_pixels = RF_Ident.LC_ASP_Num_RF_pixels{Colour_with_Sig_Pix_LC_vec};
            end
            if p.RF_Type(3) == 1 % Gaussian
                RF_Ident.LC_Gaus_FullRF_coords           = RF_Ident.LC_Gaus_RF_coords{Colour_with_Sig_Pix_LC_vec};
                RF_Ident.LC_Gaus_FullRF_Num_pixels       = RF_Ident.LC_Gaus_Num_RF_pixels{Colour_with_Sig_Pix_LC_vec};
                RF_Ident.LC_Gaus_FullRF_Gaussian_mean    = RF_Ident.LC_Gaus_Gaussian_mean{Colour_with_Sig_Pix_LC_vec};
                RF_Ident.LC_Gaus_FullRF_Gaussian_covar   = RF_Ident.LC_Gaus_Gaussian_covar{Colour_with_Sig_Pix_LC_vec};
                RF_Ident.LC_Gaus_FullRF_eval_descend_vec = RF_Ident.LC_Gaus_eval_descend_vec{Colour_with_Sig_Pix_LC_vec};
                RF_Ident.LC_Gaus_FullRF_evec_descend_mat = RF_Ident.LC_Gaus_evec_descend_mat{Colour_with_Sig_Pix_LC_vec};
                RF_Ident.LC_Gaus_FullRF_Angle_Major_Axis = RF_Ident.LC_Gaus_Angle_Major_Axis(Colour_with_Sig_Pix_LC_vec);
                RF_Ident.LC_Gaus_FullRF_Angle_Minor_Axis = RF_Ident.LC_Gaus_Angle_Minor_Axis(Colour_with_Sig_Pix_LC_vec);
                RF_Ident.LC_Gaus_FullRF_Axis_eval_Ratio  = RF_Ident.LC_Gaus_Axis_eval_Ratio(Colour_with_Sig_Pix_LC_vec);
                RF_Ident.LC_Gaus_FullRF_Rel_Ellipse_Area = RF_Ident.LC_Gaus_Rel_Ellipse_Area(Colour_with_Sig_Pix_LC_vec);
                RF_Ident.LC_Gaus_FullRF_Abs_Ellipse_Area = RF_Ident.LC_Gaus_Abs_Ellipse_Area(Colour_with_Sig_Pix_LC_vec);
                RF_Ident.LC_Gaus_FullRF_gmPDF            = RF_Ident.LC_Gaus_gmPDF{Colour_with_Sig_Pix_LC_vec};
                RF_Ident.LC_Gaus_FullRF_Thresh_height    = RF_Ident.LC_Gaus_Thresh_height{Colour_with_Sig_Pix_LC_vec};
            end
            
        else % Sig_Pix_Indicator_LC > 1
            
            if     p.RF_Type(1) == 1 % Box
                 LC_Box_RF_coords_All_centre = [];
                for i = 1:p.Spectral_Dim
                    LC_Box_RF_coords_All_centre    = [LC_Box_RF_coords_All_centre;RF_Ident.LC_Box_RF_coords_centre{i}];
                end
                LC_Box_RF_coords_All_centre        = unique(LC_Box_RF_coords_All_centre,'rows');
                % If draw box around all sig pix.
                % [RF_Ident.LC_Box_FullRF_coords,...
                %  RF_Ident.LC_Box_FullRF_Num_rows,...
                %  RF_Ident.LC_Box_FullRF_Num_cols,...
                %  RF_Ident.LC_Box_FullRF_Num_pixels] = Box_RF_fn_alt(RF_Ident.LC_ASP_FullRF_coords,p);
                % If draw box around best pix for each colour (do this as more consistent with single colour boxes).
                [RF_Ident.LC_Box_FullRF_coords,...
                    RF_Ident.LC_Box_FullRF_Num_rows,...
                    RF_Ident.LC_Box_FullRF_Num_cols,...
                    RF_Ident.LC_Box_FullRF_Num_pixels] = Box_RF_fn_alt_v3(LC_Box_RF_coords_All_centre,p); % Mod 18,12,2020 (was Box_RF_fn_alt_v2)
            end
            
            if p.RF_Type(2) == 1     % All Significant Pixels
                
                % Old code
%                RF_Ident.LC_ASP_FullRF_coords  = [];
%                 for i = 1:p.Spectral_Dim
%                     RF_Ident.LC_ASP_FullRF_coords  = [RF_Ident.LC_ASP_FullRF_coords;RF_Ident.LC_ASP_RF_coords{i}];
%                 end
%                 RF_Ident.LC_ASP_FullRF_coords      = unique(RF_Ident.LC_ASP_FullRF_coords,'rows');
%                 RF_Ident.LC_ASP_FullRF_Num_pixels  = size(RF_Ident.LC_ASP_FullRF_coords,1);
                
                % Find max LC for each colour
                max_LC_val_vec = NaN(p.Spectral_Dim,1);
                for i = 1:Sig_Pix_Indicator_LC
                    max_loop = max(LCA_Pixel_covar(RF_Ident.LC_ASP_RF_coords{Colour_with_Sig_Pix_LC_vec(i)}(:,1),RF_Ident.LC_ASP_RF_coords{Colour_with_Sig_Pix_LC_vec(i)}(:,2),Colour_with_Sig_Pix_LC_vec(i)),[],'all');
                    max_loop = max_loop(1);
                    max_LC_val_vec(Colour_with_Sig_Pix_LC_vec(i)) = max_loop;
                end
                
                % Find the colour with the best LC (if multiple choose first)
                Best_LC_Colour  = find(max_LC_val_vec == max(max_LC_val_vec));
                if length(Best_LC_Colour) > 1
                    Best_LC_Colour = Best_LC_Colour(1);
                end
                
                % Add pixels from best colour RF
                RF_Ident.LC_ASP_FullRF_coords = RF_Ident.LC_ASP_RF_coords{Best_LC_Colour};
                
                % Find coordinates of ASP RF pixels and their neighbours
                LC_ASP_Pixel_Neighbour_Coords = NaN(8*RF_Ident.LC_ASP_Num_RF_pixels{Best_LC_Colour},2);
                index_count = 1;
                for i = 1:RF_Ident.LC_ASP_Num_RF_pixels{Best_LC_Colour}
                    [Neigh_coords_loop,Num_Neigh_pixels_loop] = Neighbour_fn(RF_Ident.LC_ASP_RF_coords{Best_LC_Colour}(i,1),RF_Ident.LC_ASP_RF_coords{Best_LC_Colour}(i,2),p);
                    LC_ASP_Pixel_Neighbour_Coords(index_count:index_count+Num_Neigh_pixels_loop-1,:) = Neigh_coords_loop;
                    index_count = index_count + Num_Neigh_pixels_loop;
                end
                LC_ASP_Pixel_Neighbour_Coords(isnan(LC_ASP_Pixel_Neighbour_Coords(:,1)),:) = []; % remove NaN
                LC_ASP_Pixel_and_Neighbour_Coords = [LC_ASP_Pixel_Neighbour_Coords;RF_Ident.LC_ASP_RF_coords{Best_LC_Colour}]; % append ASP RF coords
                LC_ASP_Pixel_and_Neighbour_Coords = unique(LC_ASP_Pixel_and_Neighbour_Coords,'rows'); % remove repeated coords
                
                % Find which RFs of other colours overlap with strongest RF or its neighbouring pixels
                remaining_colours_vec = Colour_with_Sig_Pix_LC_vec;
                remaining_colours_vec(remaining_colours_vec==Best_LC_Colour) = [];
                length_remaining_colours_vec = length(remaining_colours_vec);
                for i = 1:length_remaining_colours_vec
                    intersection_loop = intersect(RF_Ident.LC_ASP_RF_coords{remaining_colours_vec(i)},LC_ASP_Pixel_and_Neighbour_Coords,'rows');
                    if ~isempty(intersection_loop)
                        RF_Ident.LC_ASP_FullRF_coords = [RF_Ident.LC_ASP_FullRF_coords;RF_Ident.LC_ASP_RF_coords{remaining_colours_vec(i)}];
                    end
                end
                RF_Ident.LC_ASP_FullRF_coords      = unique(RF_Ident.LC_ASP_FullRF_coords,'rows'); % remove repeated coords
                RF_Ident.LC_ASP_FullRF_Num_pixels  = size(RF_Ident.LC_ASP_FullRF_coords,1);
                %NB will need to add condition/rewrite in Gaussian section
                %below where use diff variables in ASP (i.e. when calc Gaussian but not ASP so ASP output has different name)
                    
            end
            
            if p.RF_Type(3) == 1     % Gaussian
                [RF_Ident.LC_Gaus_FullRF_coords,...
                    RF_Ident.LC_Gaus_FullRF_Num_pixels,...
                    RF_Ident.LC_Gaus_FullRF_Gaussian_mean,...
                    RF_Ident.LC_Gaus_FullRF_Gaussian_covar,...
                    RF_Ident.LC_Gaus_FullRF_eval_descend_vec,...
                    RF_Ident.LC_Gaus_FullRF_evec_descend_mat,...
                    RF_Ident.LC_Gaus_FullRF_Angle_Major_Axis,...
                    RF_Ident.LC_Gaus_FullRF_Angle_Minor_Axis,...
                    RF_Ident.LC_Gaus_FullRF_Axis_eval_Ratio,...
                    RF_Ident.LC_Gaus_FullRF_Rel_Ellipse_Area,...
                    RF_Ident.LC_Gaus_FullRF_Abs_Ellipse_Area,...
                    RF_Ident.LC_Gaus_FullRF_gmPDF,...
                    RF_Ident.LC_Gaus_FullRF_Thresh_height] = Gaussian_RF_fn_v3(RF_Ident.LC_ASP_FullRF_coords,RF_Ident.LC_ASP_FullRF_Num_pixels,p);
                    
            end
        end
        
    end
    
    
    if     p.RF_Ident_Meth_vec(3) == 1 % MI
        
        % Check if there are any significant pixels across all of the colours
        Sig_Pix_Indicator_MI    = 0;
        Colour_with_Sig_Pix_MI_vec = NaN(p.Spectral_Dim,1);
        for i = 1:p.Spectral_Dim
            if ~isempty(RF_Ident.MI_ASP_Num_RF_pixels{i})
                Colour_with_Sig_Pix_MI_vec(i) = i;
            end
            Sig_Pix_Indicator_MI    = Sig_Pix_Indicator_MI    + ~isempty(RF_Ident.MI_ASP_Num_RF_pixels{i});
        end
        Colour_with_Sig_Pix_MI_vec(isnan(Colour_with_Sig_Pix_MI_vec)) = [];
        
        if Sig_Pix_Indicator_MI == 0
            
            if p.RF_Type(1) == 1 % Box
                RF_Ident.MI_Box_FullRF_Num_pixels = 0;
            end
            if p.RF_Type(2) == 1 % All Significant Pixels
                RF_Ident.MI_ASP_FullRF_Num_pixels = 0;
            end
            if p.RF_Type(3) == 1 % Gaussian
                RF_Ident.MI_Gaus_FullRF_Num_pixels = 0;
            end
            
        elseif Sig_Pix_Indicator_MI == 1
            
            if p.RF_Type(1) == 1 % Box
                RF_Ident.MI_Box_FullRF_coords     = RF_Ident.MI_Box_RF_coords{Colour_with_Sig_Pix_MI_vec};
                RF_Ident.MI_Box_FullRF_Num_rows   = RF_Ident.MI_Box_Num_RF_rows{Colour_with_Sig_Pix_MI_vec};
                RF_Ident.MI_Box_FullRF_Num_cols   = RF_Ident.MI_Box_Num_RF_cols{Colour_with_Sig_Pix_MI_vec};
                RF_Ident.MI_Box_FullRF_Num_pixels = RF_Ident.MI_Box_Num_RF_pixels{Colour_with_Sig_Pix_MI_vec};
            end
            if p.RF_Type(2) == 1 % All Significant Pixels
                RF_Ident.MI_ASP_FullRF_coords     = RF_Ident.MI_ASP_RF_coords{Colour_with_Sig_Pix_MI_vec};
                RF_Ident.MI_ASP_FullRF_Num_pixels = RF_Ident.MI_ASP_Num_RF_pixels{Colour_with_Sig_Pix_MI_vec};
            end
            if p.RF_Type(3) == 1 % Gaussian
                RF_Ident.MI_Gaus_FullRF_coords           = RF_Ident.MI_Gaus_RF_coords{Colour_with_Sig_Pix_MI_vec};
                RF_Ident.MI_Gaus_FullRF_Num_pixels       = RF_Ident.MI_Gaus_Num_RF_pixels{Colour_with_Sig_Pix_MI_vec};
                RF_Ident.MI_Gaus_FullRF_Gaussian_mean    = RF_Ident.MI_Gaus_Gaussian_mean{Colour_with_Sig_Pix_MI_vec};
                RF_Ident.MI_Gaus_FullRF_Gaussian_covar   = RF_Ident.MI_Gaus_Gaussian_covar{Colour_with_Sig_Pix_MI_vec};
                RF_Ident.MI_Gaus_FullRF_eval_descend_vec = RF_Ident.MI_Gaus_eval_descend_vec{Colour_with_Sig_Pix_MI_vec};
                RF_Ident.MI_Gaus_FullRF_evec_descend_mat = RF_Ident.MI_Gaus_evec_descend_mat{Colour_with_Sig_Pix_MI_vec};
                RF_Ident.MI_Gaus_FullRF_Angle_Major_Axis = RF_Ident.MI_Gaus_Angle_Major_Axis(Colour_with_Sig_Pix_MI_vec);
                RF_Ident.MI_Gaus_FullRF_Angle_Minor_Axis = RF_Ident.MI_Gaus_Angle_Minor_Axis(Colour_with_Sig_Pix_MI_vec);
                RF_Ident.MI_Gaus_FullRF_Axis_eval_Ratio  = RF_Ident.MI_Gaus_Axis_eval_Ratio(Colour_with_Sig_Pix_MI_vec);
                RF_Ident.MI_Gaus_FullRF_Rel_Ellipse_Area = RF_Ident.MI_Gaus_Rel_Ellipse_Area(Colour_with_Sig_Pix_MI_vec);
                RF_Ident.MI_Gaus_FullRF_Abs_Ellipse_Area = RF_Ident.MI_Gaus_Abs_Ellipse_Area(Colour_with_Sig_Pix_MI_vec);
                RF_Ident.MI_Gaus_FullRF_gmPDF            = RF_Ident.MI_Gaus_gmPDF{Colour_with_Sig_Pix_MI_vec};
                RF_Ident.MI_Gaus_FullRF_Thresh_height    = RF_Ident.MI_Gaus_Thresh_height{Colour_with_Sig_Pix_MI_vec};
            end
            
        else % Sig_Pix_Indicator_MI > 1
            
            if     p.RF_Type(1) == 1 % Box
                 MI_Box_RF_coords_All_centre = [];
                for i = 1:p.Spectral_Dim
                    MI_Box_RF_coords_All_centre    = [MI_Box_RF_coords_All_centre;RF_Ident.MI_Box_RF_coords_centre{i}];
                end
                MI_Box_RF_coords_All_centre        = unique(MI_Box_RF_coords_All_centre,'rows');
                % If draw box around all sig pix.
                % [RF_Ident.MI_Box_FullRF_coords,...
                %  RF_Ident.MI_Box_FullRF_Num_rows,...
                %  RF_Ident.MI_Box_FullRF_Num_cols,...
                %  RF_Ident.MI_Box_FullRF_Num_pixels] = Box_RF_fn_alt(RF_Ident.MI_ASP_FullRF_coords,p);
                % If draw box around best pix for each colour (do this as more consistent with single colour boxes).
                [RF_Ident.MI_Box_FullRF_coords,...
                    RF_Ident.MI_Box_FullRF_Num_rows,...
                    RF_Ident.MI_Box_FullRF_Num_cols,...
                    RF_Ident.MI_Box_FullRF_Num_pixels] = Box_RF_fn_alt_v3(MI_Box_RF_coords_All_centre,p); % Mod 18,12,2020 (was Box_RF_fn_alt_v2)
            end
            
            if p.RF_Type(2) == 1     % All Significant Pixels
                
                %Old Code
%                RF_Ident.MI_ASP_FullRF_coords  = [];
%                 for i = 1:p.Spectral_Dim
%                     RF_Ident.MI_ASP_FullRF_coords  = [RF_Ident.MI_ASP_FullRF_coords;RF_Ident.MI_ASP_RF_coords{i}];
%                 end
%                 RF_Ident.MI_ASP_FullRF_coords      = unique(RF_Ident.MI_ASP_FullRF_coords,'rows');
%                 RF_Ident.MI_ASP_FullRF_Num_pixels  = size(RF_Ident.MI_ASP_FullRF_coords,1);
                
                
                % Find max MI for each colour
                max_MI_val_vec = NaN(p.Spectral_Dim,1);
                for i = 1:Sig_Pix_Indicator_MI
                    max_loop = max(MI(RF_Ident.MI_ASP_RF_coords{Colour_with_Sig_Pix_MI_vec(i)}(:,1),RF_Ident.MI_ASP_RF_coords{Colour_with_Sig_Pix_MI_vec(i)}(:,2),Colour_with_Sig_Pix_MI_vec(i)),[],'all');
                    max_loop = max_loop(1);
                    max_MI_val_vec(Colour_with_Sig_Pix_MI_vec(i)) = max_loop;
                end
                
                % Find the colour with the best MI (if multiple choose first)
                Best_MI_Colour  = find(max_MI_val_vec == max(max_MI_val_vec));
                if length(Best_MI_Colour) > 1
                    Best_MI_Colour = Best_MI_Colour(1);
                end
                
                % Add pixels from best colour RF
                RF_Ident.MI_ASP_FullRF_coords = RF_Ident.MI_ASP_RF_coords{Best_MI_Colour};
                
                % Find coordinates of ASP RF pixels and their neighbours
                MI_ASP_Pixel_Neighbour_Coords = NaN(8*RF_Ident.MI_ASP_Num_RF_pixels{Best_MI_Colour},2);
                index_count = 1;
                for i = 1:RF_Ident.MI_ASP_Num_RF_pixels{Best_MI_Colour}
                    [Neigh_coords_loop,Num_Neigh_pixels_loop] = Neighbour_fn(RF_Ident.MI_ASP_RF_coords{Best_MI_Colour}(i,1),RF_Ident.MI_ASP_RF_coords{Best_MI_Colour}(i,2),p);
                    MI_ASP_Pixel_Neighbour_Coords(index_count:index_count+Num_Neigh_pixels_loop-1,:) = Neigh_coords_loop;
                    index_count = index_count + Num_Neigh_pixels_loop;
                end
                MI_ASP_Pixel_Neighbour_Coords(isnan(MI_ASP_Pixel_Neighbour_Coords(:,1)),:) = []; % remove NaN
                MI_ASP_Pixel_and_Neighbour_Coords = [MI_ASP_Pixel_Neighbour_Coords;RF_Ident.MI_ASP_RF_coords{Best_MI_Colour}]; % append ASP RF coords
                MI_ASP_Pixel_and_Neighbour_Coords = unique(MI_ASP_Pixel_and_Neighbour_Coords,'rows'); % remove repeated coords
                
                % Find which RFs of other colours overlap with strongest RF or its neighbouring pixels
                remaining_colours_vec = Colour_with_Sig_Pix_MI_vec;
                remaining_colours_vec(remaining_colours_vec==Best_MI_Colour) = [];
                length_remaining_colours_vec = length(remaining_colours_vec);
                for i = 1:length_remaining_colours_vec
                    intersection_loop = intersect(RF_Ident.MI_ASP_RF_coords{remaining_colours_vec(i)},MI_ASP_Pixel_and_Neighbour_Coords,'rows');
                    if ~isempty(intersection_loop)
                        RF_Ident.MI_ASP_FullRF_coords = [RF_Ident.MI_ASP_FullRF_coords;RF_Ident.MI_ASP_RF_coords{remaining_colours_vec(i)}];
                    end
                end
                RF_Ident.MI_ASP_FullRF_coords      = unique(RF_Ident.MI_ASP_FullRF_coords,'rows'); % remove repeated coords
                RF_Ident.MI_ASP_FullRF_Num_pixels  = size(RF_Ident.MI_ASP_FullRF_coords,1);
                %NB will need to add condition/rewrite in Gaussian section
                %below where use diff variables in ASP (i.e. when calc Gaussian but not ASP so ASP output has different name)
                
            end
            
            if p.RF_Type(3) == 1     % Gaussian
                [RF_Ident.MI_Gaus_FullRF_coords,...
                    RF_Ident.MI_Gaus_FullRF_Num_pixels,...
                    RF_Ident.MI_Gaus_FullRF_Gaussian_mean,...
                    RF_Ident.MI_Gaus_FullRF_Gaussian_covar,...
                    RF_Ident.MI_Gaus_FullRF_eval_descend_vec,...
                    RF_Ident.MI_Gaus_FullRF_evec_descend_mat,...
                    RF_Ident.MI_Gaus_FullRF_Angle_Major_Axis,...
                    RF_Ident.MI_Gaus_FullRF_Angle_Minor_Axis,...
                    RF_Ident.MI_Gaus_FullRF_Axis_eval_Ratio,...
                    RF_Ident.MI_Gaus_FullRF_Rel_Ellipse_Area,...
                    RF_Ident.MI_Gaus_FullRF_Abs_Ellipse_Area,...
                    RF_Ident.MI_Gaus_FullRF_gmPDF,...
                    RF_Ident.MI_Gaus_FullRF_Thresh_height] = Gaussian_RF_fn_v3(RF_Ident.MI_ASP_FullRF_coords,RF_Ident.MI_ASP_FullRF_Num_pixels,p);
            end
        end
        
    end
    
    
    if     p.RF_Ident_Meth_vec(4) == 1 % SC
        
        % Check if there are any significant pixels across all of the colours
        Sig_Pix_Indicator_SC    = 0;
        Colour_with_Sig_Pix_SC_vec = NaN(p.Spectral_Dim,1);
        for i = 1:p.Spectral_Dim
            if ~isempty(RF_Ident.SC_ASP_Num_RF_pixels{i})
                Colour_with_Sig_Pix_SC_vec(i) = i;
            end
            Sig_Pix_Indicator_SC    = Sig_Pix_Indicator_SC    + ~isempty(RF_Ident.SC_ASP_Num_RF_pixels{i});
        end
        Colour_with_Sig_Pix_SC_vec(isnan(Colour_with_Sig_Pix_SC_vec)) = [];
        
        if Sig_Pix_Indicator_SC == 0
            
            if p.RF_Type(1) == 1 % Box
                RF_Ident.SC_Box_FullRF_Num_pixels = 0;
            end
            if p.RF_Type(2) == 1 % All Significant Pixels
                RF_Ident.SC_ASP_FullRF_Num_pixels = 0;
            end
            if p.RF_Type(3) == 1 % Gaussian
                RF_Ident.SC_Gaus_FullRF_Num_pixels = 0;
            end
            
        elseif Sig_Pix_Indicator_SC == 1
            
            if p.RF_Type(1) == 1 % Box
                RF_Ident.SC_Box_FullRF_coords     = RF_Ident.SC_Box_RF_coords{Colour_with_Sig_Pix_SC_vec};
                RF_Ident.SC_Box_FullRF_Num_rows   = RF_Ident.SC_Box_Num_RF_rows{Colour_with_Sig_Pix_SC_vec};
                RF_Ident.SC_Box_FullRF_Num_cols   = RF_Ident.SC_Box_Num_RF_cols{Colour_with_Sig_Pix_SC_vec};
                RF_Ident.SC_Box_FullRF_Num_pixels = RF_Ident.SC_Box_Num_RF_pixels{Colour_with_Sig_Pix_SC_vec};
            end
            if p.RF_Type(2) == 1 % All Significant Pixels
                RF_Ident.SC_ASP_FullRF_coords     = RF_Ident.SC_ASP_RF_coords{Colour_with_Sig_Pix_SC_vec};
                RF_Ident.SC_ASP_FullRF_Num_pixels = RF_Ident.SC_ASP_Num_RF_pixels{Colour_with_Sig_Pix_SC_vec};
            end
            if p.RF_Type(3) == 1 % Gaussian
                RF_Ident.SC_Gaus_FullRF_coords           = RF_Ident.SC_Gaus_RF_coords{Colour_with_Sig_Pix_SC_vec};
                RF_Ident.SC_Gaus_FullRF_Num_pixels       = RF_Ident.SC_Gaus_Num_RF_pixels{Colour_with_Sig_Pix_SC_vec};
                RF_Ident.SC_Gaus_FullRF_Gaussian_mean    = RF_Ident.SC_Gaus_Gaussian_mean{Colour_with_Sig_Pix_SC_vec};
                RF_Ident.SC_Gaus_FullRF_Gaussian_covar   = RF_Ident.SC_Gaus_Gaussian_covar{Colour_with_Sig_Pix_SC_vec};
                RF_Ident.SC_Gaus_FullRF_eval_descend_vec = RF_Ident.SC_Gaus_eval_descend_vec{Colour_with_Sig_Pix_SC_vec};
                RF_Ident.SC_Gaus_FullRF_evec_descend_mat = RF_Ident.SC_Gaus_evec_descend_mat{Colour_with_Sig_Pix_SC_vec};
                RF_Ident.SC_Gaus_FullRF_Angle_Major_Axis = RF_Ident.SC_Gaus_Angle_Major_Axis(Colour_with_Sig_Pix_SC_vec);
                RF_Ident.SC_Gaus_FullRF_Angle_Minor_Axis = RF_Ident.SC_Gaus_Angle_Minor_Axis(Colour_with_Sig_Pix_SC_vec);
                RF_Ident.SC_Gaus_FullRF_Axis_eval_Ratio  = RF_Ident.SC_Gaus_Axis_eval_Ratio(Colour_with_Sig_Pix_SC_vec);
                RF_Ident.SC_Gaus_FullRF_Rel_Ellipse_Area = RF_Ident.SC_Gaus_Rel_Ellipse_Area(Colour_with_Sig_Pix_SC_vec);
                RF_Ident.SC_Gaus_FullRF_Abs_Ellipse_Area = RF_Ident.SC_Gaus_Abs_Ellipse_Area(Colour_with_Sig_Pix_SC_vec);
                RF_Ident.SC_Gaus_FullRF_gmPDF            = RF_Ident.SC_Gaus_gmPDF{Colour_with_Sig_Pix_SC_vec};
                RF_Ident.SC_Gaus_FullRF_Thresh_height    = RF_Ident.SC_Gaus_Thresh_height{Colour_with_Sig_Pix_SC_vec};
            end
            
        else % Sig_Pix_Indicator_SC > 1
            
            if     p.RF_Type(1) == 1 % Box
                SC_Box_RF_coords_All_centre = [];
                for i = 1:p.Spectral_Dim
                    SC_Box_RF_coords_All_centre    = [SC_Box_RF_coords_All_centre;RF_Ident.SC_Box_RF_coords_centre{i}];
                end
                SC_Box_RF_coords_All_centre        = unique(SC_Box_RF_coords_All_centre,'rows');
                % If draw box around all sig pix.
                % [RF_Ident.SC_Box_FullRF_coords,...
                %  RF_Ident.SC_Box_FullRF_Num_rows,...
                %  RF_Ident.SC_Box_FullRF_Num_cols,...
                %  RF_Ident.SC_Box_FullRF_Num_pixels] = Box_RF_fn_alt(RF_Ident.SC_ASP_FullRF_coords,p);
                % If draw box around best pix for each colour (do this as more consistent with single colour boxes).
                [RF_Ident.SC_Box_FullRF_coords,...
                    RF_Ident.SC_Box_FullRF_Num_rows,...
                    RF_Ident.SC_Box_FullRF_Num_cols,...
                    RF_Ident.SC_Box_FullRF_Num_pixels] = Box_RF_fn_alt_v3(SC_Box_RF_coords_All_centre,p); % Mod 18,12,2020 (was Box_RF_fn_alt_v2)
            end
            
            if p.RF_Type(2) == 1     % All Significant Pixels
                
                % Old code
%                 RF_Ident.SC_ASP_FullRF_coords  = [];
%                 for i = 1:p.Spectral_Dim
%                     RF_Ident.SC_ASP_FullRF_coords  = [RF_Ident.SC_ASP_FullRF_coords;RF_Ident.SC_ASP_RF_coords{i}];
%                 end
%                 RF_Ident.SC_ASP_FullRF_coords      = unique(RF_Ident.SC_ASP_FullRF_coords,'rows');
%                 RF_Ident.SC_ASP_FullRF_Num_pixels  = size(RF_Ident.SC_ASP_FullRF_coords,1);
                
                % Find min SC for each colour (min here unlike others which are max)
                min_SC_val_vec = NaN(p.Spectral_Dim,1);
                for i = 1:Sig_Pix_Indicator_SC
                    min_loop = min(SCA_Pixel_covar(RF_Ident.SC_ASP_RF_coords{Colour_with_Sig_Pix_SC_vec(i)}(:,1),RF_Ident.SC_ASP_RF_coords{Colour_with_Sig_Pix_SC_vec(i)}(:,2),Colour_with_Sig_Pix_SC_vec(i)),[],'all');
                    min_loop = min_loop(1);
                    min_SC_val_vec(Colour_with_Sig_Pix_SC_vec(i)) = min_loop;
                end
                
                % Find the colour with the best SC (if multiple choose first)
                Best_SC_Colour  = find(min_SC_val_vec == min(min_SC_val_vec));
                if length(Best_SC_Colour) > 1
                    Best_SC_Colour = Best_SC_Colour(1);
                end
                
                % Add pixels from best colour RF
                RF_Ident.SC_ASP_FullRF_coords = RF_Ident.SC_ASP_RF_coords{Best_SC_Colour};
                
                % Find coordinates of ASP RF pixels and their neighbours
                SC_ASP_Pixel_Neighbour_Coords = NaN(8*RF_Ident.SC_ASP_Num_RF_pixels{Best_SC_Colour},2);
                index_count = 1;
                for i = 1:RF_Ident.SC_ASP_Num_RF_pixels{Best_SC_Colour}
                    [Neigh_coords_loop,Num_Neigh_pixels_loop] = Neighbour_fn(RF_Ident.SC_ASP_RF_coords{Best_SC_Colour}(i,1),RF_Ident.SC_ASP_RF_coords{Best_SC_Colour}(i,2),p);
                    SC_ASP_Pixel_Neighbour_Coords(index_count:index_count+Num_Neigh_pixels_loop-1,:) = Neigh_coords_loop;
                    index_count = index_count + Num_Neigh_pixels_loop;
                end
                SC_ASP_Pixel_Neighbour_Coords(isnan(SC_ASP_Pixel_Neighbour_Coords(:,1)),:) = []; % remove NaN
                SC_ASP_Pixel_and_Neighbour_Coords = [SC_ASP_Pixel_Neighbour_Coords;RF_Ident.SC_ASP_RF_coords{Best_SC_Colour}]; % append ASP RF coords
                SC_ASP_Pixel_and_Neighbour_Coords = unique(SC_ASP_Pixel_and_Neighbour_Coords,'rows'); % remove repeated coords
                
                % Find which RFs of other colours overlap with strongest RF or its neighbouring pixels
                remaining_colours_vec = Colour_with_Sig_Pix_SC_vec;
                remaining_colours_vec(remaining_colours_vec==Best_SC_Colour) = [];
                length_remaining_colours_vec = length(remaining_colours_vec);
                for i = 1:length_remaining_colours_vec
                    intersection_loop = intersect(RF_Ident.SC_ASP_RF_coords{remaining_colours_vec(i)},SC_ASP_Pixel_and_Neighbour_Coords,'rows');
                    if ~isempty(intersection_loop)
                        RF_Ident.SC_ASP_FullRF_coords = [RF_Ident.SC_ASP_FullRF_coords;RF_Ident.SC_ASP_RF_coords{remaining_colours_vec(i)}];
                    end
                end
                RF_Ident.SC_ASP_FullRF_coords      = unique(RF_Ident.SC_ASP_FullRF_coords,'rows'); % remove repeated coords
                RF_Ident.SC_ASP_FullRF_Num_pixels  = size(RF_Ident.SC_ASP_FullRF_coords,1);
                %NB will need to add condition/rewrite in Gaussian section
                %below where use diff variables in ASP (i.e. when calc Gaussian but not ASP so ASP output has different name)
                
            end
            
            if p.RF_Type(3) == 1     % Gaussian
                [RF_Ident.SC_Gaus_FullRF_coords,...
                    RF_Ident.SC_Gaus_FullRF_Num_pixels,...
                    RF_Ident.SC_Gaus_FullRF_Gaussian_mean,...
                    RF_Ident.SC_Gaus_FullRF_Gaussian_covar,...
                    RF_Ident.SC_Gaus_FullRF_eval_descend_vec,...
                    RF_Ident.SC_Gaus_FullRF_evec_descend_mat,...
                    RF_Ident.SC_Gaus_FullRF_Angle_Major_Axis,...
                    RF_Ident.SC_Gaus_FullRF_Angle_Minor_Axis,...
                    RF_Ident.SC_Gaus_FullRF_Axis_eval_Ratio,...
                    RF_Ident.SC_Gaus_FullRF_Rel_Ellipse_Area,...
                    RF_Ident.SC_Gaus_FullRF_Abs_Ellipse_Area,...
                    RF_Ident.SC_Gaus_FullRF_gmPDF,...
                    RF_Ident.SC_Gaus_FullRF_Thresh_height] = Gaussian_RF_fn_v3(RF_Ident.SC_ASP_FullRF_coords,RF_Ident.SC_ASP_FullRF_Num_pixels,p);    
            end
        end
        
    end
    
end

%% RF STA, STASD, Temporal STA, Stixel, Pixel and Temporal Stixel Covariances Code % PAR Mod 09,10,2020 (everything below here)


% Find STA, STASD, Temporal STA, Stixel, Pixel and Temporal Stixel Covariances (spectral)
for i = 1:p.Spectral_Dim
    
    if     p.RF_Ident_Meth_vec(1) == 1 % STA-SD
        
        % Box RF STA, STASD and Temp. STA
        if     p.RF_Type(1) == 1 % Box
            if ~isempty(RF_Ident.STASD_Box_Num_RF_pixels{i})
                RF_Ident.STASD_Box_RF_STA{i}    = NaN(RF_Ident.STASD_Box_Num_RF_pixels{i},p.Num_STE_bins);
                RF_Ident.STASD_Box_RF_STA_SD{i} = NaN(RF_Ident.STASD_Box_Num_RF_pixels{i},1);
                for j = 1:RF_Ident.STASD_Box_Num_RF_pixels{i}
                    RF_Ident.STASD_Box_RF_STA{i}(j,:)  =    STA(RF_Ident.STASD_Box_RF_coords{i}(j,1),RF_Ident.STASD_Box_RF_coords{i}(j,2),:,i);
                    RF_Ident.STASD_Box_RF_STA_SD{i}(j) = STA_SD(RF_Ident.STASD_Box_RF_coords{i}(j,1),RF_Ident.STASD_Box_RF_coords{i}(j,2),i);
                end
                
                if p.Time_Plot_Choice == 1     % mean
                    RF_Ident.STASD_Box_RF_STA_time{i} = NaN(p.Num_STE_bins,1);
                    for j = 1:p.Num_STE_bins
                        RF_Ident.STASD_Box_RF_STA_time{i}(j) = mean(RF_Ident.STASD_Box_RF_STA{i}(:,j));
                    end
                else % p.Time_Plot_Choice == 2 % most significant pixel
                    RF_Ident.STASD_Box_RF_STA_time{i} = squeeze(STA(RF_Ident.STASD_Box_RF_coords_centre{i}(1),RF_Ident.STASD_Box_RF_coords_centre{i}(2),:,i));
                end
            end
        end
        
        % ASP RF STA, STASD and Temp. STA
        if p.RF_Type(2) == 1     % All Significant Pixels
            if ~isempty(RF_Ident.STASD_ASP_Num_RF_pixels{i})
                RF_Ident.STASD_ASP_RF_STA{i}    = NaN(RF_Ident.STASD_ASP_Num_RF_pixels{i},p.Num_STE_bins);
                RF_Ident.STASD_ASP_RF_STA_SD{i} = NaN(RF_Ident.STASD_ASP_Num_RF_pixels{i},1);
                for j = 1:RF_Ident.STASD_ASP_Num_RF_pixels{i}
                    RF_Ident.STASD_ASP_RF_STA{i}(j,:)  =    STA(RF_Ident.STASD_ASP_RF_coords{i}(j,1),RF_Ident.STASD_ASP_RF_coords{i}(j,2),:,i);
                    RF_Ident.STASD_ASP_RF_STA_SD{i}(j) = STA_SD(RF_Ident.STASD_ASP_RF_coords{i}(j,1),RF_Ident.STASD_ASP_RF_coords{i}(j,2),i);
                end
                
                if p.Time_Plot_Choice == 1     % mean
                    RF_Ident.STASD_ASP_RF_STA_time{i} = NaN(p.Num_STE_bins,1);
                    for j = 1:p.Num_STE_bins
                        RF_Ident.STASD_ASP_RF_STA_time{i}(j) = mean(RF_Ident.STASD_ASP_RF_STA{i}(:,j));
                    end
                else % p.Time_Plot_Choice == 2 % most significant pixel
                    RF_Ident.STASD_ASP_RF_STA_time{i} = squeeze(STA(RF_Ident.STASD_ASP_RF_coords{i}(1,1),RF_Ident.STASD_ASP_RF_coords{i}(1,2),:,i));
                end
            end
        end
        
        % Gaus RF STA, STASD and Temp. STA
        if p.RF_Type(3) == 1     % Gaussian
            if ~isempty(RF_Ident.STASD_Gaus_Num_RF_pixels{i})
                RF_Ident.STASD_Gaus_RF_STA{i}    = NaN(RF_Ident.STASD_Gaus_Num_RF_pixels{i},p.Num_STE_bins);
                RF_Ident.STASD_Gaus_RF_STA_SD{i} = NaN(RF_Ident.STASD_Gaus_Num_RF_pixels{i},1);
                for j = 1:RF_Ident.STASD_Gaus_Num_RF_pixels{i}
                    RF_Ident.STASD_Gaus_RF_STA{i}(j,:)  =    STA(RF_Ident.STASD_Gaus_RF_coords{i}(j,1),RF_Ident.STASD_Gaus_RF_coords{i}(j,2),:,i);
                    RF_Ident.STASD_Gaus_RF_STA_SD{i}(j) = STA_SD(RF_Ident.STASD_Gaus_RF_coords{i}(j,1),RF_Ident.STASD_Gaus_RF_coords{i}(j,2),i);
                end
                
                if p.Time_Plot_Choice == 1     % mean
                    RF_Ident.STASD_Gaus_RF_STA_time{i} = NaN(p.Num_STE_bins,1);
                    for j = 1:p.Num_STE_bins
                        RF_Ident.STASD_Gaus_RF_STA_time{i}(j) = mean(RF_Ident.STASD_Gaus_RF_STA{i}(:,j));
                    end
                else % p.Time_Plot_Choice == 2 % most significant pixel (use ASP coords below as these hold for Gaus too)
                    RF_Ident.STASD_Gaus_RF_STA_time{i} = squeeze(STA(RF_Ident.STASD_ASP_RF_coords{i}(1,1),RF_Ident.STASD_ASP_RF_coords{i}(1,2),:,i));
                end
            end
        end
        
    end
    
    
    if     p.RF_Ident_Meth_vec(4) == 1 % SC
        
        % Box RF Stix. Cov., Pix. Cov. and Temp. Stix. Cov.
        if     p.RF_Type(1) == 1 % Box
            if ~isempty(RF_Ident.SC_Box_Num_RF_pixels{i})
                RF_Ident.SC_Box_RF_Stixel_covar{i}    = NaN(RF_Ident.SC_Box_Num_RF_pixels{i},p.Num_STE_bins);
                RF_Ident.SC_Box_RF_Pixel_covar{i} = NaN(RF_Ident.SC_Box_Num_RF_pixels{i},1);
                for j = 1:RF_Ident.SC_Box_Num_RF_pixels{i}
                    RF_Ident.SC_Box_RF_Stixel_covar{i}(j,:) = SCA_Stixel_covar(RF_Ident.SC_Box_RF_coords{i}(j,1),RF_Ident.SC_Box_RF_coords{i}(j,2),:,i);
                    RF_Ident.SC_Box_RF_Pixel_covar{i}(j)    =  SCA_Pixel_covar(RF_Ident.SC_Box_RF_coords{i}(j,1),RF_Ident.SC_Box_RF_coords{i}(j,2),i);
                end
                
                if p.Time_Plot_Choice == 1     % mean
                    RF_Ident.SC_Box_RF_Stixel_covar_time{i} = NaN(p.Num_STE_bins,1);
                    for j = 1:p.Num_STE_bins
                        RF_Ident.SC_Box_RF_Stixel_covar_time{i}(j) = mean(RF_Ident.SC_Box_RF_Stixel_covar{i}(:,j));
                    end
                else % p.Time_Plot_Choice == 2 % most significant pixel
                    RF_Ident.SC_Box_RF_Stixel_covar_time{i} = squeeze(SCA_Stixel_covar(RF_Ident.SC_Box_RF_coords_centre{i}(1),RF_Ident.SC_Box_RF_coords_centre{i}(2),:,i));
                end
            end
        end
        
        % ASP RF Stix. Cov., Pix. Cov. and Temp. Stix. Cov.
        if p.RF_Type(2) == 1     % All Significant Pixels
            if ~isempty(RF_Ident.SC_ASP_Num_RF_pixels{i})
                RF_Ident.SC_ASP_RF_Stixel_covar{i}    = NaN(RF_Ident.SC_ASP_Num_RF_pixels{i},p.Num_STE_bins);
                RF_Ident.SC_ASP_RF_Pixel_covar{i} = NaN(RF_Ident.SC_ASP_Num_RF_pixels{i},1);
                for j = 1:RF_Ident.SC_ASP_Num_RF_pixels{i}
                    RF_Ident.SC_ASP_RF_Stixel_covar{i}(j,:) = SCA_Stixel_covar(RF_Ident.SC_ASP_RF_coords{i}(j,1),RF_Ident.SC_ASP_RF_coords{i}(j,2),:,i);
                    RF_Ident.SC_ASP_RF_Pixel_covar{i}(j)    =  SCA_Pixel_covar(RF_Ident.SC_ASP_RF_coords{i}(j,1),RF_Ident.SC_ASP_RF_coords{i}(j,2),i);
                end
                
                if p.Time_Plot_Choice == 1     % mean
                    RF_Ident.SC_ASP_RF_Stixel_covar_time{i} = NaN(p.Num_STE_bins,1);
                    for j = 1:p.Num_STE_bins
                        RF_Ident.SC_ASP_RF_Stixel_covar_time{i}(j) = mean(RF_Ident.SC_ASP_RF_Stixel_covar{i}(:,j));
                    end
                else % p.Time_Plot_Choice == 2 % most significant pixel
                    RF_Ident.SC_ASP_RF_Stixel_covar_time{i} = squeeze(SCA_Stixel_covar(RF_Ident.SC_ASP_RF_coords{i}(1,1),RF_Ident.SC_ASP_RF_coords{i}(1,2),:,i));
                end
            end
        end
        
        % Gaus RF Stix. Cov., Pix. Cov. and Temp. Stix. Cov.
        if p.RF_Type(3) == 1     % Gaussian
            if ~isempty(RF_Ident.SC_Gaus_Num_RF_pixels{i})
                RF_Ident.SC_Gaus_RF_Stixel_covar{i}    = NaN(RF_Ident.SC_Gaus_Num_RF_pixels{i},p.Num_STE_bins);
                RF_Ident.SC_Gaus_RF_Pixel_covar{i} = NaN(RF_Ident.SC_Gaus_Num_RF_pixels{i},1);
                for j = 1:RF_Ident.SC_Gaus_Num_RF_pixels{i}
                    RF_Ident.SC_Gaus_RF_Stixel_covar{i}(j,:) = SCA_Stixel_covar(RF_Ident.SC_Gaus_RF_coords{i}(j,1),RF_Ident.SC_Gaus_RF_coords{i}(j,2),:,i);
                    RF_Ident.SC_Gaus_RF_Pixel_covar{i}(j)    =  SCA_Pixel_covar(RF_Ident.SC_Gaus_RF_coords{i}(j,1),RF_Ident.SC_Gaus_RF_coords{i}(j,2),i);
                end
                
                if p.Time_Plot_Choice == 1     % mean
                    RF_Ident.SC_Gaus_RF_Stixel_covar_time{i} = NaN(p.Num_STE_bins,1);
                    for j = 1:p.Num_STE_bins
                        RF_Ident.SC_Gaus_RF_Stixel_covar_time{i}(j) = mean(RF_Ident.SC_Gaus_RF_Stixel_covar{i}(:,j));
                    end
                else % p.Time_Plot_Choice == 2 % most significant pixel (use ASP coords below as these hold for Gaus too)
                    RF_Ident.SC_Gaus_RF_Stixel_covar_time{i} = squeeze(SCA_Stixel_covar(RF_Ident.SC_ASP_RF_coords{i}(1,1),RF_Ident.SC_ASP_RF_coords{i}(1,2),:,i));
                end
            end
        end
        
    end
    
end

if p.Spectral_Dim > 1
    
    %%% Choose colour with most significant pixel for Full RFs
    
    if     p.RF_Ident_Meth_vec(1) == 1 % STA-SD
        % Box RF Temporal STA (FullRF)
        if     p.RF_Type(1) == 1 % Box
            if (p.Time_Plot_Choice == 2) && (RF_Ident.STASD_Box_FullRF_Num_pixels>0) % (Mod 18,12,2020, was: '...&& (~isempty(RF_Ident.STASD_Box_FullRF_Num_pixels))')
                loop_maxSTASD = NaN(p.Spectral_Dim,1);
                for i = 1:p.Spectral_Dim
                    if ~isempty(RF_Ident.STASD_Box_Num_RF_pixels{i})
                        loop_maxSTASD(i) = max(RF_Ident.STASD_Box_RF_STA_SD{i});
                    end
                end
                [~,maxSTASD_Box_index] = max(loop_maxSTASD); % colour with highest STASD
                maxSTASD_Box_index     = maxSTASD_Box_index(1);  % if there are multiple colours with co-highest, take the first
            end
        end
        % ASP RF Temporal STA (FullRF)
        if p.RF_Type(2) == 1     % All Significant Pixels
            if (p.Time_Plot_Choice == 2) && (RF_Ident.STASD_ASP_FullRF_Num_pixels > 0) % (Mod 18,12,2020, was: '...&& (~isempty(RF_Ident.STASD_ASP_FullRF_Num_pixels))')
                loop_maxSTASD = NaN(p.Spectral_Dim,1);
                for i = 1:p.Spectral_Dim
                    if ~isempty(RF_Ident.STASD_ASP_Num_RF_pixels{i})
                        loop_maxSTASD(i) = max(RF_Ident.STASD_ASP_RF_STA_SD{i});
                    end
                end
                [~,maxSTASD_ASP_index] = max(loop_maxSTASD); % colour with highest STASD
                maxSTASD_ASP_index     = maxSTASD_ASP_index(1);  % if there are multiple colours with co-highest, take the first
            end
        end
        % Don't need below as use ASP coords
        %     % Gaus RF Temporal STA (FullRF)
        %     if p.RF_Type(3) == 1     % Gaussian
        %         if (p.Time_Plot_Choice == 2) && (~isempty(RF_Ident.STASD_Gaus_FullRF_Num_pixels))
        %             loop_maxSTASD = NaN(p.Spectral_Dim,1);
        %             for i = 1:p.Spectral_Dim
        %                 if ~isempty(RF_Ident.STASD_Gaus_Num_RF_pixels{i})
        %                     loop_maxSTASD(i) = max(RF_Ident.STASD_Gaus_RF_STA_SD{i});
        %                 end
        %             end
        %             [~,maxSTASD_index] = max(loop_maxSTASD); % colour with highest STASD
        %             maxSTASD_index     = maxSTASD_index(1);  % if there are multiple colours with co-highest, take the first
        %         end
        %     end
    end
    
    if     p.RF_Ident_Meth_vec(4) == 1 % SC
        % Box RF Temporal Stixel Covar (FullRF)
        if     p.RF_Type(1) == 1 % Box
            if (p.Time_Plot_Choice == 2) && (RF_Ident.SC_Box_FullRF_Num_pixels > 0) % (Mod 18,12,2020, was: '...&& (~isempty(RF_Ident.SC_Box_FullRF_Num_pixels))')
                loop_minSC = NaN(p.Spectral_Dim,1);
                for i = 1:p.Spectral_Dim
                    if ~isempty(RF_Ident.SC_Box_Num_RF_pixels{i})
                        loop_minSC(i) = min(RF_Ident.SC_Box_RF_Pixel_covar{i});
                    end
                end
                [~,minSC_Box_index] = min(loop_minSC); % colour with highest SC
                minSC_Box_index     = minSC_Box_index(1);  % if there are multiple colours with co-highest, take the first
            end
        end
        
        % ASP RF Temporal Stixel Covar (FullRF)
        if p.RF_Type(2) == 1     % All Significant Pixels
            if (p.Time_Plot_Choice == 2) && (RF_Ident.SC_ASP_FullRF_Num_pixels > 0) % (Mod 18,12,2020, was: '...&& (~isempty(RF_Ident.SC_ASP_FullRF_Num_pixels))')
                loop_minSC = NaN(p.Spectral_Dim,1);
                for i = 1:p.Spectral_Dim
                    if ~isempty(RF_Ident.SC_ASP_Num_RF_pixels{i})
                        loop_minSC(i) = min(RF_Ident.SC_ASP_RF_Pixel_covar{i});
                    end
                end
                [~,minSC_ASP_index] = min(loop_minSC); % colour with highest SC
                minSC_ASP_index     = minSC_ASP_index(1);  % if there are multiple colours with co-highest, take the first
            end
        end
        % Don't need below as use ASP coords
        %     % Gaus RF Temporal Stixel Covar (FullRF)
        %     if p.RF_Type(3) == 1     % Gaussian
        %         if (p.Time_Plot_Choice == 2) && (~isempty(RF_Ident.SC_Gaus_FullRF_Num_pixels))
        %             loop_minSC = NaN(p.Spectral_Dim,1);
        %             for i = 1:p.Spectral_Dim
        %                 if ~isempty(RF_Ident.SC_Gaus_Num_RF_pixels{i})
        %                     loop_minSC(i) = min(RF_Ident.SC_Gaus_RF_Pixel_covar{i});
        %                 end
        %             end
        %             [~,minSC_index] = min(loop_minSC); % colour with highest SC
        %             minSC_index     = minSC_index(1);  % if there are multiple colours with co-highest, take the first
        %         end
        %     end
    end
    
    
    % Find STA, STASD, Temporal STA, Stixel, Pixel and Temporal Stixel Covariances (Full RF)
    for i = 1:p.Spectral_Dim
        
        if     p.RF_Ident_Meth_vec(1) == 1 % STA-SD
            
            % Full RF, Box RF STA, STASD and Temp. STA
            if     p.RF_Type(1) == 1 % Box
                if RF_Ident.STASD_Box_FullRF_Num_pixels > 0 % (Mod 18,12,2020, was: '~isempty(RF_Ident.STASD_Box_FullRF_Num_pixels)')
                    RF_Ident.STASD_Box_FullRF_STA{i}    = NaN(RF_Ident.STASD_Box_FullRF_Num_pixels,p.Num_STE_bins);
                    RF_Ident.STASD_Box_FullRF_STA_SD{i} = NaN(RF_Ident.STASD_Box_FullRF_Num_pixels,1);
                    for j = 1:RF_Ident.STASD_Box_FullRF_Num_pixels
                        RF_Ident.STASD_Box_FullRF_STA{i}(j,:)  =    STA(RF_Ident.STASD_Box_FullRF_coords(j,1),RF_Ident.STASD_Box_FullRF_coords(j,2),:,i);
                        RF_Ident.STASD_Box_FullRF_STA_SD{i}(j) = STA_SD(RF_Ident.STASD_Box_FullRF_coords(j,1),RF_Ident.STASD_Box_FullRF_coords(j,2),i);
                    end
                    
                    if p.Time_Plot_Choice == 1     % mean
                        RF_Ident.STASD_Box_FullRF_STA_time{i} = NaN(p.Num_STE_bins,1);
                        for j = 1:p.Num_STE_bins
                            RF_Ident.STASD_Box_FullRF_STA_time{i}(j) = mean(RF_Ident.STASD_Box_FullRF_STA{i}(:,j));
                        end
                    else % p.Time_Plot_Choice == 2 % most significant pixel
                        RF_Ident.STASD_Box_FullRF_STA_time{i} = squeeze(STA(RF_Ident.STASD_Box_RF_coords_centre{maxSTASD_Box_index}(1),RF_Ident.STASD_Box_RF_coords_centre{maxSTASD_Box_index}(2),:,i));
                    end
                end
            end
            
            % Full RF, ASP RF STA, STASD and Temp. STA
            if p.RF_Type(2) == 1     % All Significant Pixels
                if RF_Ident.STASD_ASP_FullRF_Num_pixels > 0 % (Mod 18,12,2020, was: '~isempty(RF_Ident.STASD_ASP_FullRF_Num_pixels)')
                    RF_Ident.STASD_ASP_FullRF_STA{i}    = NaN(RF_Ident.STASD_ASP_FullRF_Num_pixels,p.Num_STE_bins);
                    RF_Ident.STASD_ASP_FullRF_STA_SD{i} = NaN(RF_Ident.STASD_ASP_FullRF_Num_pixels,1);
                    for j = 1:RF_Ident.STASD_ASP_FullRF_Num_pixels
                        RF_Ident.STASD_ASP_FullRF_STA{i}(j,:)  =    STA(RF_Ident.STASD_ASP_FullRF_coords(j,1),RF_Ident.STASD_ASP_FullRF_coords(j,2),:,i);
                        RF_Ident.STASD_ASP_FullRF_STA_SD{i}(j) = STA_SD(RF_Ident.STASD_ASP_FullRF_coords(j,1),RF_Ident.STASD_ASP_FullRF_coords(j,2),i);
                    end
                    
                    if p.Time_Plot_Choice == 1     % mean
                        RF_Ident.STASD_ASP_FullRF_STA_time{i} = NaN(p.Num_STE_bins,1);
                        for j = 1:p.Num_STE_bins
                            RF_Ident.STASD_ASP_FullRF_STA_time{i}(j) = mean(RF_Ident.STASD_ASP_FullRF_STA{i}(:,j));
                        end
                    else % p.Time_Plot_Choice == 2 % most significant pixel
                        RF_Ident.STASD_ASP_FullRF_STA_time{i} = squeeze(STA(RF_Ident.STASD_ASP_RF_coords{maxSTASD_ASP_index}(1,1),RF_Ident.STASD_ASP_RF_coords{maxSTASD_ASP_index}(1,2),:,i));
                    end
                end
            end
            
            % Full RF, Gaus RF STA, STASD and Temp. STA
            if p.RF_Type(3) == 1     % Gaussian
                if RF_Ident.STASD_Gaus_FullRF_Num_pixels > 0 % (Mod 18,12,2020, was: '~isempty(RF_Ident.STASD_Gaus_FullRF_Num_pixels)')
                    RF_Ident.STASD_Gaus_FullRF_STA{i}    = NaN(RF_Ident.STASD_Gaus_FullRF_Num_pixels,p.Num_STE_bins);
                    RF_Ident.STASD_Gaus_FullRF_STA_SD{i} = NaN(RF_Ident.STASD_Gaus_FullRF_Num_pixels,1);
                    for j = 1:RF_Ident.STASD_Gaus_FullRF_Num_pixels
                        RF_Ident.STASD_Gaus_FullRF_STA{i}(j,:)  =    STA(RF_Ident.STASD_Gaus_FullRF_coords(j,1),RF_Ident.STASD_Gaus_FullRF_coords(j,2),:,i);
                        RF_Ident.STASD_Gaus_FullRF_STA_SD{i}(j) = STA_SD(RF_Ident.STASD_Gaus_FullRF_coords(j,1),RF_Ident.STASD_Gaus_FullRF_coords(j,2),i);
                    end
                    
                    if p.Time_Plot_Choice == 1     % mean
                        RF_Ident.STASD_Gaus_FullRF_STA_time{i} = NaN(p.Num_STE_bins,1);
                        for j = 1:p.Num_STE_bins
                            RF_Ident.STASD_Gaus_FullRF_STA_time{i}(j) = mean(RF_Ident.STASD_Gaus_FullRF_STA{i}(:,j));
                        end
                    else % p.Time_Plot_Choice == 2 % most significant pixel
                        RF_Ident.STASD_Gaus_FullRF_STA_time{i} = squeeze(STA(RF_Ident.STASD_ASP_RF_coords{maxSTASD_ASP_index}(1,1),RF_Ident.STASD_ASP_RF_coords{maxSTASD_ASP_index}(1,2),:,i));
                    end
                end
            end
            
        end
        
        
        if     p.RF_Ident_Meth_vec(4) == 1 % SC
            
            % Full RF, Box RF Stix. Cov., Pix. Cov. and Temp. Stix. Cov.
            if     p.RF_Type(1) == 1 % Box
                if RF_Ident.SC_Box_FullRF_Num_pixels > 0 % (Mod 18,12,2020, was: '~isempty(RF_Ident.SC_Box_FullRF_Num_pixels)')
                    RF_Ident.SC_Box_FullRF_Stixel_covar{i}    = NaN(RF_Ident.SC_Box_FullRF_Num_pixels,p.Num_STE_bins);
                    RF_Ident.SC_Box_FullRF_Pixel_covar{i} = NaN(RF_Ident.SC_Box_FullRF_Num_pixels,1);
                    for j = 1:RF_Ident.SC_Box_FullRF_Num_pixels
                        RF_Ident.SC_Box_FullRF_Stixel_covar{i}(j,:) = SCA_Stixel_covar(RF_Ident.SC_Box_FullRF_coords(j,1),RF_Ident.SC_Box_FullRF_coords(j,2),:,i);
                        RF_Ident.SC_Box_FullRF_Pixel_covar{i}(j)    =  SCA_Pixel_covar(RF_Ident.SC_Box_FullRF_coords(j,1),RF_Ident.SC_Box_FullRF_coords(j,2),i);
                    end
                    
                    if p.Time_Plot_Choice == 1     % mean
                        RF_Ident.SC_Box_FullRF_Stixel_covar_time{i} = NaN(p.Num_STE_bins,1);
                        for j = 1:p.Num_STE_bins
                            RF_Ident.SC_Box_FullRF_Stixel_covar_time{i}(j) = mean(RF_Ident.SC_Box_FullRF_Stixel_covar{i}(:,j));
                        end
                    else % p.Time_Plot_Choice == 2 % most significant pixel
                        RF_Ident.SC_Box_FullRF_Stixel_covar_time{i} = squeeze(SCA_Stixel_covar(RF_Ident.SC_Box_RF_coords_centre{minSC_Box_index}(1),RF_Ident.SC_Box_RF_coords_centre{minSC_Box_index}(2),:,i));
                    end
                end
            end
            
            % Full RF, ASP RF Stix. Cov., Pix. Cov. and Temp. Stix. Cov.
            if p.RF_Type(2) == 1     % All Significant Pixels
                if RF_Ident.SC_ASP_FullRF_Num_pixels > 0 % (Mod 18,12,2020, was: '~isempty(RF_Ident.SC_ASP_FullRF_Num_pixels)')
                    RF_Ident.SC_ASP_FullRF_Stixel_covar{i} = NaN(RF_Ident.SC_ASP_FullRF_Num_pixels,p.Num_STE_bins);
                    RF_Ident.SC_ASP_FullRF_Pixel_covar{i}  = NaN(RF_Ident.SC_ASP_FullRF_Num_pixels,1);
                    for j = 1:RF_Ident.SC_ASP_FullRF_Num_pixels
                        RF_Ident.SC_ASP_FullRF_Stixel_covar{i}(j,:) = SCA_Stixel_covar(RF_Ident.SC_ASP_FullRF_coords(j,1),RF_Ident.SC_ASP_FullRF_coords(j,2),:,i);
                        RF_Ident.SC_ASP_FullRF_Pixel_covar{i}(j)    =  SCA_Pixel_covar(RF_Ident.SC_ASP_FullRF_coords(j,1),RF_Ident.SC_ASP_FullRF_coords(j,2),i);
                    end
                    
                    if p.Time_Plot_Choice == 1     % mean
                        RF_Ident.SC_ASP_FullRF_Stixel_covar_time{i} = NaN(p.Num_STE_bins,1);
                        for j = 1:p.Num_STE_bins
                            RF_Ident.SC_ASP_FullRF_Stixel_covar_time{i}(j) = mean(RF_Ident.SC_ASP_FullRF_Stixel_covar{i}(:,j));
                        end
                    else % p.Time_Plot_Choice == 2 % most significant pixel
                        RF_Ident.SC_ASP_FullRF_Stixel_covar_time{i} = squeeze(SCA_Stixel_covar(RF_Ident.SC_ASP_RF_coords{minSC_ASP_index}(1,1),RF_Ident.SC_ASP_RF_coords{minSC_ASP_index}(1,2),:,i));
                    end
                end
            end
            
            % Full RF, Gaus RF Stix. Cov., Pix. Cov. and Temp. Stix. Cov.
            if p.RF_Type(3) == 1     % Gaussian
                if RF_Ident.SC_Gaus_FullRF_Num_pixels > 0 % (Mod 18,12,2020, was: '~isempty(RF_Ident.SC_Gaus_FullRF_Num_pixels)')
                    RF_Ident.SC_Gaus_FullRF_Stixel_covar{i} = NaN(RF_Ident.SC_Gaus_FullRF_Num_pixels,p.Num_STE_bins);
                    RF_Ident.SC_Gaus_FullRF_Pixel_covar{i}  = NaN(RF_Ident.SC_Gaus_FullRF_Num_pixels,1);
                    for j = 1:RF_Ident.SC_Gaus_FullRF_Num_pixels
                        RF_Ident.SC_Gaus_FullRF_Stixel_covar{i}(j,:) = SCA_Stixel_covar(RF_Ident.SC_Gaus_FullRF_coords(j,1),RF_Ident.SC_Gaus_FullRF_coords(j,2),:,i);
                        RF_Ident.SC_Gaus_FullRF_Pixel_covar{i}(j)    =  SCA_Pixel_covar(RF_Ident.SC_Gaus_FullRF_coords(j,1),RF_Ident.SC_Gaus_FullRF_coords(j,2),i);
                    end
                    
                    if p.Time_Plot_Choice == 1     % mean
                        RF_Ident.SC_Gaus_FullRF_Stixel_covar_time{i} = NaN(p.Num_STE_bins,1);
                        for j = 1:p.Num_STE_bins
                            RF_Ident.SC_Gaus_FullRF_Stixel_covar_time{i}(j) = mean(RF_Ident.SC_Gaus_FullRF_Stixel_covar{i}(:,j));
                        end
                    else % p.Time_Plot_Choice == 2 % most significant pixel
                        RF_Ident.SC_Gaus_FullRF_Stixel_covar_time{i} = squeeze(SCA_Stixel_covar(RF_Ident.SC_ASP_RF_coords{minSC_ASP_index}(1,1),RF_Ident.SC_ASP_RF_coords{minSC_ASP_index}(1,2),:,i));
                    end
                end
            end
            
        end
        
    end
    
    
end

