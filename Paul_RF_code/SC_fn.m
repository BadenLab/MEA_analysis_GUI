function [Stixel_covar,Pixel_covar] = SC_fn(STE_Full,p)
% Calculate the SC array

%STE_Full = NaN(p.stim_rows,p.stim_columns,p.Num_STE_bins,p.length_spike_times);

Stixel_covar = NaN(p.stim_rows,p.stim_columns,p.Num_STE_bins);
%Pixel_covar = NaN(p.stim_rows,p.stim_columns);

for i = 1:p.Num_STE_bins
    
    for j = 1:p.stim_columns
        
        for k = 1:p.stim_rows
            
            Stixel_covar(k,j,i) = var(squeeze(STE_Full(k,j,i,:))); % or cov
            
        end
        
    end
    
end

Pixel_covar = min(Stixel_covar,[],3); % min(Stixel_covar,[],3) mean(Stixel_covar,3) max(Stixel_covar,[],3)