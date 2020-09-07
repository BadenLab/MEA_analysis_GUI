function STA = STA_fn_v2(STE_Full,mean_raw_stim_arr,p)
% Calculate the STA

%STE_Full = NaN(p.stim_rows,p.stim_columns,p.Num_STE_bins,p.length_spike_times);
% Time choice doesn't matter, fn is same either way

if     p.STA_Choice  == 1 % don't subtract average stim;
    STA = mean(STE_Full,4);
elseif p.STA_Choice  == 2 % subtract average stim;
    STA = mean(STE_Full,4) - mean_raw_stim_arr;
end