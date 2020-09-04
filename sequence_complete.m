
%% This function takes a stimulus sequence and reshapes it in a way that
%%every entry of the sequence represents 1ms of the stimulus
%%Inputs: sequence: The original sequence;
%%up_factor: The factor by which the sequence has to be multiplied to get a
%%1ms one entry sequence;
%%triggger: the trigger channel as array with the trigger events
function [sequence_out, diff_trigger] = sequence_complete (sequence, trigger)
ms = 0.001;
trigger_length = length(trigger);
diff_trigger = diff(trigger);
diff_trigger_idx = round(diff_trigger / ms);
sum_trigger = sum(diff_trigger_idx);




%Check the size of the sequence
sequence_size = size(sequence);
nr_colours = sequence_size(3);
%Test if the sequence fits to the number of trigger events
if trigger_length ~= sequence_size(1)+1
    disp('trigger channel does not fit to the size of the sequence, function return is empty')
    return
end
sequence_out = single(zeros(sum_trigger,sequence_size(2),sequence_size(3)));

% variable for the index in the new matrix;
idx = 1;
for ii = 1:sequence_size(1)
    idx_e = idx + diff_trigger_idx(ii)-1;
     for cc = 1:nr_colours
         sequence_out(idx:idx_e,:,cc) = repelem(squeeze(sequence(ii,:,cc)),diff_trigger_idx(ii),1);
     end

%      for cc = 1:nr_colours
%          sequence_out(idx:idx_e,:,cc) = repelem_log(squeeze(sequence(ii,:,cc)),diff_trigger_idx(ii));
%      end


    
%      for cc = 1:nr_colours
%          B = ones(diff_trigger_idx(ii),1);
%          sequence_out(idx:idx_e,:,cc) = kron(squeeze(sequence(ii,:,cc)),B);
%      end
    %sequence_out(idx:idx_e,:,:) = repelem(sequence(ii,:,:),diff_trigger_idx(ii),1,1);
    %repelem(sequence(ii,:,:),diff_trigger_idx(ii),1,1);
    if ii<trigger_length
    idx = idx_e+1;
    end
end
    
  
            
        
        
        
        
    
end