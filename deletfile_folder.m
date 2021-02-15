function folder_content = delete_dataset (savepath, add_info)
% This function lists all folders and files in a specific stimulus analysis
% folder and puts them into a table which is displayed on the gui to see
% which data has been produces so far for each stimulus.


M = matfile(savepath,'Writable',false);
pathname = M.pathname;
%Check whats in the folder
folder_content = dir(strcat(pathname,"Stimulus_",num2str(add_info.stim_idx)));
if size(folder_content,1) == 0 %If folder is empty return empty table
    folder_content = array2table(zeros(0,3));
    folder_content.Properties.VariableNames = {'Name','Date','Folder'};
    return
end

%if not empty, make cell with contents

folder_content = rmfield(folder_content,'datenum');
folder_content = rmfield(folder_content,'bytes');
folder_content = rmfield(folder_content,'folder');

%Remove first two entries (meaningless in the context)
folder_content(1:2) = [];


%Remove the time information from data field (because its unnecessary info)

for ii = 1:length(folder_content)
   date_before = folder_content(ii).date;
   split_idx = strfind(date_before,' ');
   folder_content(ii).date = date_before(1:split_idx-1);
end
    
%make a table
folder_content = struct2table(folder_content);
%Give the right names for the columns
folder_content.Properties.VariableNames = {'Name','Date','Folder'};

end