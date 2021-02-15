function out = delete_dataset (savepath, add_info)
%This function deletes a file or folder from a stimulus analysis file and
%also removes its entry in the DataFile

%% Failsave
%just in case delete was pressed by accident, here is the last chance to
%say "no" ;-)
 opts.Interpreter = 'tex';
        opts.Default = 'No';
        quest = ['Do you really delete the file:', add_info.analysis_file_select{1,1}];
        answer = questdlg(quest,'Boundary Condition','Yes','No',opts);
        
if ~strcmp(answer,'Yes')
    out = 1;
    return
end
%% Remove folder/file
M = matfile(savepath,'Writable',false);
pathname = M.pathname;
%Check whats in the folder
folder_content = dir(strcat(pathname,"Stimulus_",num2str(add_info.stim_idx)));
for ii = 1:length(folder_content)
    if strcmp(folder_content(ii).name,add_info.analysis_file_select{1,1})
        if ~folder_content(ii).isdir   %if its a file
            delete(strcat(folder_content(ii).folder,"\",folder_content(ii).name))
        else %if its a folder
            [status,msg] = rmdir(strcat(folder_content(ii).folder,"\",folder_content(ii).name),'s');
        end
    end
end

%% Update Data File
warning(msg);
warning('Folder is probably open in another program');
if status %Only delete from Data File if delete has worked
M = matfile(savepath,'Writable',false);
Stimulus_info = M.Stimulus_info;
L = matfile(Stimulus_info{1,add_info.stim_idx},'Writable',true);
Data = L.Data;
clear M

for ii = 1:length(Data)
    if strcmp(add_info.analysis_file_select{1,1},Data(ii).Folder)
        Data(ii) = [];
        break
    end
end

L.Data = Data;
out = 1;
end
end
        







