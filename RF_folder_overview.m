function out = RF_folder_overview (savepath,add_info)

M = matfile(savepath);

    
stimulus_location = M.Stimulus_info(1,add_info.stim_idx);
Data = load(stimulus_location{:});
Data = {Data.Data.Folder};
%Remove main from folder list as we are only interested in subfolders
Data = Data(2:end);


add_info.listbox.Savedanalysisfiles.Items = Data;


out = 1;

    
    
    
end