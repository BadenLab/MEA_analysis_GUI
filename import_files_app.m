function out = import_files_app
%% This function loads a analysis .mat file and saves it so the analysis 
%app can access it

[file, path] = uigetfile("*.mat");

complete_name = strcat(path,file);

selpath = uigetdir([],"Select save location");


loaded_files = load(complete_name);

%Change filenames in the loaded files
loaded_files.pathname = [selpath,'\'];
loaded_files.savename = [selpath,'\',file];

%% Check which stimulus information exists

if isfield(loaded_files,"Stimulus_info")
    
    info_length = length(loaded_files.Stimulus_info);
    
    for ii = 1:info_length
       location = loaded_files.Stimulus_info(1,ii); 
       location = location{1,1};
       slashes =  strfind(loaded_files.Stimulus_info{1,ii},"\");
       cut_position = NaN(1,2);
       try
       cut_position(1) = second_max(slashes);
       cut_position(2) = max(slashes);
              
       %Check if folder exists in the selected path
       old_position = [path(1:end-1),location(cut_position(1):cut_position(2))];
       if exist(old_position,'dir') == 7
           
            new_position = [selpath,location(cut_position(1):cut_position(2))];
            copyfile(old_position,new_position)
            
            new_info =[selpath,location(cut_position(1):end)];
            loaded_files.Stimulus_info{1,ii} = new_info;
       else
            loaded_files.Stimulus_info{1,ii} = [];
       end
       
       catch
           if isempty(second_max(slashes))
               loaded_files.Stimulus_info{1,ii} = [];
           else
               error('Problem with Stimulus_info variable, check in file')
           end
       end
       
        
     end
       
       
    
end


%% Save updated analysis file
new_file = strcat(selpath,"\",file);

save(new_file,'-struct',"loaded_files")

out = loaded_files;


end