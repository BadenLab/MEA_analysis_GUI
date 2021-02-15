function out = sf_organizer (stim_idx,savepath,varargin)
%% About
%This function is the main save function for all data in the app. It creates
%a new folder for a given stimulus if a folder does not yet exist or uses
%an existing folder in case a folder for that stimulus does exist. It save
%variables to given folder or subfolder and keeps track of all variables
%and subfolders created by creating a stimulus Data file which will be
%refered to in the Stimulus_Info file in the main analysis .mat file.

%The return of this function is the path of the folder, subfolder, or in
%case a variable was saved the whole path to the variable. 
%@MSeifert 2020

%% Input control
defaultOverwrite = false;
default_update_data = true;
defaultCollect = false;
p = inputParser;
validScalarPosNum = @(x) isnumeric(x) && isscalar(x);
%validMatrixShape = @(x) isnumeric(x) && (size(x,1) ==1) && (size(x,2) == 2);
addRequired(p,'stim_idx',validScalarPosNum);
addRequired(p,'savepath',@ischar);
addParameter(p,'subfoldername',@ischar);
addParameter(p,'pathname',@ischar);
addParameter(p,'variable_name',@isstring);
addParameter(p,'variable',@(x)validateattributes(x,{'nonempty'}));
addParameter(p,'filename',@isstring);
addParameter(p,'overwrite',defaultOverwrite,@islogical);
addParameter(p,'update_data',default_update_data,@islogical);
addParameter(p,'collect_files',defaultCollect,@islogical);

parse(p,stim_idx,savepath,varargin{:});

%Just for debugging
% out = p;
%% Set variables

%This will remove all variables from the input structure which are not
%beeing used
fields = fieldnames(p.Results);
variables = p.Results;
Data_idx = 1;


for ii = 1:length(fields)
   if isa(variables.(fields{ii}),'function_handle')
       variables = rmfield(variables,fields{ii});
   end
       
end

out = variables;
% Add some global variables
stim_folder_name = ['Stimulus_',num2str(stim_idx)];
datafile_name = ['Stimulus_',num2str(variables.stim_idx),'_DataFile.mat'];



%Check if this is run in parallel
%Parallel access to a matfile is not possible

% answer = is_in_parallel();
% 
% if answer == 1 && update_data == 1
%     error("Cant update stimulus data in parallel pool")
% end


%% Check pathname
%If no pathname is given, current folder is set to be the savepath

if ~isfield(variables,'pathname')
    try
    
        S = load(savepath,'pathname');
        pathname = S.pathname;
        clear S
    catch
        error('No pathname found in savefile, savefile incorrect');
        out = 0;
    end
    
else
    %In case pathname is given as input
    pathname = variables.pathname;
end
cd(pathname)

%% Check if stimulus folder exists in pathname
%First we check if the pathname is the stimulus folder (might be a mistake
%people can do

if strcmp(pathname,stim_folder_name)
    out = 0;
    error('Pathname set to stimulus folder name, check pathname input')
end

%Next, check if the folder does already exist
c_data_path = strcat(pathname, stim_folder_name,"\",datafile_name);

if isfolder(stim_folder_name)
    %disp('Stimulus folder does exist, no new folder created');
    cd(stim_folder_name)
    if variables.update_data == 1
    M = matfile(c_data_path,'Writable',true);
    Data = M.Data;
    end
       
    out = pwd;
else
    %Otherwise the folder will be created newly
    mkdir(stim_folder_name);
    cd(stim_folder_name);
    
    %If the folder didnt exist before, than also a new data file has to be
    %created to keep track of the subfolders and subfiles
    
    if variables.update_data == 1
    M = matfile(c_data_path,'Writable',true);
    end
    
    
    %Write the location of the new create matfile into the main analysis
    %file of the app
    try 
        S = load(savepath,'Stimulus_info');
        Stimulus_info = S.Stimulus_info;
    catch
        Stimulus_info = {};
    end
    Stimulus_info{stim_idx} = [pathname,stim_folder_name,'\',datafile_name];
    S = matfile(savepath,'Writable',true);
    S.Stimulus_info = Stimulus_info;
    
    
    
    %Create empty datafile structure
    Data(1).Folder = "Main";
    Data(1).Files = {};
    M.Data = Data;
    Data_idx = 1;
    out = pwd;
    
end

%% Create subfolder or switch folder
if isfield(variables, 'subfoldername')
    
    %Check if subfolder exists in stimulus folder
    if isfolder(variables.subfoldername)
        cd(variables.subfoldername)
%         variables.subfoldername
        if variables.update_data == 1
            Data_idx = find([Data.Folder] == variables.subfoldername);
        end
        sf = 1; %This keeps the information that we are now in the subfolder dir
        out = pwd;
    else
        mkdir(variables.subfoldername)
        cd(variables.subfoldername)
        if variables.update_data == 1
            Data_length = length(Data);
            Data(Data_length+1).Folder = variables.subfoldername;
            Data_idx = Data_length+1;
        end
        out = pwd;
        sf = 1;
    end
    
end


%% Save variable
%This only applies if a variable and a name for that variable is given in
%the inputs otherwise skip this part

%Check if variable is figure, which requires a different way to save it


if isfield(variables,'variable_name')&&isfield(variables,'variable')
    try
        figure_save = isfigure(variables.variable);
    catch ME
        figure_save = false;
        warning([ME.message, ' Ignore, if no error thrown.'])
    end
     
    if figure_save
        file_ending = '.fig';
    else
        file_ending = '.mat';
    end
    
    %Creat new sturcture to dynamically save the variable
    save_struc.(variables.variable_name) = variables.variable;
       
    %Check if the a filename was given
    if ~isfield(variables,'filename')
        filename = variables.variable_name;
    else
        filename = variables.filename;
    end
       
    %First check if the variable already exists
      
    if isfile([filename,file_ending])&& variables.overwrite == false
        
        % Include the desired Default answer
        opts.Interpreter = 'tex';
        opts.Default = 'No';
        quest = ['Do you want to overwrite the existing file: "', filename,file_ending];
        answer = questdlg(quest,'Boundary Condition','Yes','No',opts);
                      
        if strcmp(answer,'No')
            return
        else
            if figure_save
                savefig(variables.variable,filename);
            else
                save(filename,'-struct','save_struc','-v7.3');
            end
            if variables.update_data == 1
                Data(Data_idx).Files{end+1} = [filename,file_ending];    
            end
            out = [pwd,'\',filename,file_ending];
            
        end
        
        
    else
        if figure_save
             savefig(variables.variable,filename);
        else
             save(filename,'-struct','save_struc');
        end
        
        if variables.update_data == 1
            Data(Data_idx).Files{end+1} = [char(filename),file_ending]; 
        end
        out = [pwd,'\',filename,file_ending];
    end
end

            
        
    

if variables.update_data == 1
    for i = 1:length(Data)
        if numel(Data(i).Files) > 1 %Because unique doesnt work (and isnt required)
            %if there is only one file in the datafile
            
            Data(i).Files = unique(Data(i).Files);
        end
            
    end
    M.Data = Data;
end


if variables.collect_files == 1 && variables.update_data == 1
    
    M = matfile([pathname,'\',stim_folder_name,'\',datafile_name],'Writable',true);
    Data = M.Data;
     files = dir;
     files = {(files(3:end,:).name)};
     %Sort files depending on ending numbers
     number_end = regexp(files,'\d*','Match');
     try
     if ~any(cellfun(@isempty,number_end))
         %Only sort if all files contain numbers, otherwise dont sort
        [~,Isort] = sort(str2double([number_end{:}]));
        files = files(Isort);
     end
     
         
     catch
         
     end
    Data(Data_idx).Files = files;
    M.Data = Data;
end


out = string(out);

end