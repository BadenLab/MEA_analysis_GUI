function out = findfile_app (stim_idx,varargin)
%This function loads a given file from the stimulus analysis files. It uses
%the information stored in the Stimulus Data file.
%@MSeifert 2020

%% Input control
p = inputParser;
validScalarPosNum = @(x) isnumeric(x) && isscalar(x);
ischarstring = @(x) ischar(x) || isstring(x);
%validMatrixShape = @(x) isnumeric(x) && (size(x,1) ==1) && (size(x,2) == 2);
addRequired(p,'stim_idx',validScalarPosNum);
addRequired(p,'savepath',ischarstring);
addRequired(p,'file',ischarstring);
addParameter(p,'subfolder',ischarstring);


parse(p,stim_idx,varargin{:});

fields = fieldnames(p.Results);
variables = p.Results;

for ii = 1:length(fields)
   if isa(variables.(fields{ii}),'function_handle')
       variables = rmfield(variables,fields{ii});
   end
       
end

%Check and correct file formate
if isstring(variables.file)
    variables.file = char(variables.file);
end

%Check if the file has .mat at the end

try
if ~strcmp(variables.file(end-3:end),'.mat')
    if strcmp(variables.file(end-3:end),'.fig')
        warning("Assuming .fig formate for file")
        
    elseif contains(variables.file,'.')
        error(['Unclear file formate, check file ending or file name.'...
            ' Make sure name doesnt contain "." except for the file ending'])
    else
        variables.file = [variables.file,'.mat'];
    end
end
    
    
catch %This is in case the file is shorted than 3 characters so the indexing above will fail
    if contains(variables.file,'.')
        error(['Unclear file formate, check file ending or file name.'...
            ' Make sure name doesnt contain "." except for the file ending'])
    else
        variables.file = [variables.file,'.mat'];
    end
end

%Load the information where the data file is stored
S = load(variables.savepath,'Stimulus_info');
stimulus_info = S.Stimulus_info{variables.stim_idx};

%Load the data file
S = load(stimulus_info);
Data = S.Data;

data_length = length(Data);
%Loop over the structure length (that is over the folders in the analysis
%folder) to find the variable that is looked for
if ~isfield(variables,'subfolder')
    for ii = 1:data_length

        Files = Data(ii).Files;
        IndexF = strfind(Files,variables.file);
        try
        right_folder = nnz(find(not(cellfun('isempty',IndexF))));
        catch
            right_folder = 0;
        end
        if right_folder == 1
            Index = ii;
        end
    end

    if Index == 1
        %Reconstruct filename if folder is "main"
        IndexS = max(strfind(stimulus_info,'\'));
        stim_folder = stimulus_info(1:IndexS);

        out = strcat(stim_folder,'\',variables.file);

    else

        %Reconstruct filename if folder is not "main"
        IndexS = max(strfind(stimulus_info,'\'));
        stim_folder = stimulus_info(1:IndexS);

        out = strcat(stim_folder,Data(Index).Folder,'\',variables.file);
    end

else
%     wanted_str = contains({Data(2:end).Folder},variables.subfolder);
      wanted_str = ismember({Data(2:end).Folder},variables.subfolder);
    if ~any(wanted_str)
        warning('No file found')
        out = [];
        return
    end
    wanted_idx = find(wanted_str)+1;
    
    %Files = Data(wanted_idx).Files;
    %IndexF = contains(Files,variables.file);
    
    %Reconstruct filename if folder is not "main"
    IndexS = max(strfind(stimulus_info,'\'));
    stim_folder = stimulus_info(1:IndexS);

    out = strcat(stim_folder,Data(wanted_idx).Folder,'\',variables.file);
    
end
    
    




end