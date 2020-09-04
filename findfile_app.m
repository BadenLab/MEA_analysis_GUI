function out = findfile_app (stim_idx,savepath,file)
%This function loads a given file from the stimulus analysis files. It uses
%the information stored in the Stimulus Data file.
%@MSeifert 2020


%Check and correct file formate
if isstring(file)
    file = char(file);
end

%Check if the file has .mat at the end

try
if ~strcmp(file(end-3:end),'.mat')
    if strcmp(file(end-3:end),'.fig')
        warning("Assuming .fig formate for file")
        
    elseif contains(file,'.')
        error(['Unclear file formate, check file ending or file name.'...
            ' Make sure name doesnt contain "." except for the file ending'])
    else
        file = [file,'.mat'];
    end
end
    
    
catch %This is in case the file is shorted than 3 characters so the indexing above will fail
    if contains(file,'.')
        error(['Unclear file formate, check file ending or file name.'...
            ' Make sure name doesnt contain "." except for the file ending'])
    else
        file = [file,'.mat'];
    end
end

    


%Load the information where the data file is stored
S = load(savepath,'Stimulus_info');
stimulus_info = S.Stimulus_info{stim_idx};

%Load the data file
S = load(stimulus_info);
Data = S.Data;

data_length = length(Data);
%Loop over the structure length (that is over the folders in the analysis
%folder) to find the variable that is looked for

for ii = 1:data_length
    
    Files = Data(ii).Files;
    IndexF = strfind(Files,file);
    
    right_folder = nnz(find(not(cellfun('isempty',IndexF))));
    if right_folder == 1
        Index = ii;
    end
end

if Index == 1
    %Reconstruct filename if folder is "main"
    IndexS = max(strfind(stimulus_info,'\'));
    stim_folder = stimulus_info(1:IndexS);

    out = strcat(stim_folder,'\',file);

else

    %Reconstruct filename if folder is not "main"
    IndexS = max(strfind(stimulus_info,'\'));
    stim_folder = stimulus_info(1:IndexS);

    out = strcat(stim_folder,Data(Index).Folder,'\',file);
end

    
    
    
    
    




end