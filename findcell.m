function out = findcell (savepath,cell_idx,varargin)


p = inputParser;
default_search_all_files = false;
validScalarPosNum = @(x) isnumeric(x) && isscalar(x) && (x > 0);
addRequired(p,'savepath',@ischar);
addRequired(p,'cell_idx',validScalarPosNum);
addOptional(p,'stim_idx',validScalarPosNum);
addParameter(p,'search_all_files',default_search_all_files,@islogical);
parse(p,savepath,cell_idx,varargin{:});
%% Set variables

%This will remove all variables from the input structure which are not
%beeing used
fields = fieldnames(p.Results);
variables = p.Results;

for ii = 1:length(fields)
   if isa(variables.(fields{ii}),'function_handle')
       variables = rmfield(variables,fields{ii});
   end
       
end

%If no specific stimulus ID is given, function will search all stimuli folders
M = matfile(savepath,'Writable',false); 
if ~isfield(variables,'stim_idx')
  
   variables.stim_idx = 1:1:length(M.stim_list);
end

%Load stimulus info file
stimulus_info = M.Stimulus_info;
stim_list = M.stim_list;
search_all_files = variables.search_all_files;
%Check if all files shall be searched
if search_all_files 
    disp('All files will be searched, this may take a while')

%We search for the information in the given stimuli

for ii = variables.stim_idx
    ii
    try
    if isempty(stimulus_info{ii})
        continue %If the stimulus Data file has not been created yet
    end
    catch
        continue
    end
   
    S = load(stimulus_info{ii});
    Data = S.Data;
       
    %Check how many folders need to be searched
    nr_folders = numel([Data.Folder]);
    
    for ff = 1:nr_folders
        nr_files = numel(Data(ff).Files);
        Files = Data(ff).Files;
        for fi = 1:nr_files
            
            if contains(Files{fi},'.mat')
                M = matfile(findfile_app(ii,savepath,Files{fi}),'Writable',false);
                %Check how many variables are in the file
                M_fields = who(M);
                
                for mf = 1:length(M_fields)
                    if isfield(M.(M_fields{mf}),'cell_idx')
                        %Now we load that variable
                        variable = M.(M_fields{mf});
                        %Search for the
                        try
                            location = find([variable.cell_idx]==cell_idx);
                            variable = variable(location);
                        catch
                            location = find([variable.cell_idx]==num2str(cell_idx));
                            variable = variable(location);
                        end
                    else 
                        continue
                    end
                    variable_out.(['Stim',num2str(ii)]).(M_fields{mf}) = variable;
                    
                end
                
            end
                
            
        end
    end
    
    
    
    
    
    
end
    
elseif ~search_all_files
    
    a = 0;
    for ii = variables.stim_idx
    ii
    a = a+1;
    try
    if isempty(stimulus_info{ii})
        continue %If the stimulus Data file has not been created yet
    end
    catch
        continue
    end
   
    S = load(stimulus_info{ii});
    Data = S.Data;
    
    main_idx = find([Data.Folder] == 'Main'); %Should be first but you never know
        
    nr_files = numel(Data(main_idx).Files);
    Files = Data(main_idx).Files;
    for fi = 1:nr_files
        fi
        if contains(Files{fi},'.mat')
            M = matfile(findfile_app(ii,savepath,Files{fi}),'Writable',false);
            %Check how many variables are in the file
            M_fields = who(M);

            for mf = 1:length(M_fields)
                if isfield(M.(M_fields{mf}),'cell_idx')
                    %Now we load that variable
                    variable = M.(M_fields{mf});
                    %Search for the
                    try
                        location = find([variable.cell_idx]==cell_idx);
                        variable = variable(location);
                    catch
                        location = find([variable.cell_idx]==num2str(cell_idx));
                        variable = variable(location);
                    end
                else 
                    continue
                end
                if isempty(location)
                    continue
                end
                
                    
                variable_out.(['Stim',num2str(ii)]).(M_fields{mf}) = variable;
                variable_out.(['Stim',num2str(ii)]).stimname...
                    = stim_list{variables.stim_idx(a)};

            end

        end


    end
    end
    
    
    
    
    
    
end
    out = variable_out;
    
end