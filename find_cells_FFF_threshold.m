function Items_list = find_cells_FFF_threshold(savepath,add_info,threshold)
%% Find cells according to set threshold function
% This function find all cells that pass the set threshold and returns a
% list of their names to the gui. 
%@mars 2021

%% Get the data
S = load(findfile_app(add_info.stim_idx,savepath,'FFF_average'));
FFF_average = S.FFF_average;
clear S
%Check which cells pass the quality index
QC_pass = [FFF_average.stats_max] >= threshold;

if ~any(QC_pass)
    Items_list = {};
    warning('No cells passed the set threshold')
    return
end

cell_name = {FFF_average.cell_idx};

cell_name_pass = cell_name(QC_pass);

%Reconstruct the real displayed name
for ii = 1:length(cell_name_pass)
    
   Items_list{1,ii} = strcat('Cell_',num2str(cell_name_pass{1,ii})); 
    
end





end