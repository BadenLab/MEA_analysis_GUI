function out = filter_RF_cells (savepath,add_info)

% THis function updates the RF table showing either all cells or just the
% onces that have a receptive field

%% If only RF cells shall be shown
if add_info.settings.kernel_new.filter_CellswithRF
    Data = add_info.tables.RF_single_cell.Data;
    %Check which cells have a RF
    RF_cell_idx = any(table2array(Data(:,2:end)),2);
    Data = (Data(RF_cell_idx,:));
    add_info.tables.RF_single_cell.Data = Data;
        
else
    table_heads = {'Cell idx'};
    M = matfile(savepath,'Writable',false);
    RF_overview_log = false(length(M.cell_indices),length(add_info.settings.kernel_new.folderplot));
    for ii = 1:length(add_info.settings.kernel_new.folderplot)
        legend_name{ii} = strrep(add_info.settings.kernel_new.folderplot{:,ii},'_',' ');
        Stimulus_info = M.Stimulus_info(1,add_info.stim_idx);
        Stimulus_info = Stimulus_info{:};
        cut_path = strfind(Stimulus_info,'\');
        path_to_data = [Stimulus_info(1:cut_path(end)),add_info.settings.kernel_new.folderplot{:,ii}];
        S = load([path_to_data,'\RF_overview.mat']);
        RF_overview = S.RF_overview;
        clear S

        %Split off cell idx
        cell_idx = [RF_overview.cell_idx];
        RF_overview = rmfield(RF_overview,'cell_idx');
        RF_overview = rmfield(RF_overview,'file');
        RF_overview = rmfield(RF_overview,'p');
        %convert to array
        RF_overview = reshape(struct2array(RF_overview),...
            numel(fieldnames(RF_overview)),[]);

        %Check if any RF were found
        table_heads{1,ii+1} = legend_name{ii}; 
        RF_overview_log(:,ii) = any(RF_overview,1)';
        table_dummy = table(cell_idx',RF_overview_log);
        table_dummy = splitvars(table_dummy);
        table_dummy.Properties.VariableNames = table_heads;
        
    end
 add_info.tables.RF_single_cell.Data = table_dummy;

end
out = 1;
end
