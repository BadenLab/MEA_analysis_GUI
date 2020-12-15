function out = plotoverview_kernel_new (savepath, add_info)

%This function searches for overview variables for all the different
%possibilities of calculating the RF from the gui and plots an overview, so
%that the user knows which cells actually show significant kernels and can
%further inspect single cells. Paul's RF analysis code has to have run
%before otherwise nothing is plotted.
if length(add_info.stim_idx) > 4
    error("Only 4 different settings can be compared at a time")
end

table_heads = {'Cell idx'}; 

M = matfile(savepath,'Writable',false);
RF_overview_log = false(length(M.cell_indices),length(add_info.settings.kernel_new.folderplot));

for ii = 1:length(add_info.settings.kernel_new.folderplot)
Stimulus_info = M.Stimulus_info(1,add_info.stim_idx);
Stimulus_info = Stimulus_info{:};
cut_path = strfind(Stimulus_info,'\');

path_to_data = [Stimulus_info(1:cut_path(end)),add_info.settings.kernel_new.folderplot{:,ii}];

legend_name{ii} = strrep(add_info.settings.kernel_new.folderplot{:,ii},'_',' ');

%% Get the data from the files
try
S = load([path_to_data,'\RF_overview.mat']);
RF_overview = S.RF_overview;
clear S

%Split off cell idx
cell_idx = [RF_overview.cell_idx];
RF_overview = rmfield(RF_overview,'cell_idx');
RF_overview = rmfield(RF_overview,'file');
%convert to array
RF_overview = reshape(struct2array(RF_overview),...
    numel(fieldnames(RF_overview)),[]);

%Check if any RF were found
RF_overview_log(:,ii) = any(RF_overview,1)';

catch
    error("Couldnt load file, check folder");
end

table_heads{1,ii+1} = legend_name{ii}; 

end
%Plotting
%dummy array for barplot
RF_overview_log_plot = ones(size(RF_overview_log,1),size(RF_overview_log,2));
%Delete old plot
Cell_overview_panel = add_info.panels.DatasetoverviewPanel;
delete(Cell_overview_panel.Children);
%plot new plot as grouped barplot
ax1 = subplot(1,1,1,'parent',Cell_overview_panel);
b = bar(ax1,cell_idx,RF_overview_log_plot,'stacked','EdgeColor','none','FaceColor','flat');
%Dummy bar for the legend
%b_c = bar(ax1,cell_idx,RF_overview_log_plot,'stacked','EdgeColor','none','FaceColor','flat');
%Colour for the bar
colour_array = [0 0.45 0.75; 0.85 0.32 0.099; 0.93 0.69 0.125; 0.49 0.18 0.56];

for i = 1:length(b)
    for ii = 1:length(RF_overview_log)
        if RF_overview_log(ii,i)
            ax1.Children(i).CData(ii,:) = colour_array(i,:);
        else
            ax1.Children(i).CData(ii,:) = [1, 1, 1]; %white
        end
    end
end

legend(ax1,legend_name,'Location', 'southoutside');

ax1.XLabel.String = "Cell Idx";
ax1.YTickLabel = [];

%Update single cell overview table
table_dummy = table(cell_idx',RF_overview_log);
table_dummy = splitvars(table_dummy);
table_dummy.Properties.VariableNames = table_heads;

add_info.tables.RF_single_cell.Data = table_dummy;
add_info.tables.RF_single_cell.ColumnName = table_heads;
out = 1;


%% Return information about receptive field type
%Only if only one recording was selected
if length(add_info.settings.kernel_new.folderplot) == 1
   if size(RF_overview,1) == 1
       set(add_info.panels.RFtypePanel_2.Children(3,1),'enable', 'on')
       set(add_info.panels.RFtypePanel_2.Children(1:2,1),'enable', 'off')
   elseif size(RF_overview,1) == 2
       set(add_info.panels.RFtypePanel_2.Children(2:3,1),'enable', 'on')
       set(add_info.panels.RFtypePanel_2.Children(1,1),'enable', 'off')
   else
       set(add_info.panels.RFtypePanel_2.Children(:,1),'enable', 'on')
   end
          
    
end





end