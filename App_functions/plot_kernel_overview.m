
function out = plot_kernel_overview (savepath, add_info)
%This script counts how many cells are alike regarding their kernels based
%on results calculated by the kernel analysis script
stim_idx = add_info.stim_idx;
nr_colours = 4;


try
     M = matfile(findfile_app(stim_idx,savepath,'cell_kernel_overview'));
     
     test = size(M.cell_kernel_overview);
     if sum(test) == 0
         out = 0; % Return 0 because field doesnt exist 
         disp('Kernel info is empty, check quality criteria')
         return 
     end
     clear M
catch 
    %If the structure doesnt exist, function gets returned here
    disp('Kernel_info not found, check if quality criteria has been calculated')
    out = 0;
    return
end


%%
S = load(findfile_app(stim_idx,savepath,'cell_kernel_overview'));
cell_kernel_overview = S.cell_kernel_overview;


%Delet the cells that didnt pass the peak test
a = 1;
for ii = 1:numel(cell_kernel_overview)
    try
       if isnan(cell_kernel_overview(a).luminance_kernel)
           cell_kernel_overview(a) = [];
       else 
           a = a+1;
       end
    end
end
    
real_cell_idx = [cell_kernel_overview.cell_idx];    
nr_cells = size(cell_kernel_overview,2);
%1. Step: Look how many luminance kernels we have this is easily done using
%the luminance_kernel array.
luminance_idx = find([cell_kernel_overview.luminance_kernel]);
luminance_cell_idx = real_cell_idx(luminance_idx);
colour_idx = find(~[cell_kernel_overview.luminance_kernel]);


result_table = table;
%calculate the precentage of luminance kernels
luminance_prop = size(luminance_idx,2)/nr_cells;
result_table.luminance_kernel_prop = luminance_prop;

luminance_array = [cell_kernel_overview(luminance_idx).active_colours];
luminance_array_text = {cell_kernel_overview(luminance_idx).kernel_type};
luminance_on = contains(luminance_array_text,'On');
luminance_off = ~luminance_on;

luminance_on_idx = luminance_cell_idx(luminance_on);
luminance_off_idx = luminance_cell_idx(luminance_off);
if ~isempty(luminance_array)
luminance_array_split = NaN(numel(luminance_array),2);
luminance_array_split(1:nnz(luminance_array(luminance_off)),1) = luminance_array(luminance_off);
luminance_array_split(1:nnz(luminance_array(luminance_on)),2) = luminance_array(luminance_on);
nr_kernel_lum = zeros(size(luminance_array_split,2),numel(luminance_array));
unique_luminance = NaN(numel(luminance_array),size(luminance_array_split,2));

for ii = 1:size(luminance_array_split,2)
    
    unique_luminance_temp = unique(luminance_array_split(:,ii));
    unique_luminance_temp = unique_luminance_temp(1:lastNaN(unique_luminance_temp));
    unique_lum_sz = numel(unique_luminance_temp);
    unique_luminance(1:unique_lum_sz,ii) = unique_luminance_temp;
    nr_comb = lastNaN(unique_luminance(:,ii));
    
   
    for kk = 1:nr_comb
        nr_kernel_lum(ii,kk) = nnz(luminance_array_split(:,ii) == unique_luminance(kk,ii));
        if ii == 1
            Cell_indices(kk).luminance_off = luminance_off_idx(1,...
                find(luminance_array_split(:,ii) == unique_luminance(kk,ii))); 
        else
            Cell_indices(kk).luminance_on = luminance_on_idx(1,...
                find(luminance_array_split(:,ii) == unique_luminance(kk,ii))); 
        end
            
      
    end
    
end


end

colour_array = [cell_kernel_overview(colour_idx).active_colours];
colour_array_text = {cell_kernel_overview(colour_idx).kernel_type};
colour_array_cell_idx = [cell_kernel_overview(colour_idx).cell_idx];

%Next we search for complex_on_off cells
complex_OO_cells = find(strcmp(colour_array_text, 'Complex_ON_OFF'));
complex_OO_idx = colour_array_cell_idx(complex_OO_cells);
if ~isempty(complex_OO_cells)
complex_array(1:nnz(colour_array(complex_OO_cells)),1) = colour_array(complex_OO_cells);
nr_kernel_complex_OO = zeros(1,numel(complex_array));


unique_complex_OO = unique(complex_array(:,1));
unique_complex_OO = unique_complex_OO(1:lastNaN(unique_complex_OO));
nr_comb = lastNaN(unique_complex_OO);

for kk = 1:nr_comb
    nr_kernel_complex_OO(1,kk) = nnz(complex_array == unique_complex_OO(kk));
    Cell_indices(kk).complex_OO = complex_OO_idx(1,...
                find(complex_array == unique_complex_OO(kk)));

end
else
    nr_kernel_complex_OO = 0;
end



%Next we search for opponent cells. Same approach

complex_oppt_cells = find(strcmp(colour_array_text, 'Complex_opp'));
complex_opp_idx = colour_array_cell_idx(complex_oppt_cells);
if ~isempty(complex_oppt_cells)
complex_array(1:nnz(colour_array(complex_oppt_cells)),1) = colour_array(complex_oppt_cells);
nr_kernel_complex = zeros(1,numel(complex_array));


unique_complex = unique(complex_array(:,1));
unique_complex = unique_complex(1:lastNaN(unique_complex));
nr_comb = lastNaN(unique_complex);

for kk = 1:nr_comb
    nr_kernel_complex(1,kk) = nnz(complex_array == unique_complex(kk));
    Cell_indices(kk).complex_oppt = complex_opp_idx(1,...
                find(complex_array == unique_complex(kk)));

end
else
     nr_kernel_complex = 0;
end
    

%Next we search UV_Blue_OP cells

UVBlue_oppt_cells = find(strcmp(colour_array_text, 'UV_Blue_OP'));
UVBlue_idx = colour_array_cell_idx(UVBlue_oppt_cells);
if ~isempty(UVBlue_oppt_cells)
UVBlue_array(1:nnz(colour_array(UVBlue_oppt_cells)),1) = colour_array(UVBlue_oppt_cells);
nr_kernel_UVBlue = zeros(1,numel(UVBlue_array));

unique_UVBlue = unique(UVBlue_array(:,1));
nr_comb = lastNaN(unique_UVBlue);

for kk = 1:nr_comb
    nr_kernel_UVBlue(1,kk) = nnz(UVBlue_array == unique_UVBlue(kk));
    Cell_indices(kk).UVBlue = UVBlue_idx(1,...
                find(UVBlue_array == unique_UVBlue(kk)));

end
else
    nr_kernel_UVBlue = 0;
end

%Next we search Blue_Green_OP cells

BlueGreen_oppt_cells = find(strcmp(colour_array_text, 'Blue_Green_OP'));
BlueGreen_idx = colour_array_cell_idx(BlueGreen_oppt_cells);
if ~isempty(BlueGreen_oppt_cells)
BlueGreen_array(1:nnz(colour_array(BlueGreen_oppt_cells)),1) = colour_array(BlueGreen_oppt_cells);
nr_kernel_BlueGreen = zeros(1,numel(BlueGreen_array));

unique_BlueGreen = unique(BlueGreen_array(:,1));
nr_comb = lastNaN(unique_BlueGreen);

for kk = 1:nr_comb
    nr_kernel_BlueGreen(1,kk) = nnz(BlueGreen_array == unique_BlueGreen(kk));
    Cell_indices(kk).BlueGreen = BlueGreen_idx(1,...
                find(BlueGreen_array == unique_BlueGreen(kk)));

end
else
    nr_kernel_BlueGreen = 0;
end

%Next we search for Green Red Op cells

GreenRed_oppt_cells = find(strcmp(colour_array_text, 'Green_Red_OP'));
GreenRed_idx = colour_array_cell_idx(GreenRed_oppt_cells);
if ~isempty(GreenRed_oppt_cells)
GreenRed_array(1:nnz(colour_array(GreenRed_oppt_cells)),1) = colour_array(GreenRed_oppt_cells);
nr_kernel_GreenRed = zeros(1,numel(GreenRed_array));

unique_GreenRed = unique(GreenRed_array(:,1));
nr_comb = lastNaN(unique_GreenRed);

for kk = 1:nr_comb
    nr_kernel_GreenRed(1,kk) = nnz(GreenRed_array == unique_GreenRed(kk));
    Cell_indices(kk).GreenRed = GreenRed_idx(1,...
                find(GreenRed_array == unique_GreenRed(kk)));

end
else
    nr_kernel_GreenRed = 0;
end


%No we have all the numbers we need. Next stage is plotting the result.
%We will use the rectangle function to draw an overview over which colours
%are active as done before by Tom. 
%First we need to crunch the numbers. Luminance cells are first on the
%left. We need one array to contain all the numbers that shall be plotted
%in a bar plot. Because we than need to sort the cells by their numbers,
%each celltype needs an unique index attached to them so that it can be found
%later in the array. Best is, we store everything in a table

plot_table = table('Size',[nr_cells,4],'VariableTypes',{'cell','double','double','double'},...
    'VariableNames',{'Type','Colour','Number','Idx'});

%Fill it, first with off cells

nr_off = find(nr_kernel_lum(1,:),1,'last');
if ~isempty(nr_off) 
plot_table.Type(1:nr_off) = {'OFF'};
%Fill in the active colours
plot_table.Colour(1:nr_off) = unique_luminance(1:nr_off,1);
%Numbers
plot_table.Number(1:nr_off) = nr_kernel_lum(1,1:nr_off);
for i = 1:nr_off
plot_table.Indices(i) = {Cell_indices(i).luminance_off};
end
nr = nr_off +1; %Index for total number
%Next, on cells
else
    nr = 1;
end

nr_on = find(nr_kernel_lum(2,:),1,'last');
if ~isempty(nr_on)
plot_table.Type(nr:nr+nr_on-1) = {'On'};
plot_table.Colour(nr:nr+nr_on-1) = unique_luminance(1:nr_on,2);
plot_table.Number(nr:nr+nr_on-1) = nr_kernel_lum(2,1:nr_on);
for i = 1:nr_on
plot_table.Indices(i+nr-1) = {Cell_indices(i).luminance_on};
end
nr = nr+nr_on;
end


%Next, complex On OFF

nr_complex_OO = find(nr_kernel_complex_OO,1,'last');
if ~isempty(nr_complex_OO) 
plot_table.Type(nr:nr+nr_complex_OO-1) = {'Complex_ON_OFF'};
plot_table.Colour(nr:nr+nr_complex_OO-1) = unique_complex_OO(1:nr_complex_OO,1);
plot_table.Number(nr:nr+nr_complex_OO-1) = nr_kernel_complex_OO(1,1:nr_complex_OO);
for i = 1:nr_complex_OO
plot_table.Indices(i+nr-1) = {Cell_indices(i).complex_OO};
end
nr = nr+nr_complex_OO;
end

%Blue_UV OPP

nr_UVBlue = find(nr_kernel_UVBlue,1,'last');
if ~isempty(nr_UVBlue)
plot_table.Type(nr:nr+nr_UVBlue-1) = {'UV_Blue_OP'};
plot_table.Colour(nr:nr+nr_UVBlue-1) = unique_UVBlue(1:nr_UVBlue,1);
plot_table.Number(nr:nr+nr_UVBlue-1) = nr_kernel_UVBlue(1,1:nr_UVBlue);
for i = 1:nr_UVBlue
plot_table.Indices(i+nr-1) = {Cell_indices(i).UVBlue};
end
nr = nr+nr_UVBlue;
end

%Blue Green OPP

nr_BlueGreen = find(nr_kernel_BlueGreen,1,'last');
if ~isempty(nr_BlueGreen)
plot_table.Type(nr:nr+nr_BlueGreen-1) = {'Blue_Green_OP'};
plot_table.Colour(nr:nr+nr_BlueGreen-1) = unique_BlueGreen(1:nr_BlueGreen,1);
plot_table.Number(nr:nr+nr_BlueGreen-1) = nr_kernel_BlueGreen(1,1:nr_BlueGreen);
for i = 1:nr_BlueGreen
plot_table.Indices(i+nr-1) = {Cell_indices(i).BlueGreen};
end
nr = nr+nr_BlueGreen;
end

%Green Red OPP

nr_GreenRed = find(nr_kernel_GreenRed,1,'last');
if ~isempty(nr_GreenRed)
plot_table.Type(nr:nr+nr_GreenRed-1) = {'Green_Red_OP'};
plot_table.Colour(nr:nr+nr_GreenRed-1) = unique_GreenRed(1:nr_GreenRed,1);
plot_table.Number(nr:nr+nr_GreenRed-1) = nr_kernel_GreenRed(1,1:nr_GreenRed);
for i = 1:nr_GreenRed
plot_table.Indices(i+nr-1) = {Cell_indices(i).GreenRed};
end
nr = nr+nr_GreenRed;
end

%Complex Colour
nr_complex = find(nr_kernel_complex,1,'last');
if ~isempty(nr_complex)
plot_table.Type(nr:nr+nr_complex-1) = {'Complex_opp'};
plot_table.Colour(nr:nr+nr_complex-1) = unique_complex(1:nr_complex,1);
plot_table.Number(nr:nr+nr_complex-1) = nr_kernel_complex(1,1:nr_complex);
for i = 1:nr_complex
plot_table.Indices(i+nr-1) = {Cell_indices(i).complex_oppt};
end

end


%Delete protruding zeros
table_last = find(plot_table.Colour,1,'last');
plot_table(table_last+1:end,:) = [];
%Create index
table_idx = (1:1:size(plot_table,1));
plot_table.Idx = table_idx';
%Sort table by number
plot_table_sorted = sortrows(plot_table,'Number','descend');

%% Plotting 
%First Barplot with the numbers
barplot1 = figure;
ax1 = subplot(2,1,1);
hold on
labelnr = (1:1:size(plot_table_sorted,1));
labelnr_str = sprintfc('%d',labelnr);
for i = 1:size(plot_table_sorted,1)
    h = bar(i,plot_table_sorted.Number(i));
    if strcmp(plot_table_sorted.Type(i),'On')
        set(h,'FaceColor','#D3D3D3');
    elseif strcmp(plot_table_sorted.Type(i),'OFF')
        set(h,'FaceColor','k');
    elseif strcmp(plot_table_sorted.Type(i),'Complex_ON_OFF')
        set(h,'FaceColor','#808080')
    else
        set(h,'FaceColor','#FFA500');
    end
   
end
ax1.XTick = labelnr;
set(gca,'xticklabel',labelnr_str);
hold off

%Second Plot with rectangles indicating the active colours
ax2 = subplot(2,1,2);
axis([0 size(plot_table_sorted,1) 1 5])
hold on
a = 1;
for i = 1:size(plot_table_sorted,1)

    colour_string = plot_table_sorted.Colour(i);
    colour_string = mat2str(colour_string);
    if contains(colour_string,"1")
        rec(a) = rectangle('Position',[i 1 0.5 1]);
        rec(a).FaceColor = 'm';
        a = a+1;
       
    end
    if contains(colour_string,"2")
        rec(a) = rectangle('Position',[i 1 0.5 1]);
        rec(a).FaceColor = 'k';
        a = a+1;
       
    end
    if contains(colour_string, "3")
        rec(a) = rectangle('Position',[i 2 0.5 1]);
        rec(a).FaceColor = 'b';
        a = a+1;
        
    end
    if contains(colour_string,"4")
        rec(a) = rectangle('Position',[i 2 0.5 1]);
        rec(a).FaceColor = 'k';
        a = a+1;
       
    end
    if contains(colour_string,"5")
        rec(a) = rectangle('Position',[i 3 0.5 1]);
        rec(a).FaceColor = 'g';
        a = a+1;
        
    end
    if contains(colour_string,"6")
        rec(a) = rectangle('Position',[i 3 0.5 1]);
        rec(a).FaceColor = 'k';
        a = a+1;
       
    end
    if contains(colour_string,"7")
        rec(a) = rectangle('Position',[i 4 0.5 1]);
        rec(a).FaceColor = 'r';
        a = a+1;
        
    end
    if contains(colour_string,"8")
        rec(a) = rectangle('Position',[i 4 0.5 1]);
        rec(a).FaceColor = 'k';
        a = a+1;
       
    end
     
           
      
   a = 1;    
   clear rec
end


hold off
ax2.Position = ax2.Position + [0 0.25 0 -0.2];
linkaxes([ax1,ax2],'x')

ylabel(ax1,'Cells')
ax2.XTickLabel = [];
ylabel(ax2,'Active Colours');


%Next we need to plot the average kernels 
%Preallocate structure for that we need to preload a Kernel to get the size
%of the Kernels
S = load(findfile_app(stim_idx,savepath,'Kernel_location'));
kernel_location = S.Kernel_location;

S = load(kernel_location(1));
% kernel_location = Stimulus_info(stim_idx).Kernel_location;
% S = load(kernel_location(1));
Kernel = S.Kernels;
Kernel_size = size(Kernel);
mean_traces = zeros(size(plot_table_sorted,1),Kernel_size(1),Kernel_size(3));
nr_subplots = numSubplots(size(plot_table_sorted,1));
kernel_figure = figure;
set(gcf,'color','w');
for i = 1:size(plot_table_sorted,1)
    %preallocate matrix
    
    
   Cell_indices = plot_table_sorted.Indices{i,1};
   plot_traces = zeros(length(Cell_indices),Kernel_size(1),Kernel_size(3));
   for ii = 1:length(Cell_indices)
      table_position = [cell_kernel_overview.cell_idx] == Cell_indices(ii);
      channel_idx = cell_kernel_overview(table_position).channel_idx;
      %Load the Kernel into memory
      load_idx = str2double(kernel_location(:,2)) == Cell_indices(ii);
      S = load(kernel_location(load_idx,1));
      Kernel = S.Kernels;
      plot_traces(ii,:,:)  = squeeze(Kernel(:,channel_idx,:));
      
            
    
     
      
   end
   
   
   mean_traces(i,:,:) = mean(plot_traces,1,'omitnan');
   mean_trace = squeeze(mean_traces(i,:,:));
   % normalize
   %Normalize trace
    ac_length = length(mean_trace);
    ac_length_half = ac_length/2;
    ac_future_idx = ceil(ac_length*(2/3));
    xvalues = (-ac_length_half+1:1:ac_length_half);
    active_channel_future = (squeeze(mean_trace(ac_future_idx:end,:)));
    active_channel_future_mean = mean(active_channel_future,2);
    ac_mean = mean(active_channel_future_mean);
    ac_std = std(active_channel_future_mean);
    mean_traces_norm = znormalise(mean_trace,ac_mean,ac_std);

   
   
   
   %mean_traces_norm = squeeze(mat2gray(mean_traces(i,:,:)));
   
   
   
   colour_string = {'m','b','g','r'};
   for cc = 1:Kernel_size(3)
      subplot(nr_subplots(1),nr_subplots(2),i);
      plot(xvalues,mean_traces_norm(:,cc),colour_string{cc})
      hold on
       
   end
   set(gca,'box','off')
   h = gca;
   h.YAxis.Visible = 'off';
   xline(0);
   set(gca,'ytick',[]);
   title(num2str(i));
   hold off
    
end


%% Calculate plot values
ON_dominance = NaN(Kernel_size(3), length(cell_kernel_overview));
OFF_dominance = NaN(Kernel_size(3),length(cell_kernel_overview));

for ii = 1:length(cell_kernel_overview)
    
   %extract the table
   detailed_info = cell_kernel_overview(ii).detailed_info;
   ON_dominance(:,ii) = detailed_info.ON_dominance;
   OFF_dominance(:,ii) = detailed_info.OFF_dominance;
   for cc = 1:nr_colours
       if detailed_info.kernel_type_log{cc} == 1
           ON_peaks(cc,ii) = detailed_info.On_peaks(cc);
           OFF_peaks(cc,ii) = NaN;
           OFF_freq(cc,ii) = NaN;
           ON_freq(cc,ii) = detailed_info.On_frequency(cc);
       elseif detailed_info.kernel_type_log{cc} == -1
           ON_peaks(cc,ii) = NaN;
           OFF_peaks(cc,ii) = detailed_info.Off_peaks(cc);
           OFF_freq(cc,ii) = detailed_info.Off_frequency(cc);
           ON_freq(cc,ii) = NaN;
       else
           ON_peaks(cc,ii) = NaN;
           OFF_peaks(cc,ii) = NaN;
           OFF_freq(cc,ii) = NaN;
           ON_freq(cc,ii) = NaN;
       end
           
   end
    
end

nr_histbins = 50;
hist_step = 1/nr_histbins;
x_bar = (0:hist_step:1);


%Check whats the max value for the peaks and define the edges for the
%histogram
max_peaks = max(ON_peaks,[],'all');
min_peaks = max(OFF_peaks,[],'all');
max_peaks_all = ceil(max([max_peaks,min_peaks]))+20;
bin_amplitudes = (0:2:max_peaks_all);

%Same for the Frequency data
max_peaks_fq = max(ON_freq,[],'all');
min_peaks_fq = max(OFF_freq,[],'all');
max_freq_all = ceil(max([max_peaks_fq,min_peaks_fq]))+20;
bin_amplitudes_fq = (0:1:max_freq_all);

hist_off = NaN(Kernel_size(3),nr_histbins);
hist_on = NaN(Kernel_size(3),nr_histbins);
hist_peaks_off = NaN(Kernel_size(3),length(bin_amplitudes)-1);
hist_peaks_on = NaN(Kernel_size(3),length(bin_amplitudes)-1);
hist_freq_off = NaN(Kernel_size(3),length(bin_amplitudes_fq)-1);
hist_freq_on = NaN(Kernel_size(3),length(bin_amplitudes_fq)-1);
for ii = 1:Kernel_size(3)
   
    hist_off(ii,:) = histcounts(OFF_dominance(ii,:),nr_histbins);
    hist_on(ii,:) = histcounts(ON_dominance(ii,:),nr_histbins);
    hist_peaks_off(ii,:) = histcounts(-OFF_peaks(ii,:),sort(-bin_amplitudes));
    hist_peaks_on(ii,:) = histcounts(ON_peaks(ii,:),bin_amplitudes);
    hist_freq_off(ii,:) = histcounts(OFF_freq(ii,:),sort(bin_amplitudes_fq));
    hist_freq_on(ii,:) = histcounts(ON_freq(ii,:),sort(bin_amplitudes_fq));
    
end


%Transform the respective arrays for plotting into GPU arrays, this
%increases the speed and performance of figures

hist_off = gpuArray(hist_off);
hist_on = gpuArray(hist_on);
hist_peaks_off = gpuArray(hist_peaks_off);
hist_peaks_on = gpuArray(hist_peaks_on);
hist_freq_off = gpuArray(hist_freq_off);
hist_freq_on = gpuArray(hist_freq_on);

%% Plotting

ax_sub_hist_dominance = cell(Kernel_size(3),2);
ax_sub_hist_peaks = cell(Kernel_size(3),2);
ax_sub_hist_freq = cell(Kernel_size(3),2);
ax_hist_dominance = cell(Kernel_size(3),2);
ax_hist_peaks = cell(Kernel_size(3),2);
ax_hist_freq = cell(Kernel_size(3),2);
dominance_figure = figure;
peaks_figure = figure;
freq_figure = figure;
a = 1;
for ii = 1:Kernel_size(3)
    %Here we create subplots in the following order: dominance counts, peaks,
    %frequency
    %OFF
    set(0,'CurrentFigure',dominance_figure)
    ax_sub_hist_dominance{ii,1} = subplot(Kernel_size(3),2,a);
    ax_hist_dominance{ii,1} = histogram('BinCounts',hist_off(ii,:),'BinEdges',x_bar,'FaceColor',colour_string{ii});
    set(gca,'Color','#D3D3D3')
    set(0,'CurrentFigure',peaks_figure)
    ax_sub_hist_peaks{ii,1} = subplot(Kernel_size(3),2,a);
    ax_hist_peaks{ii,1} = histogram('BinCounts',hist_peaks_off(ii,:),'BinEdges',sort(-bin_amplitudes),'FaceColor',colour_string{ii});
    set(gca,'Color','#D3D3D3')
    set(0,'CurrentFigure',freq_figure)
    ax_sub_hist_freq{ii,1} = subplot(Kernel_size(3),2,a);
    ax_hist_freq{ii,1} = histogram('BinCounts',hist_freq_off(ii,:),'BinEdges',bin_amplitudes_fq,'FaceColor',colour_string{ii},'EdgeColor','none');
    set(gca,'Color','#D3D3D3')
    %ON
    a = a+1;
    set(0,'CurrentFigure',dominance_figure)
    ax_sub_hist_dominance{ii,2} = subplot(Kernel_size(3),2,a);
    ax_hist_dominance{ii,2} = histogram('BinCounts',hist_on(ii,:),'BinEdges',x_bar,'FaceColor',colour_string{ii});
    set(0,'CurrentFigure',peaks_figure)
    ax_sub_hist_peaks{ii,2} = subplot(Kernel_size(3),2,a);
    ax_hist_peaks{ii,2} = histogram('BinCounts',hist_peaks_on(ii,:),'BinEdges',bin_amplitudes,'FaceColor',colour_string{ii});
    set(0,'CurrentFigure',freq_figure)
    ax_sub_hist_freq{ii,2} = subplot(Kernel_size(3),2,a);
    ax_hist_freq{ii,2} = histogram('BinCounts',hist_freq_on(ii,:),'BinEdges',bin_amplitudes_fq,'FaceColor',colour_string{ii},'EdgeColor','none');
    a = a+1;    
end

%Make the graph look nicer
%Add labels
set(ax_sub_hist_dominance{1,1}.Title,'String','OFF')
set(ax_sub_hist_dominance{1,2}.Title,'String','ON')


set(ax_sub_hist_peaks{1,1}.Title,'String','OFF')
set(ax_sub_hist_peaks{1,2}.Title,'String','ON')

set(ax_sub_hist_freq{1,1}.Title,'String','OFF')
set(ax_sub_hist_freq{1,2}.Title,'String','ON')

set(ax_sub_hist_dominance{end,1}.XLabel,'String','Dominance (1 = stronges colour)')
set(ax_sub_hist_dominance{end,2}.XLabel,'String','Dominance (1 = stronges colour)')
set(ax_sub_hist_dominance{end,1}.YLabel,'String','Cells')

set(ax_sub_hist_peaks{end,1}.XLabel,'String','Amplitude (SD)')
set(ax_sub_hist_peaks{end,2}.XLabel,'String','Amplitude (SD)')
set(ax_sub_hist_peaks{end,1}.YLabel,'String','Cells')


set(ax_sub_hist_freq{end,1}.XLabel,'String','Frequency of response [Hz]')
set(ax_sub_hist_freq{end,2}.XLabel,'String','Frequency of response [Hz]')
set(ax_sub_hist_freq{end,1}.YLabel,'String','Cells')


%Link axes
linkaxes([ax_sub_hist_dominance{:}],'xy');
linkaxes([ax_sub_hist_peaks{:,1}],'x');
linkaxes([ax_sub_hist_peaks{:,2}],'x');
linkaxes([ax_sub_hist_peaks{:,:}],'y');
linkaxes([ax_sub_hist_freq{:}],'xy');
ax_sub_hist_freq{1,1}.XLim = [0 20];
%Change xlim
%Since axes are linked we just need to change xlim for one axis



%Add title
sgtitle(dominance_figure,'Colour dominance across cells')
sgtitle(peaks_figure,'Response amplitude')
sgtitle(freq_figure,'Response frequency')



for ii = 1:Kernel_size(3)
    if ii ~= Kernel_size(3)
    ax_sub_hist_dominance{ii,1}.XTick = [];
    ax_sub_hist_dominance{ii,2}.XTick = [];
    
    ax_sub_hist_peaks{ii,1}.XTick = [];
    ax_sub_hist_peaks{ii,2}.XTick = [];
    
    ax_sub_hist_freq{ii,1}.XTick = [];
    ax_sub_hist_freq{ii,2}.XTick = [];
    
    end
    ax_sub_hist_dominance{ii,2}.YTick = [];
    ax_sub_hist_peaks{ii,2}.YTick = [];
    ax_sub_hist_freq{ii,2}.YTick = [];
    
            
end

autoArrangeFigures(1,3,1);



%% Save figures
%This part saves the created figures as a .fig file and collects the
%savename afterwards so the Datafile for that stimulus is updated

%First check if a subfolder for figures exists in the stimulus folder

out = sf_organizer(stim_idx,savepath,'subfoldername',"Figures");
cd(out);
savefig([barplot1,kernel_figure,dominance_figure,peaks_figure,freq_figure],'Kernel_overview');

[~] = sf_organizer(stim_idx,savepath,'subfoldername',"Figures",'collect_files',true);


%% PLotting average receptive field sizes
%Coming soon

%Check if receptive field size information is available

if isfield(cell_kernel_overview,'Receptive_field_overview')

    
%Add the receptive field sizes for all cells and all clusters
 error_figure = figure;
 size_figure = figure;
 
for i = 1:height(plot_table_sorted)
    %preallocate matrix
       
   Cell_indices = plot_table_sorted.Indices{i,1};
   rf_size = zeros(length(Cell_indices),Kernel_size(3));  
   rf_amplitude = zeros(length(Cell_indices),Kernel_size(3));  
   rf_error = zeros(length(Cell_indices),Kernel_size(3));   
   
   for ii = 1:length(Cell_indices)
      table_position = [cell_kernel_overview.cell_idx] == Cell_indices(ii);
      fitresult = cell_kernel_overview(table_position).Receptive_field_overview.fitresult;
      fiterror = cell_kernel_overview(table_position).Receptive_field_overview.fiterr;
      for cc = 1:Kernel_size(3)
          %This calculates the size assuming ellipse shape
          rf_size(ii,cc) = fitresult.Size_x(cc)*fitresult.Size_y(cc)*pi;
          rf_amplitude(ii,cc) = fitresult.Amplitude(cc);
          rf_error(ii,cc) = fiterror(cc,3)*fiterror(cc,4)*pi;
      end
      
               
   end
   mean_rf_color = mean(rf_size,1,'omitnan');
   std_rf_color = std(rf_size,[],1,'omitnan');
   
   
   
   plot_table_sorted.rf_size{i} = rf_size;
   plot_table_sorted.mean_rf_color{i} = mean_rf_color;
   plot_table_sorted.std_rf_color{i} = std_rf_color;
   plot_table_sorted.rf_amplitude{i} = rf_amplitude;
   plot_table_sorted.rf_error{i} = rf_error;
  
   set(0,'CurrentFigure',size_figure)
   hold on
   sp(i) = subplot(nr_subplots(1),nr_subplots(2),i);
%    rf_idx = false(4,1);
%    for kk = 1:4
%        if ~isnan(plot_table_sorted.std_rf_color{i}(kk)) 
%             rf_idx(kk) = true;
%        end
%    end
%      
%    [h1,L1,~,~] = violin(plot_table_sorted.rf_size{i,1}(:,rf_idx),'bw',0.3,'facecolor','k');
   h = boxplot(plot_table_sorted.rf_size{i,1},'Notch','on','Colors',[colour_string{:}]);
%    nr_colors = find(rf_idx);
%     for kk = 1:nnz(rf_idx)
%         
%        h1(1).FaceColor = colour_string{nr_colors(kk)};
%        
%     end 
   
   
   set(0,'CurrentFigure',error_figure)
   hold on
   subplot(nr_subplots(1),nr_subplots(2),i);
   
   
%    [h2,L2,~,~] = violin(plot_table_sorted.rf_error{i,1}(:,rf_idx),'bw',0.3,'FaceColor','k');
  g = boxplot(plot_table_sorted.rf_error{i,1},'Notch','on','Colors',[colour_string{:}]);
   
%      for kk = 1:nnz(rf_idx)
%         
%        h2(1).FaceColor = colour_string{nr_colors(kk)};
%        
%     end 
%    
       
end
   

set(0,'CurrentFigure',size_figure)
sgtitle('RF Sizes');
set(0,'CurrentFigure',error_figure)
sgtitle('Error of fit');
linkaxes([sp],'y')

%% Save figures
%This part saves the created figures as a .fig file and collects the
%savename afterwards so the Datafile for that stimulus is updated

%First check if a subfolder for figures exists in the stimulus folder

out = sf_organizer(stim_idx,savepath,'subfoldername',"Figures");
cd(out);
savefig([size_figure,error_figure],'RF_overview');

[~] = sf_organizer(stim_idx,savepath,'subfoldername',"Figures",'collect_files',true);



end








    
end









