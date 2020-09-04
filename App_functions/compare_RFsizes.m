function out = compare_RFsizes (savepath,add_info)
%This script compares different metrics of receptive fields (size etc)


stim_idx = add_info.stim_idx;
nr_stimuli = length(stim_idx);
mm_per_pixel = [0.024,0.151,0.054];

limit = 10;
colour_string = {'m','b','g','r'};
S = load(savepath,'stim_list');

stim_list = strrep(S.stim_list,'_',' '); %Remove _ from string for plotting

%% Loop over stimuli selected
for ss = 1:nr_stimuli
    
%Load example Kernel to get the size
S = load(findfile_app(stim_idx(ss),savepath,'Kernel_location'));
kernel_location = S.Kernel_location;
S = load(kernel_location(1,1));
Kernels = S.Kernels;

Kernel_size = size(Kernels);


%Load fiel with the metadata
S = load(findfile_app(stim_idx(ss),savepath,'cell_kernel_overview'));

cell_kernel_overview = S.cell_kernel_overview;
nr_cells = length(cell_kernel_overview);


rf_size = NaN(nr_cells,Kernel_size(3),2);
rf_speed = NaN(nr_cells,Kernel_size(3),2);
for ii = 1:nr_cells
    %find active colours and get receptive field sizes for on and off
    %cells independently
    kernel_type_log = [cell_kernel_overview(ii).detailed_info.kernel_type_log{:}];
    RF = cell_kernel_overview(ii).Receptive_field_overview;



for ll = 1:Kernel_size(3)
    if kernel_type_log(ll) == 1
        rf_size(ii,ll,1) = RF.fitresult.Size_x(ll) * RF.fitresult.Size_y(ll)*pi;
        

    elseif kernel_type_log(ll) == -1
        rf_size(ii,ll,2) = RF.fitresult.Size_x(ll) * RF.fitresult.Size_y(ll)*pi;
    else
        continue 
    end

end





end


histo_sizes(ss) = figure;
edges = (0:0.5:limit);
x_edges = edges*mm_per_pixel(ss);

for ii = 1:Kernel_size(3)

    [N,~] = histcounts(rf_size(:,ii,1),edges);
    %Normalize
    N = N/max(N);
    subplot(1,2,1)
    ax(ii) = plot(x_edges(2:end),N,'Color',colour_string{ii});
    hold on
    KS_test(ss,ii,1) = kstest(rf_size(:,ii));
end
title('ON')
ylim([0,1.3])
xlim([10e-3,1])
set(gca, 'XScale', 'log')
grid on
hold off



for ii = 1:Kernel_size(3)

    [N1,~] = histcounts(rf_size(:,ii,2),edges);
    %Normalize
    N1 = N1/max(N1);
    subplot(1,2,2)
    ax1(ii) = plot(x_edges(2:end),N1,'Color',colour_string{ii});
    hold on
    KS_test(ss,ii,2) = kstest(rf_size(:,ii));
    
end
ylim([0,1.3])
xlim([10e-3,1])
set(gca, 'XScale', 'log')
title('OFF')
grid on

hold off
sgtitle(['RF sizes all cells, stimulus ',stim_list{stim_idx(ss)}]);


%% Calculate correlations between colors in the same cell

%Calculate number of possible colour combinations
combinations = combnk(1:1:Kernel_size(3),2);
combinations_temp = combinations;
combinations(3,:) = combinations_temp(4,:);
combinations(4,:) = combinations_temp(3,:);

combinations_text = string(combinations);

combinations_text = strrep(combinations_text,'1','UV');
combinations_text = strrep(combinations_text,'2','Blue');
combinations_text = strrep(combinations_text,'3','Green');
combinations_text = strrep(combinations_text,'4','Red');

combinations_text = strcat(combinations_text(:,1),' vs. ',combinations_text(:,2));

nr_combinations = length(combinations);


%Create matrix to store combinations
rf_combinations = NaN(nr_cells,nr_combinations,2);
OnOfftext = {'On', 'Off'};
for kk = 1:2
figure
for ii = 1:nr_combinations
    
        
        rf_combinations(:,ii,:) = rf_size(:,combinations(ii,1),:)./rf_size(:,combinations(ii,2),:);
        subplot(3,3,ii)
        histogram(rf_combinations(:,ii,kk),edges,'Normalization','pdf');
        xlabel(combinations_text(ii));
        xlim([0 5])
        ylim([0 1])
        xticks([1 2 3 4 5])
end
sgtitle(['Relative RF sizes ',OnOfftext{kk},' cells, Stimulus ' ,stim_list{stim_idx(ss)}]); 
end


rf_combinations_mean = squeeze(nanmean(rf_combinations,1));

for kk = 1:2
    bx(kk) = figure;
    boxplot(-log2(rf_combinations(:,:,kk)),combinations_text);
    title(['Relative RF sizes ',OnOfftext{kk},' cells, Stimulus ' ,stim_list{stim_idx(ss)}]); 
    ylim([-3,3])
    yticks([-3 -2 -1 0 1 2 3])
    
   
end



out = 1;
end
autoArrangeFigures
end