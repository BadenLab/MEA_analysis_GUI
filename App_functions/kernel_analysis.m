function out = kernel_analysis(savepath, add_info)
nr_peaks = 2;
stim_idx = add_info.stim_idx;
thr = 4;

try
     M = matfile(findfile_app(stim_idx,savepath,'Kernel_info'));
     
     test = size(M.Kernel_info);
     S = load(savepath,'cell_indices');
     cell_indices = S.cell_indices;
     if sum(test) == 0
         out = 0; % Return 0 because field doesnt exist 
         disp('Kernel info is empty, check quality criteria')
         return 
     end
     
catch 
    %If the structure doesnt exist, function gets returned here
    disp('Kernel_info not found, check if quality criteria has been calculated')
    out = 0;
    return
end
    
    S = load(findfile_app(stim_idx,savepath,'Kernel_location'));
    Kernel_locations = S.Kernel_location;
    S = load(findfile_app(stim_idx,savepath,'Kernel_info'));
    Kernel_info = S.Kernel_info;
    
    
    %Create log array for cells that have passed the QC
    true_kernel_log = zeros(size(Kernel_info,2),1);
    for ii = 1:size(Kernel_info,2)
    true_kernel_log(ii) = any(Kernel_info(ii).true_kernel_log);
    end
    
    %% Loop
    true_kernel_idx = find(true_kernel_log);
    main_channels = NaN(numel(true_kernel_idx),1);
    a = 1;
     %preallocate the strcuture
    for ss = 1:numel(true_kernel_idx)
        cell_kernel_overview(ss).luminance_kernel = NaN;

    end
    
    for ii = true_kernel_idx'
       ii
       % Get the most active channel
       try
       S = load(Kernel_locations(ii,1));
       catch
           continue
       end
       Kernels = S.Kernels;
       Kernels_size = size(Kernels);
       %Create array of channels identified by the QC function
       [~,~,Channels_idx] = find(Kernel_info(ii).true_kernel_idx);
       Channels_idx = unique(Channels_idx);
       clear S
       
       if numel(Channels_idx) ~= 1
            main_channels(a) = most_active(Kernels,Channels_idx);
       else 
           main_channels(a) = Channels_idx;
       end
        a = a+1;
        
        %Next step is to identify what kind of kernel we have, is there
        %opponency or not, which colour is off, which is on,
        %which colour is most active?
        
        %First we identify if a colour active and is on or of
        
        %Load the respective channel
        try
        active_channel = squeeze(Kernels(:,main_channels(a-1),:));
        catch
            continue
        end
        %Create 4 different versions of the kernel, one normalized to
        %center, one mirrored around the x axis (to be able to use the
        %findpeak function)
        
        %Take future part of the kernel to calculate the mean and std of
        %the noise
        ac_length = length(active_channel);
        ac_future_idx = ceil(ac_length*(2/3));
        active_channel_future = (active_channel(ac_future_idx:end,:));
        active_channel_future_mean = mean(active_channel_future,2);
        ac_mean = mean(active_channel_future_mean);
        ac_std = std(active_channel_future_mean);
        active_channel_norm = znormalise(active_channel,ac_mean,ac_std);
        
        
        %active_channel_norm = normalize(active_channel,1,'center');
        %active_channel_norm = normalize(active_channel,1,'zscore');
        active_channel_abs = -active_channel_norm;
        
        %smooth the traces to avoid sharp peaks being detected as peaks
        active_channel_norm = movmean(active_channel_norm,10,1);
        active_channel_abs = movmean(active_channel_abs,10,1);
        
        %Set threshold for what counts as a peak
        half_way = size(active_channel,1)/2;
        active_channel_std = std(active_channel_norm(half_way:end,:),1);
        active_channel_std_thr = active_channel_std*thr; %threshold = times std
        pks_save = NaN(Kernels_size(3)*2,2);
        locs_save = NaN(Kernels_size(3)*2,2);
        width_save = NaN(Kernels_size(3)*2,2); 
        %Cut the beginning of the active channel, because we are just
        %interested in what happend before the spike not after
        active_channel_norm = active_channel_norm(1:Kernels_size(1)/2,:);
        channel_std = zeros(Kernels_size(3),1);
        b=1;
        for ki = 1:Kernels_size(3)
        [pks,locs,width] = findpeaks(active_channel_norm(:,ki),'NPeaks',nr_peaks,...
            'MinPeakProminence',active_channel_std_thr(ki),'MinPeakDistance',10 ,'Annotate','extents');
        [pks_neg,locs_neg,width_neg] = findpeaks(active_channel_abs(1:250,ki),'NPeaks',nr_peaks,...
            'MinPeakProminence',active_channel_std_thr(ki),'MinPeakDistance',10, 'Annotate','extents');
        nr_pks = numel(pks);
        nr_pks_neg = numel(pks_neg);
        nr_width = numel(width);
        nr_width_neg = numel(width_neg);
        %Calculate the standard deviation for each colour, so we can
        %calculate the peak amplitude in units of std.
        
        
        
        try
            pks_save(b:b+nr_pks-1,1) = pks;
            width_save(b:b+nr_width-1,1) = width;
            locs_save(b:b+nr_pks-1,1) = locs;
            
        catch
           
        end
        try
            
            locs_save(b:b+nr_pks_neg-1,2) = locs_neg;
            pks_save(b:b+nr_pks_neg-1,2) = pks_neg;
            width_save(b:b+nr_width_neg-1,2) = width_neg;
            b = b+2;
        catch
            b = b+2;
        end
        end
        
        
        
        
        
        
        %Find out which of the colours is the most active 
        [M,I_sign] = max(pks_save,[],2,'omitnan');
        [~,I] = max(M,[],'omitnan');
        %Which colour is this?
        max_colour = ceil(I/nr_peaks);
        %Store in logical array
        colour_dominant = false(Kernels_size(3),1);
        colour_dominant(max_colour) = true;
        
        %Calculate how dominant this colour is over others
        colour_dominance = pks_save/pks_save(I,I_sign(I));
       
        %Check which of the two possible peaks is the largest
        p1 = 1;
        p2 = nr_peaks;
        colour_dominance_sg = NaN(Kernels_size(3),nr_peaks);
        colour_peaks_sg = NaN(Kernels_size(3),nr_peaks);
        colour_width = NaN(Kernels_size(3),nr_peaks);
        peaks_position = NaN(Kernels_size(3),nr_peaks);
        for cc = 1:Kernels_size(3)
            
        [colour_dominance_sg(cc,:),I1] = max(colour_dominance(p1:p2,:),[],1,'omitnan');
        I1 = I1 + [-1+p1, -1+p1];
%       [colour_peaks_sg(cc,:),~] = max(pks_save(p1:p2,:),[],1,'omitnan');
        colour_peaks_sg(cc,1) =  pks_save(I1(1),1);
        colour_peaks_sg(cc,2) =  pks_save(I1(2),2);
        colour_width(cc,1)  =  width_save(I1(1),1);
        colour_width(cc,2)  = width_save(I1(2),2);
        peaks_position(cc,1) = locs_save(I1(1),1);
        peaks_position(cc,2) = locs_save(I1(1),2);
        p1 = p2+1;
        p2 = p2+nr_peaks;
        end
       
        %Create table to store all the information
        colour = {'UV';'Blue';'Green';'Red'};
        ON_dominance = colour_dominance_sg(:,1);
        OFF_dominance = colour_dominance_sg(:,2);
        On_peaks = colour_peaks_sg(:,1);
        Off_peaks = colour_peaks_sg(:,2);
        On_width = colour_width(:,1);
        Off_width = colour_width(:,2);
        position_on = peaks_position(:,1);
        position_off = peaks_position(:,2);
        thousand = 1000;
        On_frequency = thousand./(2*On_width);
        Off_frequency = thousand./(2*Off_width);
        table_all = table(colour,colour_dominant,ON_dominance,OFF_dominance,...
            On_peaks,Off_peaks,On_width,Off_width,On_frequency,Off_frequency,...
            position_on, position_off);
        
        
        %Next we check which kind of cell we have, On, ON-OFF, OFF...
        %For this we need the locs calculated earlier
        p1 = 1;
        p2 = nr_peaks;
        kernel_type = cell(Kernels_size(3),1);
        kernel_type_log = num2cell(zeros(Kernels_size(3),1));
        
        for cc = 1:Kernels_size(3)
         
            
            locs_temp = locs_save(p1:p2,:);
            locs_temp(isnan(locs_temp)) = 0;
            locs_temp = locs_temp(:);
            if all(locs_temp == 0) 
                p1 = p2+1;
                p2 = p2+nr_peaks;
                continue
            end
            %Add logical indices to it so we can sort it
            locs_temp(1:nr_peaks,2) = 1;
            locs_temp(nr_peaks+1:end,2) = 0;
            [~,IN] = sort(locs_temp(:,1));
            locs_sorted = locs_temp(IN,:);
            %delete empty locs
            locs_sorted(locs_sorted(:,1)==0,:) = [];
            
            %Create descriptive matrix
            locs_desc = cell(size(locs_sorted,1),1);
            locs_desc(locs_sorted(:,2) == 1) = {'On'};
            locs_desc(locs_sorted(:,2) == 0) = {'OFF'};
                     
            locs_desc = sprintf('%s_', locs_desc{:});
            locs_desc(end)=[];
            %store as numbers as well
            locs_log = cell(size(locs_sorted,1),1);
            locs_log(locs_sorted(:,2) == 1) = {'1'};
            locs_log(locs_sorted(:,2) == 0) = {'-1'};
                     
            locs_log = sprintf('%s', locs_log{:});
%             locs_log(end)=[];
            
            if length(locs_log) >= 2
                if contains(locs_log(end-1),'-')
                    locs_log = -1;
                else
                locs_log = str2double(locs_log(end));
                end
            else 
                locs_log = str2double(locs_log);
            end
            
            kernel_type{cc} = locs_desc;
            kernel_type_log{cc} = locs_log;
            %save in table
            p1 = p2+1;
            p2 = p2+nr_peaks;
        end
            %Update the table with the results
            table_all.kernel_type = kernel_type;
            table_all.kernel_type_log = kernel_type_log;
            
            
           
            Idx1 = find(~cellfun('isempty',table_all.kernel_type));
            table_all.active_colour = zeros(Kernels_size(3),1);
            table_all.active_colour(Idx1) = 1;
            out(a-1).table_all = table_all;
            out(a-1).kernel_idx = true_kernel_idx(a-1);
            
            
            %Next we want to check the kernels for colour opponencies
            %We will only find cells with on/off opponencies and complex
            %cells which show ON_OFF_ON kernels etc and normal on or off
            %cells which show same kernels for all colours
            
            %If we mutiply the numbers, stored in kernel_type_log we get
            %only 2 as output whenever a cell is opponent and we only take
            %the unique kernel types
            kernel_type_log_1 = cell2mat(kernel_type_log);
            opponent_log = num2cell(zeros(Kernels_size(3),1));                     
            %Get the luminance only sensitive kernels (all values are the
            %same)
            %Filter out complex cells
            [kernel_unique] = unique(kernel_type_log_1);
            kernel_unique(kernel_unique == 0) = [];
            
           
                                     
            
            if numel(kernel_unique) == 1
                
               
                cell_kernel_overview(a-1).luminance_kernel = 1;
                Idx = find(~cellfun('isempty', kernel_type));
                if kernel_unique == 1
                cell_kernel_overview(a-1).kernel_type = 'On';
                cell_kernel_overview(a-1).kernel_type_log = 1;
                elseif kernel_unique == -1
                cell_kernel_overview(a-1).kernel_type = 'OFF';
                cell_kernel_overview(a-1).kernel_type_log = -1;
                end
            elseif numel(kernel_unique) == 2
                cell_kernel_overview(a-1).luminance_kernel = 0;
                %All kernel values (up to two numbers) that can be devided by 2 are on off, all
                %others are on, cell with numbers with 3 places are complex
                %kernels
                to_test = cell2mat(table_all.kernel_type_log);
                zero_crossing_col = zero_crossings(to_test);
                
                if zero_crossing_col == 1
                    cell_kernel_overview(a-1).kernel_type = 'UV_Blue_OP';
                    cell_kernel_overview(a-1).kernel_type_log = zero_crossing_col;
                    
                elseif zero_crossing_col == 2
                    cell_kernel_overview(a-1).kernel_type = 'Blue_Green_OP';
                    cell_kernel_overview(a-1).kernel_type_log = zero_crossing_col;
                elseif zero_crossing_col == 3
                    cell_kernel_overview(a-1).kernel_type = 'Green_Red_OP';
                    cell_kernel_overview(a-1).kernel_type_log = zero_crossing_col;
                elseif nnz(zero_crossing_col) == 0
                    cell_kernel_overview(a-1).kernel_type = 'Complex_ON_OFF';
                    cell_kernel_overview(a-1).kernel_type_log = zero_crossing_col;
                elseif nnz(zero_crossing_col) >= 2
                    cell_kernel_overview(a-1).kernel_type = 'Complex_opp';
                    cell_kernel_overview(a-1).kernel_type_log = zero_crossing_col;
                    
                end
            elseif numel(kernel_unique) == 3
                cell_kernel_overview(a-1).luminance_kernel = 0;
                to_test = cell2mat(table_all.kernel_type_log);
                zero_crossing_col = zero_crossings(to_test);
                cell_kernel_overview(a-1).kernel_type = 'Complex_ON_OFF';
                cell_kernel_overview(a-1).kernel_type_log = zero_crossing_col;
            end
            cell_kernel_overview(a-1).cell_idx = str2double(Kernel_locations(ii,2));
            %We save the colours as one number here, that makes it easier
            %later to compare different numbers with each other to extract
            %specific kernel types
           
            colour_polarity = cell2mat(kernel_type_log);
          
            
            
            %Here we create a unique colour code for on and off
            %combinations of different colours as follows:
            %           ON  OFF
            %UV         1   2
            %Blue       3   4
            %Green      5   6
            %Red        7   8
            colour_code = zeros(1,length(Idx1));
            
            for aa = 1:length(Idx1)
               colour_idx = Idx1(aa);
               if colour_polarity (colour_idx) == 1
                   colour_code(aa) = colour_idx+colour_idx-1;
               elseif colour_polarity (colour_idx) == -1
                   colour_code(aa) = colour_idx+colour_idx;
               end
               
            end
            
            cell_kernel_overview(a-1).active_colours = str2double(sprintf('%1d',colour_code));
            cell_kernel_overview(a-1).channel_idx = main_channels(a-1);
            cell_kernel_overview(a-1).detailed_info = table_all;
         

    
    
            
        
            
       
            
            
                
                    
                
            
            
         
                        
                    
                
                
               
                

                
            
            

            
            
            
            
            
            
            
        
            
     
    end
        
    
out = sf_organizer(stim_idx,savepath,...
    'variable_name','cell_kernel_overview','variable',cell_kernel_overview); 
     
       
%        Stimulus_info(stim_idx).cell_kernel_overview = cell_kernel_overview;
%        M = matfile(savepath,'Writable',true);
%        M.Stimulus_info = Stimulus_info;
       
    end
    
    
    
    
    
    
