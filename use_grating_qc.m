function out = use_grating_qc (savepath,add_info)

%Loading the datafile with the information for the gratings
% We only load the datafiles that are needed based on what the user
% selected
if add_info.settings.Gratings.OTest || add_info.settings.Gratings.RTest
    File = load(findfile_app(add_info.stim_idx,savepath,'Data_circular'));
    Data_circular = File.Data_circular;
    clear File
    nr_cells = length(Data_circular);
end


if add_info.settings.Gratings.QI || add_info.settings.Gratings.FI
    File = load(findfile_app(add_info.stim_idx,savepath,'Grating_QC'));
    Grating_QC = File.Grating_QC;
    clear File
    nr_cells = length(Grating_QC);
end
  


%We have a unknown number of true if statements (depending on the user
%input) so lets use the variable a to count how many of the if statements
%are true. We can also use a as variable to index into an array which keeps
%information about which cells pass which quality criteria, in the end we
%just need to multiply all possible conditions to see which cells pass all
%criteria the user has selected.
a = 0;

true_Gratings = false(nr_cells,1);


%% Quality index
if add_info.settings.Gratings.QI
    
    a = a+1;
    %Check which cells cross the threshold given by the user
    true_Gratings(:,a) = [Grating_QC.QC_pass]>= add_info.settings.Gratings.QI_thr;
  
end


%% Frequency index
if add_info.settings.Gratings.FI
    a = a+1;
    true_Gratings(:,a) = [Grating_QC.freq_pass];
end

%% RTest

if add_info.settings.Gratings.RTest
    a = a+1;
        
    %Checking which cells have passed the quality criteria
    true_Gratings(:,a) = [Data_circular.circ_rtest_sig] == 1;
        
end
% RTest Window
% 
% if add_info.settings.Gratings.RTest_window
%     
%     %Coming soon
%     
% %     %Checking which cells have passed the quality criteria
% %     true_Gratings = [Data_circular.	] == 1;
% %     Gratings_nr = nnz(true_Gratings);
% %     %Put the name of those cells into a cell array
% %     if Gratings_nr ~=  0
% %     RTest_name_str = cell(1,Gratings_nr);   
% %     for ii = 1:Gratings_nr
% %         RTest_name_str{ii} = ['Cell_', num2str(Data_circular(true_Gratings(ii)).cell_idx)];
% %     end
% %     end
%     
%     
%     
%     
% end

%% OTest
if add_info.settings.Gratings.OTest
    a = a+1;
        
    %Checking which cells have passed the quality criteria
    true_Gratings(:,a) = [Data_circular.circ_otest_sig] == 1;
        
end

% OTest window
% 
% if add_info.settings.Gratings.RTest_window
%     
%     %Coming soon
%     
% %     %Checking which cells have passed the quality criteria
% %     true_Gratings = [Data_circular.	] == 1;
% %     Gratings_nr = nnz(true_Gratings);
% %     %Put the name of those cells into a cell array
% %     if Gratings_nr ~=  0
% %     RTest_name_str = cell(1,Gratings_nr);   
% %     for ii = 1:Gratings_nr
% %         RTest_name_str{ii} = ['Cell_', num2str(Data_circular(true_Gratings(ii)).cell_idx)];
% %     end
% %     end
%     
%     
%     
%     
% end
    
        

%% Summary
%Here we check which cells pass all tests and construct their name for the
%Items list

true_Gratings = all(true_Gratings,2);

if nnz(true_Gratings) == 0
    out = {}; %returns empty list
else
    
    cell_idxs = [Grating_QC.cell_idx];
    cell_idxs = cell_idxs(1,true_Gratings);
    out = cell(1,length(cell_idxs));
    for ii = 1:length(cell_idxs)
        out{1,ii} = strcat('Cell_',num2str(cell_idxs(1,ii)));
    end
end

    
       
           
end
                       