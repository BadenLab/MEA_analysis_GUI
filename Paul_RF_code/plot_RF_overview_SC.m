function out = plot_RF_overview_SC (savepath,add_info)
cell_idx = add_info.kernel_overview_cell;
RF_panel = add_info.settings.kernel_new.RF_overview_panel;
nr_options = 4; %Is there a way to not hard code this?
a = 1;
delete(RF_panel.Children);
%try to find the respective files with the information
try
    S = load(findfile_app(add_info.stim_idx,savepath,"RF_Ident_SS_CI.mat")...
      ,"RF_overview");
    SS_CI_file = S.RF_overview;
    indices = [SS_CI_file.cell_idx];
    index = indices == cell_idx;
    SS_CI_file = SS_CI_file(index);
    SS_CI_file_name = fieldnames(SS_CI_file);
    SS_CI_file_name(end) = [];
    SS_CI_file_name = strrep(SS_CI_file_name,'_',' ');
    SS_CI_file_name = strrep(SS_CI_file_name,'STASD ','');
    SS_CI_file = struct2array(SS_CI_file);
    SS_CI_file(end) = [];
    
    
    %Plotting
    
    ax(a) = subplot(nr_options/2,nr_options/2,a,'parent',RF_panel);
    bar(ax(a),SS_CI_file,'k');
    ax(a).XTickLabel = SS_CI_file_name;
    ax(a).YTickLabel = [];
    title(ax(a),"SS CI")
    xtickangle(ax(a),45)
    
    %Activate checkboxes.
    if any(SS_CI_file)
        set(add_info.settings.kernel_new.panels.ConfIntCheckBox_SS_plot, 'enable', 'on')
        if find(contains(SS_CI_file_name,'Box'))
            set(add_info.settings.kernel_new.panels.BoxCheckBox_plot, 'enable', 'on')
        end

        if find(contains(SS_CI_file_name,'ASP'))
            set(add_info.settings.kernel_new.panels.AllsignificantpixelCheckBox_plot, 'enable', 'on')
        end

        if find(contains(SS_CI_file_name,'Gaus'))
            set(add_info.add_info.settings.kernel_new.panels.GaussianCheckBox_plot, 'enable', 'on')
        end
    end
    
    a = a+1;
catch
    SS_CI_file = NaN;
    a = a+1;
end

try
    S = load(findfile_app(add_info.stim_idx,savepath,"RF_Ident_SS_Std.mat")...
      ,"RF_overview");
    SS_Std_file = S.RF_overview;
    indices = [SS_Std_file.cell_idx];
    index = indices == cell_idx;
    SS_Std_file = SS_Std_file(index);
    SS_Std_file_name = fieldnames(SS_Std_file);
    SS_Std_file_name(end) = [];
    SS_Std_file_name = strrep(SS_Std_file_name,'_',' ');
    SS_Std_file_name = strrep(SS_Std_file_name,'STASD ','');
    SS_Std_file = struct2array(SS_Std_file);
    SS_Std_file(end) = [];
    %Plotting
    
    ax(a) = subplot(nr_options/2,nr_options/2,a,'parent',RF_panel);
    bar(ax(a),SS_Std_file,'k');
    ax(a).XTickLabel = SS_Std_file_name;
    ax(a).YTickLabel = [];
    title(ax(a),"SS SDThres")
    xtickangle(ax(a),45)
    
    if any(SS_Std_file)
        set(add_info.settings.kernel_new.panels.NumStdDevCheckBox_SS_plot, 'enable', 'on')
        if find(contains(SS_Std_file_name,'Box'))
            set(add_info.settings.kernel_new.panels.BoxCheckBox_plot, 'enable', 'on')
        end

        if find(contains(SS_Std_file_name,'ASP'))
            set(add_info.settings.kernel_new.panels.AllsignificantpixelCheckBox_plot, 'enable', 'on')
        end

        if find(contains(SS_Std_file_name,'Gaus'))
            set(add_info.settings.kernel_new.panels.GaussianCheckBox_plot, 'enable', 'on')
        end
    end
    
    
    a = a+1;
catch
    SS_Std_file = NaN;
    a = a+1;
end

try
    S = load(findfile_app(add_info.stim_idx,savepath,"RF_Ident_SC_CI.mat")...
      ,"RF_overview");
    SC_Cl_file = S.RF_overview;
    indices = [SC_Cl_file.cell_idx];
    index = indices == cell_idx;
    SC_Cl_file = SC_Cl_file(index);
    SC_Cl_file_name = fieldnames(SC_Cl_file);
    SC_Cl_file_name(end) = [];
    SC_Cl_file_name = strrep(SC_Cl_file_name,'_',' ');
    SC_Cl_file_name = strrep(SC_Cl_file_name,'SC ','');
    SC_Cl_file = struct2array(SC_Cl_file);
    SC_Cl_file(end) = [];
    
    
    ax(a) = subplot(nr_options/2,nr_options/2,a,'parent',RF_panel);
    bar(ax(a),SC_Cl_file,'k');
    ax(a).XTickLabel = SC_Cl_file_name;
    ax(a).YTickLabel = [];
    title(ax(a),"SC CI")
    xtickangle(ax(a),45)
    
    if any(SC_Cl_file)
        set(add_info.settings.kernel_new.panels.ConfIntCheckBox_SC_plot, 'enable', 'on')
        if find(contains(SC_Cl_file_name,'Box'))
            set(add_info.settings.kernel_new.panels.BoxCheckBox_plot, 'enable', 'on')
        end

        if find(contains(SC_Cl_file_name,'ASP'))
            set(add_info.settings.kernel_new.panels.AllsignificantpixelCheckBox_plot, 'enable', 'on')
        end

        if find(contains(SC_Cl_file_name,'Gaus'))
            set(add_info.add_info.settings.kernel_new.panels.GaussianCheckBox_plot, 'enable', 'on')
        end
    end
    
    a = a+1;
catch
    SC_Cl_file = NaN;
    a = a+1;
end

try
    S = load(findfile_app(add_info.stim_idx,savepath,"RF_Ident_SC_Std.mat")...
      ,"RF_overview");
    SC_Std_file = S.RF_overview;
    indices = [SC_Std_file.cell_idx];
    index = indices == cell_idx;
    SC_Std_file = SC_Std_file(index);
    SC_Std_file_name = fieldnames(SC_Std_file);
    SC_Std_file_name(end) = [];
    SC_Std_file_name = strrep(SC_Std_file_name,'_',' ');
    SC_Std_file_name = strrep(SC_Std_file_name,'SC ','');
    SC_Std_file = struct2array(SC_Std_file);
    SC_Std_file(end) = [];
    
    
    ax(a) = subplot(nr_options/2,nr_options/2,a,'parent',RF_panel);
    bar(ax(a),SC_Std_file,'k');
    ax(a).XTickLabel = SC_Std_file_name;
    ax(a).YTickLabel = [];
    title(ax(a),"SC STD")
    xtickangle(ax(a),45)
    
    if any(SC_Std_file)
        set(add_info.settings.kernel_new.panels.ConfIntCheckBox_SS_plot, 'enable', 'on')
        if find(contains(SC_Std_file_name,'Box'))
            set(add_info.settings.kernel_new.panels.BoxCheckBox_plot, 'enable', 'on')
        end

        if find(contains(SC_Std_file_name,'ASP'))
            set(add_info.settings.kernel_new.panels.AllsignificantpixelCheckBox_plot, 'enable', 'on')
        end

        if find(contains(SC_Std_file_name,'Gaus'))
            set(add_info.add_info.settings.kernel_new.panels.GaussianCheckBox_plot, 'enable', 'on')
        end
    end
        
    a = a+1;
catch
    SC_Std_file = NaN;
    a = a+1;
end


out = 1;





end