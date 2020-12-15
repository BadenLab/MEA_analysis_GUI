classdef MEA_analysis_exported < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        UIFigure                        matlab.ui.Figure
        FileMenu                        matlab.ui.container.Menu
        LoadsortedspikesMenu            matlab.ui.container.Menu
        LoaddatasetMenu                 matlab.ui.container.Menu
        ExportdataMenu                  matlab.ui.container.Menu
        ExportsubsetofspikesMenu        matlab.ui.container.Menu
        ExportanalysisfilesMenu         matlab.ui.container.Menu
        ExportfullrecordingMenu         matlab.ui.container.Menu
        ImportdataMenu                  matlab.ui.container.Menu
        EditMenu                        matlab.ui.container.Menu
        CopypanelcontentMenu            matlab.ui.container.Menu
        SettingsMenu                    matlab.ui.container.Menu
        NoiseFileMenu                   matlab.ui.container.Menu
        KernelMenu                      matlab.ui.container.Menu
        MovieFrequencyMenu              matlab.ui.container.Menu
        ScriptFolderMenu                matlab.ui.container.Menu
        AssignenewMenu                  matlab.ui.container.Menu
        ParallelprocessingMenu          matlab.ui.container.Menu
        OnMenu                          matlab.ui.container.Menu
        OFFMenu                         matlab.ui.container.Menu
        WindowsMenu                     matlab.ui.container.Menu
        CloseallMenu                    matlab.ui.container.Menu
        TabGroup                        matlab.ui.container.TabGroup
        MainAnalysisTab                 matlab.ui.container.Tab
        Mainpagegrid                    matlab.ui.container.GridLayout
        UITable                         matlab.ui.control.Table
        RunButton                       matlab.ui.control.Button
        RefreshListButton               matlab.ui.control.Button
        ChangeFolderButton              matlab.ui.control.Button
        Stim_Axes                       matlab.ui.control.UIAxes
        EditVariableButton              matlab.ui.control.Button
        RemoveselectedfileButton        matlab.ui.control.Button
        CreatenewstimulusfromzoomedButton  matlab.ui.control.Button
        DataListBox                     matlab.ui.control.ListBox
        StimuliListBoxLabel             matlab.ui.control.Label
        StimuliListBox                  matlab.ui.control.ListBox
        ScriptsListBoxLabel             matlab.ui.control.Label
        ScriptsListBox                  matlab.ui.control.ListBox
        AllspikesPanel                  matlab.ui.container.Panel
        QualitycheckPanel               matlab.ui.container.Panel
        DataLabel                       matlab.ui.control.Label
        KernelsTab                      matlab.ui.container.Tab
        UsequalityindexCheckBox         matlab.ui.control.CheckBox
        ShowButton                      matlab.ui.control.Button
        Kernel_Axes                     matlab.ui.control.UIAxes
        Kernel_Axes_2                   matlab.ui.control.UIAxes
        Kernel_Axes_4                   matlab.ui.control.UIAxes
        Kernel_Axes_3                   matlab.ui.control.UIAxes
        TimeSliderLabel                 matlab.ui.control.Label
        TimeSlider                      matlab.ui.control.Slider
        PlayButton                      matlab.ui.control.Button
        StopButton                      matlab.ui.control.Button
        QualityCriteriaPanel            matlab.ui.container.Panel
        ResetButton                     matlab.ui.control.Button
        CalculateButton                 matlab.ui.control.Button
        StdThresholdSliderLabel         matlab.ui.control.Label
        StdThresholdSlider              matlab.ui.control.Slider
        KernelsListBoxLabel             matlab.ui.control.Label
        KernelsListBox                  matlab.ui.control.ListBox
        ShowtimekernelsCheckBox         matlab.ui.control.CheckBox
        NrneighboursSliderLabel         matlab.ui.control.Label
        NrneighboursSlider              matlab.ui.control.Slider
        LinkAxesButton                  matlab.ui.control.Button
        RejectCellQIButton              matlab.ui.control.Button
        KernelsnewTab                   matlab.ui.container.Tab
        Kernel_new_grid                 matlab.ui.container.GridLayout
        RFIdentificationMethodPanel     matlab.ui.container.Panel
        STASDCheckBox                   matlab.ui.control.CheckBox
        SelfCovarianceCheckBox          matlab.ui.control.CheckBox
        RFOptionsPanel                  matlab.ui.container.Panel
        RFTypePanel                     matlab.ui.container.Panel
        BoxCheckBox                     matlab.ui.control.CheckBox
        AllsignificantpixelCheckBox     matlab.ui.control.CheckBox
        GaussianCheckBox                matlab.ui.control.CheckBox
        NumSDGaussSpinnerLabel          matlab.ui.control.Label
        NumSDGaussSpinner               matlab.ui.control.Spinner
        NumRingBoxesSpinnerLabel        matlab.ui.control.Label
        NumRingBoxesSpinner             matlab.ui.control.Spinner
        QualityControlThresholdPanel    matlab.ui.container.Panel
        STASDPanel                      matlab.ui.container.Panel
        QCTypePanel                     matlab.ui.container.Panel
        NumStdDevCheckBox               matlab.ui.control.CheckBox
        ConfIntCheckBox                 matlab.ui.control.CheckBox
        SDThreshEditFieldLabel          matlab.ui.control.Label
        SDThreshEditField               matlab.ui.control.EditField
        CIUpperEditFieldLabel           matlab.ui.control.Label
        CIUpperEditField                matlab.ui.control.EditField
        SelfCovariancePanel             matlab.ui.container.Panel
        QCTypePanel_2                   matlab.ui.container.Panel
        NumStdDevCheckBox_2             matlab.ui.control.CheckBox
        ConfIntCheckBox_2               matlab.ui.control.CheckBox
        SDThreshEditField_2Label        matlab.ui.control.Label
        SDThreshEditField_2             matlab.ui.control.EditField
        CILowerEditFieldLabel           matlab.ui.control.Label
        CILowerEditField                matlab.ui.control.EditField
        StimulusFrequencyEditFieldLabel  matlab.ui.control.Label
        StimulusFrequencyEditField      matlab.ui.control.EditField
        HzLabel                         matlab.ui.control.Label
        FindRFsButton                   matlab.ui.control.Button
        PlottingResultsPanel            matlab.ui.container.Panel
        PlotsinglecellButton            matlab.ui.control.Button
        CelltoplotSpinnerLabel          matlab.ui.control.Label
        CelltoplotSpinner               matlab.ui.control.Spinner
        STASDLabel                      matlab.ui.control.Label
        SelfCovarianceLabel             matlab.ui.control.Label
        NumStdDevCheckBox_SC_plot       matlab.ui.control.CheckBox
        NumStdDevCheckBox_SS_plot       matlab.ui.control.CheckBox
        ConfIntCheckBox_SC_plot         matlab.ui.control.CheckBox
        ConfIntCheckBox_SS_plot         matlab.ui.control.CheckBox
        RFTypePanel_2                   matlab.ui.container.Panel
        BoxCheckBox_plot                matlab.ui.control.CheckBox
        AllsignificantpixelCheckBox_plot  matlab.ui.control.CheckBox
        GaussianCheckBox_plot           matlab.ui.control.CheckBox
        HeatmapButtonGroup              matlab.ui.container.ButtonGroup
        GrayButton                      matlab.ui.control.RadioButton
        ColourButton                    matlab.ui.control.RadioButton
        PlotstatisticsforallcellsButton  matlab.ui.control.Button
        GeneraloptionsPanel             matlab.ui.container.Panel
        NumSTEbinsSpinnerLabel          matlab.ui.control.Label
        NumSTEbinsSpinner               matlab.ui.control.Spinner
        PlotRFinformationCheckBox       matlab.ui.control.CheckBox
        PlotrasterplotPanel             matlab.ui.container.Panel
        CelltoplotSpinner_2Label        matlab.ui.control.Label
        CelltoplotSpinner_2             matlab.ui.control.Spinner
        PlotsinglecellButton_2          matlab.ui.control.Button
        CellRFoverviewPanel             matlab.ui.container.Panel
        RFcellsoverviewPanel            matlab.ui.container.Panel
        GratingsTab                     matlab.ui.container.Tab
        GridLayout                      matlab.ui.container.GridLayout
        PolarPlotsPanel                 matlab.ui.container.Panel
        Panel2                          matlab.ui.container.Panel
        SignificancecriteriaPanel       matlab.ui.container.Panel
        RTestCheckBox                   matlab.ui.control.CheckBox
        OTestCheckBox                   matlab.ui.control.CheckBox
        RTestTimeWindowCheckBox         matlab.ui.control.CheckBox
        OTestTimeWindowCheckBox         matlab.ui.control.CheckBox
        ShowButton_2                    matlab.ui.control.Button
        GratingsListBox                 matlab.ui.control.ListBox
        Grating_cells_qcPanel           matlab.ui.container.Panel
        Grating_qc_summaryPanel         matlab.ui.container.Panel
        Grating_frequencyPanel          matlab.ui.container.Panel
        FFFTab                          matlab.ui.container.Tab
        GridLayout2                     matlab.ui.container.GridLayout
        PlotsPanel                      matlab.ui.container.Panel
        CellListBoxLabel                matlab.ui.control.Label
        CellListBox                     matlab.ui.control.ListBox
        QualityOverviewPanel            matlab.ui.container.Panel
        AlltracesoverviewPanel          matlab.ui.container.Panel
        QualityoverviewdetailedPanel    matlab.ui.container.Panel
        ClusteringTab                   matlab.ui.container.Tab
        StatusLampLabel                 matlab.ui.control.Label
        StatusLamp                      matlab.ui.control.Lamp
        SelectedstimulusEditFieldLabel  matlab.ui.control.Label
        SelectedstimulusEditField       matlab.ui.control.EditField
    end

    
    properties (Access = private)
        
        %%This part declares variables which can be accessed globaly using the app.variable nomenclature
      
        spike_info %Stores information about the main analysis file.
        dataselection %Stores the index of the recording selected in the Listbox containing loaded files
        nr_r %Keeps track of how many files have been loaded in total
        working_path %Working path for the app
        data_info = cell(1,3) %Cell with information shown in the first list box (Filename, cluster, spikes)
        selected_stimuli %Index of the stimulu(i) selected in the stimulus list 
        scriptfolder %= "C:\Users\Marvin Seifert\Documents\MATLAB\App_functions" %Default script folder which gets loaded at startup
        appfolder  % = "C:\Users\Marvin Seifert\Documents\MATLAB\App" %Folder with the main scripts for the app
        selected_script %Stores the name of the script selected in the scripts list box
        to_function %Structure which contains all the information which is input to a respective script once "Run" is clicked
        spiketimestamps %Stores loaded spiketimestamps once run is hit until selected stimulus or recording is changed
        add_info %Stores additional information in a structure and is use as input into the analysis function
        Kernel_info %Structure which saves data for the Kernel tab
        stop_movie %Stores 1 or 0 if button pushed or not
        settings_location %The location of the settings.mat file which is loaded at startup
        out%return of the script run by the run_script function
        stimulus_info %Stores information about the selected stimulus if existing
        Grating_info %Saves data for the Gratings tab like Kernel_info does for Kernels
        Tab_reload %True or false, if true the tabfunction runs when tabs are changed, if false, tabfunction does 
        %not run (Tabfunction shall only run if stimulus selected is changed)
        FFF_figures %Stores handles of figures in the FFF tab
        Cuda_status %Boolean, saves if CUDA driver/ device exists on the computer
        
        all_cells_figures %Stores figures to be potted into the main tab
        
        %Stim Ch
        stim_ch%Data of the stim_ch which is displayed under the stimulus list
        downsample_f%factor by which the stim_ch trace is downsampled to save memory space
        xlineob_1%vertical line which indicates the begining of a stimulus in the stim trace
        xlineob_2%vertical line which indicates the ending of a stimulus in the stim trace
        
        
        test1
        test2
        test3
        test4
             
        %Data Manipulation
        Data_Variable
        
        %Axis Kernel Plot
        k_ax1
        k_ax2
        k_ax3
        k_ax4
        
    end
    
    methods (Access = private)
    %This section contains different functions which can be used within the app (they are not accessable in normal Matlab)    
        
        function run_script(app)
            %This function runs the script with the name defined in the
            %variable (app.)selected_stimuli. It takes only one input which
            %is app but relies on app.selected_stimuli to be present.
             app.StatusLamp.Color = 'red'; %Change lamp to red
             drawnow(); %Update GUI
             
                %In case the script cant be found it returns an error message
            
            try
%                 file = app.data_info{app.dataselection, 1}; 
            selected_stimuli = app.selected_stimuli; %Chronologic number of the selected stimulus
             
            catch ME %Error message returned to the Matlab Command Window
                message = sprintf('Error running script. Select Recording and stimulus first');
                uiwait(warndlg(message));
                app.StatusLamp.Color = 'green';
                return
                
            end
            try 
                %Try update the app.add_info structure which will be handed
                %to the function called. 
                
                app.add_info.stim_begin  = app.spike_info.stim_begin(selected_stimuli);%Only pick the times for the stimuli choosen
                app.add_info.stim_end= app.spike_info.stim_end(selected_stimuli);
                               
            catch ME
                 message = sprintf('Error loading variable "stim_begin" and or "stim_end" \n%s', ME.message);
                 uiwait(warndlg(message));
                 app.StatusLamp.Color = 'green';
                 return
            end
            

            
            try
                function_to_run = erase(app.selected_script,'.m'); %delets the file type description
                %Now we create the call to the function as a string, this
                %includes input arguments, next we will use eval to
                %evaluate that string as if it was written as code
                function_in = strcat(function_to_run,"(app.data_info{app.dataselection,1},app.add_info)"); 
                app.add_info
                app.add_info.out = eval(function_in);
%                 app.add_info.out %uncomment to see what the function
%                 called returns
                
                
            catch ME
                disp( getReport( ME, 'extended', 'hyperlinks', 'on' ) )
                app.StatusLamp.Color = 'green';
                return
                
            end
               
            
            end
            
        
        
        function Tabfunction(app, Tab)
        %This function runs whenever a new tab is selected in the GUI. It
        %is like a startup function for different tabs
        if app.Tab_reload == true %Only run if tab shall be reloaded
            if strcmp(Tab, 'Kernels') == 1
                try
                   
                    
                    %Load the Data file for the selected stimulus
                    S = load(app.spike_info.Stimulus_info{app.add_info.stim_idx});
                    Data = S.Data;
                    clear S
                    
                    %Check if a folder with the name "Kernel" exists
                    Kernel_idx = find([Data.Folder] == "Kernel"); 
                    Kernel_name_str = Data(Kernel_idx).Files;
                    Kernel_length = length(Kernel_name_str);
                    
                    %Sort data based on cell indices
                    TF = zeros(1,Kernel_length);
                    for ii = 1:Kernel_length
                        
                        log = isstrprop(Kernel_name_str{ii},'digit');
                        string_name = Kernel_name_str{ii};
                        TF(ii) = string(string_name(log));
                    end
                    [~, a_order] = sort(TF);
                    Kernel_name_str = Kernel_name_str(a_order);
                    
                    %Update the Kernel ListBox
                    app.KernelsListBox.Items = Kernel_name_str;
                                        
                    
                catch  ME
                   %Errors here mostly mean that no Kernels have been
                   %calculated
                    disp( getReport( ME, 'extended', 'hyperlinks', 'on' ) )
                end
                
                
                %This part updates the FFF tab
            elseif strcmp(Tab, 'FFF') == 1
                try  %Since FFF are much more straightforward, we just list all indices of all cells in the file
                     File = load(app.data_info{app.dataselection,1},'-mat','all_clusters','cell_indices'); 
                     Cell_nr = File.all_clusters;
                     Cell_nr_str = File.cell_indices;
                     Cell_name_str = {1,Cell_nr};
                     for ii = 1:Cell_nr
                         Cell_name_str{1, ii} = strcat('Cell_',num2str(Cell_nr_str(ii)));
                     end
                     %update the Cell Listbox
                     app.CellListBox.Items = Cell_name_str;
                 catch  ME
                   %For debugging
                    disp( getReport( ME, 'extended', 'hyperlinks', 'on' ) )
                end
                
                %Delet existing figures
                try
                    delete(app.QualityOverviewPanel.Children);
                catch ME
                    disp( getReport( ME, 'extended', 'hyperlinks', 'on' ) )
                end
                try
                    delete(app.AlltracesoverviewPanel.Children);
                catch ME
                    disp( getReport( ME, 'extended', 'hyperlinks', 'on' ) )
                end
                try
                    delete(app.QualityoverviewdetailedPanel.Children);
                catch ME
                    disp( getReport( ME, 'extended', 'hyperlinks', 'on' ) )
                end
                
                %Try to load figures
                try
                    dummy_figure = openfig(findfile_app(app.add_info.stim_idx,app.spike_info.savename,"overview_FFF.fig"),'invisible');
                    dummy_ax = findobj(dummy_figure,'Type','Axes');
                    app.FFF_figures.Overview = copyobj(dummy_ax, app.QualityOverviewPanel);    
                    delete(dummy_figure)
                    1
                    
                                        
                catch ME
                    disp( getReport( ME, 'extended', 'hyperlinks', 'on' ) )
                end
                
                try
                    dummy_figure = openfig(findfile_app(app.add_info.stim_idx,app.spike_info.savename,"all_traces_plot.fig"),'invisible');
                    dummy_ax = findobj(dummy_figure,'Type','Axes');
                    app.FFF_figures.Alltraces = copyobj(dummy_ax, app.AlltracesoverviewPanel);    
                    delete(dummy_figure)
                    colormap(app.FFF_figures.Alltraces,flipud(gray));
                                      
                catch ME
                    disp( getReport( ME, 'extended', 'hyperlinks', 'on' ) )
                end
                
                try
                    dummy_figure = openfig(findfile_app(app.add_info.stim_idx,app.spike_info.savename,"detail_stats_FFF.fig"),'invisible');
                    dummy_ax = findobj(dummy_figure,'Type','Axes');
                    colormap (flipud(gray));
                    app.FFF_figures.QualityOV = copyobj(dummy_ax, app.QualityoverviewdetailedPanel);    
                    delete(dummy_figure)
                    colormap(app.FFF_figures.QualityOV,flipud(gray));              
                catch ME
                    disp( getReport( ME, 'extended', 'hyperlinks', 'on' ) )
                end
                
                
                
                
                
                
               %This part updates the Gartings Tab, equal to the FFF tab
            elseif strcmp(Tab, 'Gratings') == 1
                if length(app.add_info.stim_idx) > 1
                    warning("Multiple stimuli selected, this tab only works with one stimulus selected")
                    return
                end
                
                try
                    File = load(app.data_info{app.dataselection,1},'-mat','all_clusters','cell_indices'); 
                    Cell_nr = File.all_clusters;
                    Cell_nr_str = File.cell_indices;
                    Cell_name_str = {1,Cell_nr};
                    for ii = 1:Cell_nr
                        Cell_name_str{1, ii} = strcat('Cell_',num2str(Cell_nr_str(ii)));
                    end
                    app.GratingsListBox.Items = Cell_name_str;
                    
                                       
                    
                catch ME
                     disp( getReport( ME, 'extended', 'hyperlinks', 'on' ) )
                end
                
                try
                    %Try to load figures with the QC information
                    dummy_figure = openfig(findfile_app(app.add_info.stim_idx,app.spike_info.savename,"Grating_qc_figure.fig"),'invisible');
                    dummy_ax = findobj(dummy_figure,'Type','Axes');
                    copyobj(dummy_ax, app.Grating_cells_qcPanel);    
                    delete(dummy_figure)
                    
                catch ME
                     disp( getReport( ME, 'extended', 'hyperlinks', 'on' ) )
                end
                
                 try
                    %Try to load figures with the QC information
                    dummy_figure = openfig(findfile_app(app.add_info.stim_idx,app.spike_info.savename,"grating_qc_summary.fig"),'invisible');
                    dummy_ax = findobj(dummy_figure,'Type','Axes');
                    copyobj(dummy_ax, app.Grating_qc_summaryPanel);    
                    delete(dummy_figure)
                    
                catch ME
                     disp( getReport( ME, 'extended', 'hyperlinks', 'on' ) )
                 end
                 
                 try
                    %Try to load figures with the QC information
                    dummy_figure = openfig(findfile_app(app.add_info.stim_idx,app.spike_info.savename,"frequency_qc_summary.fig"),'invisible');
                    dummy_ax = findobj(dummy_figure,'Type','Axes');
                    copyobj(dummy_ax, app.Grating_frequencyPanel);    
                    delete(dummy_figure)
                    
                catch ME
                     disp( getReport( ME, 'extended', 'hyperlinks', 'on' ) )
                 end
                 
                 
                 
                 
            elseif strcmp(Tab, 'Kernels new') == 1
                
                try
                    app.add_info.settings.kernel_new.cells_overview_panel = app.RFcellsoverviewPanel; 
                    app.selected_script = "plotoverview_kernel_new.m";
                    run_script(app)
                catch ME
                    disp( getReport( ME, 'extended', 'hyperlinks', 'on' ) )
                end
                
                
            end
            
        end
        app.Tab_reload = false;
        end
        
        
        
        function screen_position(app)
            %This function sets the screen position of the app during
            %startup
             screenSize = get(groot,'ScreenSize');
             screenWidth = screenSize(3);
             screenHeight = screenSize(4);
             left = screenWidth*0.1;
             bottom = screenHeight*0.1;
             width = screenWidth*0.8;
             height = screenHeight*0.8;
             drawnow;
             app.UIFigure.Position = [left bottom width height]; 
          
        end
        
        function edit_Variable(btn,app,uit,name)
            %This function saves a new value assigned to a variable in the
            %varaible list box when edit variable and in "save" is clicked
            %in the pop up window
            new_Variable = uit.Data;
            MF = matfile(app.data_info{app.dataselection,1},"Writable",true);
            MF.(name) = new_Variable;
            
            
        end
        
        
        
        function savesettings(app) %Saves new settings to settings file after updated in the GUI
            
            settings = app.add_info.settings;
            save(app.settings_location,"settings",'-mat')
        end
    end
    

    % Callbacks that handle component events
    methods (Access = private)

        % Code that executes after component creation
        function startupFcn(app)
            %%This function assignes values to variables declared earlier
            %%and runs functions before the actual app is started
            
            %Witch lamp to red
            app.StatusLamp.Color = 'red';
            drawnow() %Updates the gui
            %
            
            %Set the resolution of the app
            screen_position(app)
            
            % Define main path and relative paths from settings
            [app.appfolder,~,~] =  fileparts(mfilename('fullpath')); %Get the folder in which the app is running
            
            app.settings_location = strcat(app.appfolder,'\Settings.mat'); %Sets the location of the file with the settings information
            
            % Settings
            %This loads the settings from the settings file or creates a
            %new settings file when no existing file is found
            try
                S = load(app.settings_location,'-mat'); %Load the data from the file specified above
                app.add_info.settings = S.settings;
                app.scriptfolder = S.settings.location.scripts;
                
            catch
                disp('Could not find settings file, settings set to default')
                path = uigetdir(app.appfolder,'Select main script folder');
                if isequal(path,0)
                   disp('User selected Cancel')
                else
                   disp(['User selected ', path])
                   app.add_info.settings.location.scripts = path;
                   app.scriptfolder = path;
                   addpath(genpath(app.appfolder)); %This adds the app folder to the matlab path
                end
                savesettings(app)
                
            end
           
            %General
            app.dataselection = uint8(0);
            app.UITable.ColumnName = ["Filename", "Clusters", "Spikes"]; %Set default column names for the datafile table
            app.UITable.ColumnFormat = {'char', 'numeric', 'numeric'}; %Define the formate of the columns
            app.UITable.ColumnEditable = true(1,3); %all columns editable
            app.nr_r = uint8(1);
            app.StimuliListBox.Multiselect = 'on'; %This enables multiselections in the stimulus list
            app.StimuliListBox.Items = {};
            app.Cuda_status = logical(gpuDeviceCount);
            app.appfolder = pwd;
            %addpath(app.scriptfolder); %Adds the scriptfolder to the matlab search path
            addpath(app.appfolder); %Adds the folder with the main functions for the app
            app.UITable.Enable = 'off';
            
            %Set the script List
            folder_content = dir(app.scriptfolder); %List all contents in the script folder
            app.ScriptsListBox.Items = {folder_content(3:end).name}; %Adds the content into the list in the GUI
            app.ScriptsListBox.Multiselect = 'off'; %Disables multiselection in the list
            app.selected_script = string(app.ScriptsListBox.Value); %Save the name of the first script in the list at startup
            
            %Tabs general
            app.Tab_reload = true;
            
            %Kernel Tab
            app.Kernel_info.Cell = 0; %Sets Kernel count to zero
            app.TimeSlider.Limits = [-100 100];
            app.Kernel_info.kernel_time = 250;
            app.stop_movie = 0;
            app.Kernel_info.std_threshold = 2;
            app.Kernel_info.neighbours = 1; %Determinds how many time series are shown
            app.Kernel_info.time_kernels = 0; %Stores value if time kernels shall be plotted or not
            app.Kernel_info.movie_frequency = 0.05;
            app.ShowtimekernelsCheckBox.Enable = 'off';
            %Stim_Ch
            app.downsample_f = 10;
            app.Stim_Axes
                        
            %Gratings tab
            app.Grating_info.rtest = 0;
            app.Grating_info.otest = 0;
                       
            app.PolarPlotsPanel.AutoResizeChildren = 'off';
            app.Panel2.AutoResizeChildren = 'off';
            app.Grating_info.rtest = 0;
            app.Grating_info.otest = 0;%Stores inital parameters for Quality Criteria Tick Boxes in Gratings Tab
%             app.test1 = subplot(1,4,1,'Parent',app.GratingsTab);
%             app.test2 = subplot(1,4,2,'Parent',app.GratingsTab);
%             app.test3 = subplot(1,4,3,'Parent',app.GratingsTab);
%             app.test4 = subplot(1,4,4,'Parent',app.GratingsTab);
%             plot(app.test1,[1 2 3 4 5])
%             plot(app.test2,[9 8 7 6 5])
%             plot(app.test3,[9 8 7 6 5])
%             plot(app.test4,[9 8 7 6 5])
%                   

            %FFF tab
            app.PlotsPanel.AutoResizeChildren = 'off';
            
            %Kernel_new tab
            set(findall(app.STASDPanel, '-property', 'enable'), 'enable', 'off')     
            set(findall(app.SelfCovariancePanel, '-property', 'enable'), 'enable', 'off')
            set(findall(app.GaussianCheckBox, '-property', 'enable'), 'enable', 'off')
            app.CellRFoverviewPanel.AutoResizeChildren = 'off';
            app.RFcellsoverviewPanel.AutoResizeChildren = 'off';

            
            
            app.StatusLamp.Color = 'green';
            drawnow() %Updates the gui
        end

        % Menu selected function: LoadsortedspikesMenu
        function LoadsortedspikesMenuSelected(app, event)
            %This function runs when "Load sorted spikes" in the menu is clicked 
            %It runs the function "Loadspikes_app" which is stored in the
            %apps main folder. Once this function has run the Uitable is
            %filled with the information about the cells and spikes that
            %have been loaded
            
            app.StatusLamp.Color = 'red';
            drawnow();
                        
            app.Tab_reload = true; %If true tabs will return to default look          
            try
                app.spike_info = Loadspikes_app; %The out file in loadspikes_app must contain all the information which is saved in the analysis mat file
                               
                %Write the essential information into the apps memory
                app.working_path(app.nr_r) = string(app.spike_info.pathname);
                app.data_info{app.nr_r, 1} = app.spike_info.savename;
                app.data_info{app.nr_r, 2} = app.spike_info.all_clusters;
                app.data_info{app.nr_r, 3} = app.spike_info.all_spikes;
                
                %Update Table in the GUI (Table is switched off as long as
                %no recording is selected
                if app.nr_r == 1 
                    app.UITable.Enable = 'on';
                end
                app.UITable.Data = app.data_info;
                
                %Update number of recordings
                app.nr_r = app.nr_r + 1;
            catch ME %Error
                    message = sprintf('Error loading spiketimes\n%s', ME.message);
                    uiwait(warndlg(message));
            end
            
            
            app.StatusLamp.Color = 'green';
            drawnow();
            figure(app.UIFigure);
        end

        % Menu selected function: LoaddatasetMenu
        function LoaddatasetMenuSelected(app, event)
            %This function loads an existing dataset previously created and
            %saved by the app. It than adds the information into the
            %Uitable. It returns an error if the data cant be loaded, or essential information 
            %has not been saved into the savefile
            
            app.StatusLamp.Color = 'red';
            drawnow();
            
            app.Tab_reload = true; %If true tabs will return to default look                 
            try
                
                [hdfname,pathname] = uigetfile('*.mat','Select .mat data file'); %Ask user to choose a file to load into the app
                data_file=strcat(pathname,hdfname);
                data = load(data_file,'-mat','pathname', 'savename', 'all_spikes', 'all_clusters'); %loads the file
                
                %Check if file is already loaded
                Idx = strfind(app.data_info{:,1}, data.savename);
                
                if Idx == 1
                    error('File has already been loaded')
                end
                
                app.working_path(app.nr_r) = string(data.pathname); %updates the essential information in the app
                app.data_info{app.nr_r, 1} = data.savename;
                app.data_info{app.nr_r, 2} = data.all_clusters;
                app.data_info{app.nr_r, 3} = data.all_spikes;
%                 app.savenames(app.nr_r) = data.savename;
%                 app.spike_nr(app.nr_r) = data.all_spikes;
%                 app.cluster_nr(app.nr_r) = data.all_clusters;
                app.UITable.Data = app.data_info;
                app.spike_info = load(data_file,'-mat');
                
                if app.nr_r == 1 
                    app.UITable.Enable = 'on'; %Update Table in the GUI (Table is switched off as long as
                %no recording is selected
                end
                %Update number of recordings
                app.nr_r = app.nr_r+1;
            catch ME %Error handling
                    message = sprintf('Error loading file \n%s', ME.message);
                    uiwait(warndlg(message));
            end
            
                    
                       
            
            app.StatusLamp.Color = 'green';
            drawnow();
            figure(app.UIFigure);
        end

        % Cell selection callback: UITable
        function UITableCellSelection(app, event)
            %This is the callback function for selection of a dataset in
            %the UItable
                       
            %%Update the GUI depending on the user selection in the table
            indices = event.Indices;
            app.dataselection = indices(1,1);
            
            %update spike info, which contains all information about the
            %selected recording
            app.spike_info = load(app.data_info{app.dataselection,1},'-mat'); 
            %Items = load(app.data_info{app.dataselection,1},'-mat','stim_list');%Load stim list from the file selected
            app.StimuliListBox.Items = app.spike_info.stim_list; %Change stim list in the GUI
            app.StimuliListBox.ItemsData = 1:numel(app.spike_info.stim_list); %This gives an index to every stimulus item in the list so that the stimulus can be picked later
            %Change the contents of the variables window (depending on
            %which file was selected)
            listOfVariables = fieldnames(app.spike_info);
            app.DataListBox.Items = listOfVariables;
                        
            %%Update the stimulus trace figure
            %S = load(app.data_info{app.dataselection,1},'-mat','Ch'); %Load the stimulus channel from the mat file
            if app.Cuda_status
                app.stim_ch = gpuArray(app.spike_info.Ch.Ch01_02);%Load stimulus channel into gpu
                xvalue = (0:app.downsample_f:numel(app.stim_ch));%create values for x-axis
                app.stim_ch = downsample(app.stim_ch,app.downsample_f);%Trace is downsampled for performance reasons
                %sampling_Freq = S.Ch.SamplingFrequency;
            
            else
                
                app.stim_ch = app.spike_info.Ch.Ch01_02;
                Ch_length = length(app.stim_ch);
                xvalue = (0:app.downsample_f:Ch_length);%create values for x-axis
                ds_array = xvalue+1;
                if length(ds_array)>Ch_length
                    ds_array = ds_array(1,1:Ch_length);
                end
                try
                    app.stim_ch = app.stim_ch(1,ds_array);
                catch ME
                    ds_array = ds_array(1,1:end-1);
                    app.stim_ch = app.stim_ch(1,ds_array);
                end
               
                
            end
                if size(xvalue,2) > size(app.stim_ch)
                    difference = size(xvalue,2) - size(app.stim_ch,2);
                    idx_begin = size(xvalue,2)-difference;
                    xvalue(:,idx_begin+1:end) = [];
                end
                plot(app.Stim_Axes,xvalue,app.stim_ch);%Plot the trace to the GUI
            
                
                
             %Load figures into Main Tab
            try
                delete(app.AllspikesPanel.Children);
            catch ME
                disp( getReport( ME, 'extended', 'hyperlinks', 'on' ) )
            end
            
            %Try to load figures
            try
                dummy_figure = openfig([app.spike_info.pathname,'hdf_figure.fig'],'invisible');
                dummy_ax = findobj(dummy_figure,'Type','Axes');
                app.all_cells_figures.allspikes = copyobj(dummy_ax, app.AllspikesPanel);    
                delete(dummy_figure)
                1
                
                                    
            catch ME
                disp( getReport( ME, 'extended', 'hyperlinks', 'on' ) )
            end
            
            
            try
                delete(app.QualitycheckPanel.Children);
            catch ME
                disp( getReport( ME, 'extended', 'hyperlinks', 'on' ) )
            end
            
             %Try to load figures
            try
                dummy_figure = openfig([app.spike_info.pathname,'freq_figure.fig'],'invisible');
                dummy_ax = findobj(dummy_figure,'Type','Axes');
                app.all_cells_figures.freq = copyobj(dummy_ax, app.QualitycheckPanel);    
                delete(dummy_figure)
                1
                
                                    
            catch ME
                disp( getReport( ME, 'extended', 'hyperlinks', 'on' ) )
            end
            
                      
           
            %app.dataselection
                        
        end

        % Value changed function: StimuliListBox
        function StimuliListBoxValueChanged(app, event)
            %This is the callback function in case one or more stimuli are selected
            %in the stimulus listbox
            app.add_info.stim_idx = event.Value; %Get the index
            app.SelectedstimulusEditField.Value = app.spike_info.stim_list{app.add_info.stim_idx};
            app.Tab_reload = true; %Set tabs to default
            %Loading the stimulus channel and the stim_begin and stim_end
            %times for the selected stimuli
            app.selected_stimuli = app.StimuliListBox.Value;
            
            try
            S = load(app.data_info{app.dataselection,1},'-mat','Stimulus_info');
            app.stimulus_info = S.Stimulus_info(app.selected_stimuli);
            clear S
            catch ME %Error handling
                %message = sprintf('Error loading file \n%s', ME.message);
                %uiwait(warndlg(message));
            end
                       
            s_s = app.selected_stimuli;
            nr_select = numel(s_s);
            
            
            %Update the stim channel figure
            sampling_Freq = app.spike_info.Ch.SamplingFrequency;
                        
                       
            %Display the stim_begin and stim_end times in the Stim_channel
            %figure
            try
                %This erases existing lines
                 set(app.xlineob_1,'visible','off')
                 set(app.xlineob_2,'visible','off')
            catch ME
                
            end
           %Create new lines
            for ii = 1:nr_select
                stim_begin = app.spike_info.stim_begin(s_s(ii))*sampling_Freq;
                stim_end = app.spike_info.stim_end(s_s(ii))*sampling_Freq;
                
                try
                    app.xlineob_1(ii) = xline(app.Stim_Axes,stim_begin,'g',{num2str(stim_begin)});
                    app.xlineob_2(ii) = xline(app.Stim_Axes,stim_end,'r',{num2str(stim_end)});
                catch 
                    try
                        app.xlineob_2(ii) = xline(app.Stim_Axes,stim_end,'r',{num2str(stim_end)});
                    catch 
                        continue
                    end
                    continue
                end
            end
            
            
                       
                        
        end

        % Button pushed function: ChangeFolderButton
        function ChangeFolderButtonPushed(app, event)
            %This function sets a new folder for the scripts that can be
            %run
            app.scriptfolder = uigetdir;
            folder_content = dir(app.scriptfolder); %List all contents in the script folder
            app.ScriptsListBox.Items = {folder_content(3:end).name}; %Adds the content into the list in the GUI
        end

        % Button pushed function: RunButton
        function RunButtonPushed(app, event)
                        
            %This runs the activated script and hands over the respective
            %stimulus information
            app.StatusLamp.Color = 'red';
            app.selected_script = app.ScriptsListBox.Value;
            app.selected_script
            run_script(app)
            app.StatusLamp.Color = 'green'; 
        end

        % Value changed function: ScriptsListBox
        function ScriptsListBoxValueChanged(app, event)
            app.selected_script = string(app.ScriptsListBox.Value); %Save the function name selected as string
        end

        % Button pushed function: RefreshListButton
        function RefreshListButtonPushed(app, event)
            folder_content = dir(app.scriptfolder); %List all contents in the script folder
            app.ScriptsListBox.Items = {folder_content(3:end).name}; %Adds the content into the list in the GUI
        end

        % Button pushed function: ShowButton
        function ShowButtonPushed(app, event)
            %This function shows spatial kernels calculated before
          
            
            
          cla(app.Kernel_Axes,'reset')
          cla(app.Kernel_Axes_2,'reset')
          cla(app.Kernel_Axes_3,'reset')
          cla(app.Kernel_Axes_4,'reset')
            
            app.StatusLamp.Color = 'red';
            %app.Kernel_info
        
            
            try
                app.add_info.Kernel_info = app.Kernel_info;
                app.selected_script = "Get_spatial_kernel_app.m";
                run_script(app);
                app.Kernel_info.Kernels = app.add_info.out;
                
                %%Plot kernels
                
                app.TimeSlider.Limits = [-app.Kernel_info.Kernels.bins/2 app.Kernel_info.Kernels.bins/2];
                %Creates an imshow from the kernel file
                imshow(app.Kernel_info.Kernels.uv_kernel_grey(:,:,app.Kernel_info.kernel_time),'parent',app.Kernel_Axes);
                colormap(app.Kernel_Axes,Colormap.dense);
%                 app.Kernel_Axes.Position(3) = 250;
%                 app.Kernel_Axes.Position(4) = 250; 
                imshow(app.Kernel_info.Kernels.blue_kernel_grey(:,:,app.Kernel_info.kernel_time),'parent',app.Kernel_Axes_2);
                colormap(app.Kernel_Axes_2,flipud(Colormap.ice));
%                 app.Kernel_Axes_2.Position(3) = 250;
%                 app.Kernel_Axes_2.Position(4) = 250; 
                imshow(app.Kernel_info.Kernels.green_kernel_grey(:,:,app.Kernel_info.kernel_time),'parent',app.Kernel_Axes_3);
                colormap(app.Kernel_Axes_3,Colormap.algae);
%                 app.Kernel_Axes_3.Position(3) = 250;
%                 app.Kernel_Axes_3.Position(4) = 250; 
                imshow(app.Kernel_info.Kernels.red_kernel_grey(:,:,app.Kernel_info.kernel_time),'parent',app.Kernel_Axes_4);
                colormap(app.Kernel_Axes_4,Colormap.amp);
%                 app.Kernel_Axes_4.Position(3) = 250;
%                 app.Kernel_Axes_4.Position(4) = 250;

                %% Get the limits of the axes
                
                
                app.k_ax1.xlim = app.Kernel_Axes.XLim;
                app.k_ax1.ylim = app.Kernel_Axes.YLim;
                
                app.k_ax2.xlim = app.Kernel_Axes_2.XLim;
                app.k_ax2.ylim = app.Kernel_Axes_2.YLim;
                
                app.k_ax3.xlim = app.Kernel_Axes_3.XLim;
                app.k_ax3.ylim = app.Kernel_Axes_3.YLim;
                
                app.k_ax4.xlim = app.Kernel_Axes_4.XLim;
                app.k_ax4.ylim = app.Kernel_Axes_4.YLim;
                
                
            catch ME
                message = sprintf('Error loading Kernel \n%s', ME.message);
                uiwait(warndlg(message));
            end
            
                        
            %Check if time kernels shall be plotted
            if app.Kernel_info.time_kernels == 1
                try
                    app.add_info.Kernel_info = app.Kernel_info;
                    app.selected_script = "Plot_time_kernels.m";
                    run_script(app);
                    
                catch ME
                    message = sprintf('Couldnt plot spatial kernels \n%s', ME.message);
                    uiwait(warndlg(message));
                end
                app.Kernel_info.time_kernels = 0;
            end
                
            app.StatusLamp.Color = 'green';
        end

        % Value changed function: UsequalityindexCheckBox
        function UsequalityindexCheckBoxValueChanged(app, event)
            %This function erases all cells which have no passed the
            %kernels quality criteria from the kernel list box
            Box_tick = app.UsequalityindexCheckBox.Value;
            
            if Box_tick == 0
                app.Tab_reload = true;           
                Tabfunction(app,"Kernels")
                app.ShowtimekernelsCheckBox.Value = 0;
                app.ShowtimekernelsCheckBox.Enable = 'off';
            else
                app.ShowtimekernelsCheckBox.Enable = 'on';
                %Load Kernel location file                
                File = load(findfile_app(app.selected_stimuli,app.spike_info.savename,'Kernel_info'));
                Kernel_info = File.Kernel_info;
                
                Kernel_nr = numel(Kernel_info);
                Kernel_name_str = {};
                for ii = 1:Kernel_nr
                    true_Kernel_log(ii) = logical(nnz(Kernel_info(ii).true_kernel_log));
                end
                
                S = load(findfile_app(app.add_info.stim_idx,app.spike_info.savename,'Kernel_location'));
                kernel_location = S.Kernel_location(true_Kernel_log,1);
                
                %Create short name
                for ii = 1:length(kernel_location)
                    cell_string = char(kernel_location(ii,1));
                    slash_cut = strfind(cell_string,'\');
                    slash_cut = slash_cut(end)+1;
                   
                    Kernel_name_str{1,ii} = cell_string(slash_cut:end);
                                        
                end
               
                
%                 S = load(app.spike_info.Stimulus_info{app.add_info.stim_idx});
%                 Data = S.Data;
%                 clear S
%                                
%                 Kernel_idx = find([Data.Folder] == "Kernel");
%                 Kernel_name_str = Data(Kernel_idx).Files;
%                 Kernel_name_str = {Kernel_name_str{1,kernel_cells}};
%                 
%                 %Sort data
%                 Kernel_length = length(Kernel_name_str);
%                 %Sort data
%                 TF = zeros(1,Kernel_length);
%                 for ii = 1:Kernel_length
%                     
%                     log = isstrprop(Kernel_name_str{ii},'digit')
%                     string_name = Kernel_name_str{ii};
%                     TF(ii) = string(string_name(log));
%                 end
%                 
%                 [~, a_order] = sort(TF);
%                 Kernel_name_str = Kernel_name_str(a_order);
%                 
%                 
%                 
% %                 true_Kernel_idx = find(true_Kernel_log);
% %                 a = 1;
% %                 for ii = true_Kernel_idx
% %                     Kernel_name_str{1, a} = strcat('Cell_',num2str(ii));
% %                     a = a+1;                         
% %                 end
                app.KernelsListBox.Items = Kernel_name_str;
                app.KernelsListBox.Value = Kernel_name_str{1,1};
                app.Kernel_info.Cell = str2double(regexp(app.KernelsListBox.Value,'\d*','Match'));
            end
            
            
        end

        % Value changed function: TimeSlider
        function TimeSliderValueChanged(app, event)
            %This function changes the spatial kernel imshow to the time
            %specified in the time slider
            
            app.Kernel_info.kernel_time = ceil(app.TimeSlider.Value + app.Kernel_info.Kernels.bins/2);
            app.TimeSlider.Limits = [-app.Kernel_info.Kernels.bins/2 app.Kernel_info.Kernels.bins/2];
            
%             %Load colourmaps
%             %colormap hot
%             c(:,:,1) = colormap(Colormap.amp);
%             %colormap summer
%             c(:,:,2) = colormap(Colormap.algae);
%             %colormap winter
%             c(:,:,3) = colormap(Colormap.ice);
%             c(:,:,3) = flipud(c(:,:,3));
%             %colormap cool
%             c(:,:,4) = colormap(Colormap.dense);
                
            imshow(app.Kernel_info.Kernels.uv_kernel_grey(:,:,app.Kernel_info.kernel_time),'parent',app.Kernel_Axes);
            colormap(app.Kernel_Axes,Colormap.dense);
            imshow(app.Kernel_info.Kernels.blue_kernel_grey(:,:,app.Kernel_info.kernel_time),'parent',app.Kernel_Axes_2);
            colormap(app.Kernel_Axes_2,flipud(Colormap.ice));
            imshow(app.Kernel_info.Kernels.green_kernel_grey(:,:,app.Kernel_info.kernel_time),'parent',app.Kernel_Axes_3);
            colormap(app.Kernel_Axes_3,Colormap.algae);
            imshow(app.Kernel_info.Kernels.red_kernel_grey(:,:,app.Kernel_info.kernel_time),'parent',app.Kernel_Axes_4);
            colormap(app.Kernel_Axes_4,Colormap.amp);
            
            app.k_ax1.xlim = app.Kernel_Axes.XLim;
            app.k_ax1.ylim = app.Kernel_Axes.YLim;
            
            app.k_ax2.xlim = app.Kernel_Axes_2.XLim;
            app.k_ax2.ylim = app.Kernel_Axes_2.YLim;
            
            app.k_ax3.xlim = app.Kernel_Axes_3.XLim;
            app.k_ax3.ylim = app.Kernel_Axes_3.YLim;
            
            app.k_ax4.xlim = app.Kernel_Axes_4.XLim;
            app.k_ax4.ylim = app.Kernel_Axes_4.YLim;
            
        end

        % Button pushed function: PlayButton
        function PlayButtonPushed(app, event)
            %This function runs the imshow like a movie by updating the
            %timeslider value in steps according to the frequency specified
            
            %app.TimeSlider.Value = (-app.Kernel_info.Kernels.bins/2)+1;
            %this sets the value to the begining each time play is pressed
           
            for ii = 1:app.Kernel_info.Kernels.bins
                TimeSliderValueChanged(app)
                app.TimeSlider.Value = app.TimeSlider.Value+1;
                pause(app.Kernel_info.movie_frequency)
                if app.stop_movie == 1
                    app.stop_movie = 0;
                    break
                end
            end
        end

        % Button pushed function: StopButton
        function StopButtonPushed(app, event)
            app.stop_movie = 1; %Stops Kernel movie 
        end

        % Value changed function: KernelsListBox
        function KernelsListBoxValueChanged(app, event)
           app.Kernel_info.previous = event.PreviousValue;
           Name = app.KernelsListBox.Value;
           app.Kernel_info.Cell = str2double(regexp(Name,'\d*','Match'));
           
           if app.ShowtimekernelsCheckBox.Value == 1
           app.Kernel_info.time_kernels = 1;
           end
           ShowButtonPushed(app, event)
          
        end

        % Selection change function: TabGroup
        function TabGroupSelectionChanged(app, event)
            selectedTab = app.TabGroup.SelectedTab;
            Tabfunction(app, selectedTab.Title);
            
        end

        % Value changed function: NrneighboursSlider
        function NrneighboursSliderValueChanged(app, event)
            app.Kernel_info.neighbours = app.NrneighboursSlider.Value;
            
        end

        % Value changed function: ShowtimekernelsCheckBox
        function ShowtimekernelsCheckBoxValueChanged(app, event)
            app.Kernel_info.time_kernels = app.ShowtimekernelsCheckBox.Value;
            
        end

        % Button pushed function: CalculateButton
        function CalculateButtonPushed(app, event)
            try
                
                app.add_info.Kernel_info =  app.Kernel_info;
                app.selected_script = "kernel_quality_std_app.m";
                run_script(app);
            catch ME
                message = sprintf('Error loading Kernel \n%s', ME.message);
                uiwait(warndlg(message));
                
            end
                
        end

        % Menu selected function: NoiseFileMenu
        function NoiseFileMenuSelected(app, event)
            %Choose noise file for the kernel analysis, location is saved
            %in settings, so available also after app restarts
            
            app.StatusLamp.Color = 'red';
            drawnow();
            
            try
                [file,path] = uigetfile({'*.h5';'*.hdf5'});
                app.add_info.settings.location.noise = strcat(path,file);
                savesettings(app);
            catch ME
                message = sprintf('Error changing noise file \n%s', ME.message);
                uiwait(warndlg(message));
            end
            
            app.StatusLamp.Color = 'green';
            drawnow();
                        
        end

        % Value changed function: DataListBox
        function DataListBoxValueChanged(app, event)
            app.Data_Variable = app.DataListBox.Value;
            
        end

        % Button pushed function: EditVariableButton
        function EditVariableButtonPushed(app, event)
            %Try to get the name of the selected variable
            try 
                name = string(app.Data_Variable);
            catch ME
                message = sprintf('No Variable selected \n%s', ME.message);
                uiwait(warndlg(message));
            end
            %Load the data from that variable
            app.data_info{app.dataselection,1}
            S = load(app.data_info{app.dataselection,1},'-mat',name);
            eval_string = strcat('S.',name);
            variable = eval(eval_string);
            fig = uifigure('Position',[100 100 350 275]);
            uit = uitable(fig,'Position',[1 100 350 100],'Data',variable);
            set(uit,'ColumnEditable',true(1,size(variable,2)))
            btn = uibutton(fig,'push','Text','Save','Position',[1,250,50,25],'ButtonPushedFcn',@(btn,event)app.edit_Variable(app,uit,name));
            
            
        end

        % Button pushed function: RemoveselectedfileButton
        function RemoveselectedfileButtonPushed(app, event)
            %Get data from the table
            Data = app.UITable.Data;
            Data(app.dataselection,:) = [];
            app.UITable.Data = Data;
            app.dataselection = app.dataselection -1;
            app.nr_r = app.nr_r-1;
            try
            %Reset all the tables:
            app.StimuliListBox.Items = {}; %Change stim list in the GUI
            app.StimuliListBox.ItemsData = [];
            app.DataListBox.Items = {};
            app.Stim_Axes.cla;
            app.KernelsListBox.Items = {};
            app.KernelsListBox.ItemsData = [];
            app.Kernel_Axes.cla;
            app.Kernel_Axes_2.cla;
            app.Kernel_Axes_3.cla;
            app.Kernel_Axes_4.cla;
            catch
                
            end
            if app.nr_r == 1 
                    app.UITable.Enable = 'off'
                  
            end
                    
        end

        % Menu selected function: MovieFrequencyMenu
        function MovieFrequencyMenuSelected(app, event)
            prompt = {'Enter frame rate for kernel movie in Hz (default = 20)'};
            dlgtitle = 'Frame rate';
            definput = {'20'};
            answer = inputdlg(prompt,dlgtitle,[1 40],definput);
            app.Kernel_info.movie_frequency = 1/(str2double(answer));
            app.Kernel_info.movie_frequency
        end

        % Button pushed function: LinkAxesButton
        function LinkAxesButtonPushed(app, event)
          if app.k_ax1.xlim ~= app.Kernel_Axes.XLim
              
              app.Kernel_Axes_2.XLim = app.Kernel_Axes.XLim;
              app.k_ax2.xlim = app.Kernel_Axes.XLim;
              
              app.Kernel_Axes_3.XLim = app.Kernel_Axes.XLim;
              app.k_ax3.xlim = app.Kernel_Axes.XLim;
              
              app.Kernel_Axes_4.XLim = app.Kernel_Axes.XLim;
              app.k_ax4.xlim = app.Kernel_Axes.XLim;
          elseif app.k_ax2.xlim ~= app.Kernel_Axes_2.XLim
              
              app.Kernel_Axes.XLim = app.Kernel_Axes_2.XLim;
              app.k_ax1.xlim = app.Kernel_Axes_2.XLim;
              
              app.Kernel_Axes_3.XLim = app.Kernel_Axes_2.XLim;
              app.k_ax3.xlim = app.Kernel_Axes_2.XLim;
              
              app.Kernel_Axes_4.XLim = app.Kernel_Axes_2.XLim;
              app.k_ax4.xlim = app.Kernel_Axes_2.XLim;
          elseif app.k_ax3.xlim ~= app.Kernel_Axes_3.XLim
              
              app.Kernel_Axes.XLim = app.Kernel_Axes_3.XLim;
              app.k_ax1.xlim = app.Kernel_Axes_3.XLim;
              
              app.Kernel_Axes_2.XLim = app.Kernel_Axes_3.XLim;
              app.k_ax2.xlim = app.Kernel_Axes_3.XLim;
              
              app.Kernel_Axes_4.XLim = app.Kernel_Axes_3.XLim;
              app.k_ax4.xlim = app.Kernel_Axes_3.XLim;
              
          elseif app.k_ax4.xlim ~= app.Kernel_Axes_4.XLim
              
               app.Kernel_Axes.XLim = app.Kernel_Axes_4.XLim;
              app.k_ax1.xlim = app.Kernel_Axes_4.XLim;
              
              app.Kernel_Axes_2.XLim = app.Kernel_Axes_4.XLim;
              app.k_ax2.xlim = app.Kernel_Axes_4.XLim;
              
              app.Kernel_Axes_3.XLim = app.Kernel_Axes_4.XLim;
              app.k_ax3.xlim = app.Kernel_Axes_4.XLim;
              
          end
          
          
          if app.k_ax1.ylim ~= app.Kernel_Axes.YLim
              
              app.Kernel_Axes_2.YLim = app.Kernel_Axes.YLim;
              app.k_ax2.ylim = app.Kernel_Axes.YLim;
              
              app.Kernel_Axes_3.YLim = app.Kernel_Axes.YLim;
              app.k_ax3.ylim = app.Kernel_Axes.YLim;
              
              app.Kernel_Axes_4.YLim = app.Kernel_Axes.YLim;
              app.k_ax4.ylim = app.Kernel_Axes.YLim;
          elseif app.k_ax2.ylim ~= app.Kernel_Axes_2.YLim
              
              app.Kernel_Axes.YLim = app.Kernel_Axes_2.YLim;
              app.k_ax1.ylim = app.Kernel_Axes_2.YLim;
              
              app.Kernel_Axes_3.YLim = app.Kernel_Axes_2.YLim;
              app.k_ax3.ylim = app.Kernel_Axes_2.YLim;
              
              app.Kernel_Axes_4.YLim = app.Kernel_Axes_2.YLim;
              app.k_ax4.ylim = app.Kernel_Axes_2.YLim;
          elseif app.k_ax3.ylim ~= app.Kernel_Axes_3.YLim
              
              app.Kernel_Axes.YLim = app.Kernel_Axes_3.YLim;
              app.k_ax1.ylim = app.Kernel_Axes_3.YLim;
              
              app.Kernel_Axes_2.YLim = app.Kernel_Axes_3.YLim;
              app.k_ax2.ylim = app.Kernel_Axes_3.YLim;
              
              app.Kernel_Axes_4.YLim = app.Kernel_Axes_3.YLim;
              app.k_ax4.ylim = app.Kernel_Axes_3.YLim;
              
          elseif app.k_ax4.ylim ~= app.Kernel_Axes_4.YLim
              
               app.Kernel_Axes.YLim = app.Kernel_Axes_4.YLim;
              app.k_ax1.ylim = app.Kernel_Axes_4.YLim;
              
              app.Kernel_Axes_2.YLim = app.Kernel_Axes_4.YLim;
              app.k_ax2.ylim = app.Kernel_Axes_4.YLim;
              
              app.Kernel_Axes_3.YLim = app.Kernel_Axes_4.YLim;
              app.k_ax3.ylim = app.Kernel_Axes_4.YLim;
              
          end
          
          
            
              
              
              
        end

        % Value changed function: CellListBox
        function CellListBoxValueChanged(app, event)
            app.StatusLamp.Color = 'red';
            value = app.CellListBox.Value;
            name_begin = strfind(value,"_")+1;
            cell_number = str2double(value(name_begin:end));
            app.add_info.FFF_select = cell_number;
            app.add_info.FFF_panel = app.PlotsPanel;
            
            
                      
            
            try
                app.selected_script = "plot_FFF_app.m";
                run_script(app);
                
            catch ME
                message = sprintf('Error loading FFF data \n%s', ME.message);
                uiwait(warndlg(message));
            end
            app.StatusLamp.Color = 'green';
            
            
            
        end

        % Value changed function: StdThresholdSlider
        function StdThresholdSliderValueChanged(app, event)
            app.Kernel_info.std_threshold = app.StdThresholdSlider.Value;
            
        end

        % Button pushed function: RejectCellQIButton
        function RejectCellQIButtonPushed(app, event)
            try
            S = load(app.data_info{app.dataselection,1},'-mat','Stimulus_info');
            Stimulus_info = S.Stimulus_info;
            Kernel_info = Stimulus_info(app.add_info.stim_idx).Kernel_info;
            idx = app.Kernel_info.Cell;
            size1 = size(Kernel_info(idx).true_kernel_idx);
            Kernel_info(idx).true_kernel_idx = zeros(size1);
            
            size2 = size(Kernel_info(idx).nr_true_kernel);
            Kernel_info(idx).nr_true_kernel = zeros(size2);
            
            size3 = size(Kernel_info(idx).true_kernel_log);
            Kernel_info(idx).true_kernel_log = zeros(size3);
            
            %Save to matfile
            Stimulus_info(app.add_info.stim_idx).Kernel_info = Kernel_info;
            app.Kernel_info = Kernel_info;
            M = matfile(app.data_info{app.dataselection,1},"Writable",true);
            M.Stimulus_info = Stimulus_info;
            
            catch ME
            message = sprintf('Error rejecting Kernel QI, check if QI is valid for this Kernel \n%s', ME.message);
            uiwait(warndlg(message));
            end
            UsequalityindexCheckBoxValueChanged(app)
            app.KernelsListBox.Value = app.Kernel_info.previous;
            
        end

        % Menu selected function: ExportsubsetofspikesMenu
        function ExportsubsetofspikesMenuSelected(app, event)
            app.StatusLamp.Color = 'r';
            app.data_info{app.dataselection,1}
            app.add_info
            try
                out = extract_dataset((app.data_info{app.dataselection,1}),app.add_info);
            catch ME
                message = sprintf('Error exporting dataset, select recording and stimulus first \n%s', ME.message);
                uiwait(warndlg(message));
            end
            app.StatusLamp.Color = 'g';
        end

        % Cell edit callback: UITable
        function UITableCellEdit(app, event)
            indices = event.Indices;
            newData = event.NewData;
            
        end

        % Menu selected function: CloseallMenu
        function CloseallMenuSelected(app, event)
            close all
        end

        % Button pushed function: CreatenewstimulusfromzoomedButton
        function CreatenewstimulusfromzoomedButtonPushed(app, event)
            %Here we get the XLim values for a zoomed in part of the stim
            %axes graph to create a new stimulus. This new created stimulus
            %will appear in the stim_list on the last position.
                      
            stim_begin = int64(app.Stim_Axes.XLim(1));
            stim_end = int64(app.Stim_Axes.XLim(2));
            
            %find next proximate trigger
            Ch = gather(app.spike_info.Ch.Ch01_02);
                   
            Ch = Ch(stim_begin:stim_end);
            
            %stim_begin = double(stim_begin);
            %stim_end = double(stim_end);
            Ch_t = Ch>3000;
            stim_end = find(Ch_t,1,'last')+stim_begin-1;
            stim_begin = find(Ch_t,1,'first')+stim_begin-1;
            
            %Add window behind last trigger
            for ii = 1:length(stim_begin)
     
                [~,locs_stim] = findpeaks(double(Ch_t),'MinPeakProminence',1,'MinPeakDistance',178);
                if isempty(locs_stim)
                    continue
                else
                stim_diff = diff(locs_stim);
                trig_max = max(stim_diff);
                stim_end = stim_end+trig_max; %Problem if the trigger channel is shorter...
                end
         
            end
            
            stim_begin = stim_begin/app.spike_info.Ch.SamplingFrequency;
            stim_end = stim_end/app.spike_info.Ch.SamplingFrequency;
            
            
            %Add to existing
            app.spike_info.stim_begin(end+1) = stim_begin;
            app.spike_info.stim_end(end+1) = stim_end;
            
            
            
            prompt = "Please enter a name for the new stimulus";
            dlgtitle = 'Name of new stimulus';
            dims = [1 35];
            definput = {['Stimulus',num2str(length(app.spike_info.stim_begin))]};
            stim_name = inputdlg(prompt,dlgtitle,dims,definput);
            stim_name{1}
            
            app.spike_info.stim_list{end+1} = stim_name{1};
            app.spike_info.stim_list
            %Update the stimulus list in app
            app.StimuliListBox.Items = app.spike_info.stim_list;
            app.StimuliListBox.ItemsData = 1:numel(app.spike_info.stim_list); %This gives an index to every Item in the list so that the stimulus can be picked later
            %Update the savefile
            M = matfile(app.data_info{app.dataselection,1},"Writable",true);
            M.stim_begin = app.spike_info.stim_begin;
            M.stim_end = app.spike_info.stim_end;
            M.stim_list = app.spike_info.stim_list;
            
            
            
            
        end

        % Value changed function: GratingsListBox
        function GratingsListBoxValueChanged(app, event)
           app.Kernel_info.previous = event.PreviousValue;
           app.GratingsListBox.Value
           Name = app.GratingsListBox.Value;
           app.Grating_info.Cell = [];
           for ii = 1:length(Name) %This allows for multiple cells to be selected at once
                app.Grating_info.Cell(ii) = str2double(regexp(Name{ii},'\d*','Match'));
           end
          
           ShowButton_2Pushed(app, event)    
           
        end

        % Button pushed function: ShowButton_2
        function ShowButton_2Pushed(app, event)
                
            
            app.StatusLamp.Color = 'red';
            drawnow() %Updates the gui
            %app.Kernel_info
        
            
            try
                app.add_info.Grating_info = app.Grating_info;
                app.add_info.Grating_info.panel_id = app.PolarPlotsPanel;
                app.add_info.Grating_info.panel_id2 = app.Panel2;
                app.selected_script = "draw_gratings.m";
                run_script(app);
                              
                
            catch ME
                message = sprintf('Error loading Kernel \n%s', ME.message);
                uiwait(warndlg(message));
            end
                
            app.StatusLamp.Color = 'green';
        end

        % Value changed function: RTestCheckBox
        function RTestCheckBoxValueChanged(app, event)
            value = app.RTestCheckBox.Value;
            if value
                app.Grating_info.rtest = 1;
                    if app.Grating_info.otest == 0
                        %Loading the datafile with the information for the gratings
                        File = load(findfile_app(app.selected_stimuli,app.spike_info.savename,'Data_circular'));
                        Data_circular = File.Data_circular;
                        Gratings_nr = length(Data_circular);
                        true_Gratings = find([Data_circular.circ_rtest_sig] == 1);
                        Gratings_name_str = {};
                        Gratings_nr = length(true_Gratings);
                        for ii = 1:Gratings_nr
                            Gratings_name_str{ii} = ['Cell_', num2str(Data_circular(true_Gratings(ii)).cell_idx)];
                        end
                    elseif app.Grating_info.otest == 1
                        %Loading the datafile with the information for the gratings
                        File = load(findfile_app(app.selected_stimuli,app.spike_info.savename,'Data_circular'));
                        Data_circular = File.Data_circular;
                        Gratings_nr = length(Data_circular);
                        true_Gratings1 = ([Data_circular.circ_rtest_sig] == 1);
                        true_Gratings2 = ([Data_circular.circ_otest_sig] == 1);
                        true_Gratings = logical(true_Gratings1.*true_Gratings2);
                        
                        true_Gratings_cell = [Data_circular(true_Gratings).cell_idx];
                        
                        Gratings_name_str = {};
                        Gratings_nr = nnz(true_Gratings);
                        for ii = 1:Gratings_nr
                            Gratings_name_str{ii} = ['Cell_', num2str(true_Gratings_cell(ii))];
                        end
                    end
                        
                
                
                
                
                
                app.GratingsListBox.Items = Gratings_name_str;
                try
                    app.GratingsListBox.Value = Gratings_name_str{1,1};
                    app.Grating_info.Cell = str2double(regexp(app.GratingsListBox.Value,'\d*','Match'));
                catch ME
                    
                end
            else
                app.Grating_info.rtest = 0;
                Tabfunction(app,"Gratings")
            end
                       
            
            
                       
        end

        % Value changed function: OTestCheckBox
        function OTestCheckBoxValueChanged(app, event)
            value = app.OTestCheckBox.Value;
            if value
                app.Grating_info.otest = 1;
            else
                app.Grating_info.otest = 0;
            end
        end

        % Menu selected function: AssignenewMenu
        function AssignenewMenuSelected(app, event)
            %Lamp update
            app.StatusLamp.Color = 'red'; %Change lamp to red
            drawnow(); %Update GUI
            %          
            
            path = uigetdir(app.appfolder,'Select main script folder');
            if isequal(path,0)
               disp('User selected Cancel')
            else
               disp(['User selected ', path])
               app.add_info.settings.location.scripts = path;
               app.scriptfolder = path;
            end
            
            %Update GUI and settings file
            savesettings(app)
            RefreshListButtonPushed(app, event)
            
            %Lamp update
            app.StatusLamp.Color = 'green'; %Change lamp to red
            drawnow(); %Update GUI
            %
        end

        % Value changed function: AllsignificantpixelCheckBox, 
        % BoxCheckBox, CILowerEditField, CIUpperEditField, 
        % ConfIntCheckBox, ConfIntCheckBox_2, GaussianCheckBox, 
        % NumRingBoxesSpinner, NumSDGaussSpinner, 
        % NumSTEbinsSpinner, NumStdDevCheckBox, 
        % NumStdDevCheckBox_2, PlotRFinformationCheckBox, 
        % SDThreshEditField, SDThreshEditField_2, STASDCheckBox, 
        % SelfCovarianceCheckBox, StimulusFrequencyEditField
        function NumStdDevCheckBoxValueChanged(app, event)
            %Frequency
            kernel_new.Hz = str2double(app.StimulusFrequencyEditField.Value);
            %RF Identification Method
            kernel_new.SS = app.STASDCheckBox.Value;
            if kernel_new.SS
                set(findall(app.QCTypePanel, '-property', 'enable'), 'enable', 'on') %Turns panel on           
            else
                set(findall(app.STASDPanel, '-property', 'enable'), 'enable', 'off')%Turns panel off     
            end
            
            kernel_new.CI = app.SelfCovarianceCheckBox.Value;
            
            if kernel_new.CI 
                set(findall(app.QCTypePanel_2, '-property', 'enable'), 'enable', 'on') %Turns panel on           
            else
                set(findall(app.SelfCovariancePanel, '-property', 'enable'), 'enable', 'off')%Turns panel off     
            end
            %STA-SD panel
            kernel_new.SS_Std = app.NumStdDevCheckBox.Value;
            if kernel_new.SS_Std
               set(app.SDThreshEditField, 'enable', 'on')
               set(app.SDThreshEditFieldLabel, 'enable', 'on')
               
               set(app.ConfIntCheckBox, 'enable', 'off')
            else
               set(app.SDThreshEditField, 'enable', 'off')
               set(app.SDThreshEditFieldLabel, 'enable', 'off')
              
            end
                            
            kernel_new.SS_CI = app.ConfIntCheckBox.Value;
            if kernel_new.SS_CI
               set(app.CIUpperEditField, 'enable', 'on')
               set(app.CIUpperEditFieldLabel, 'enable', 'on')
               
               set(app.NumStdDevCheckBox, 'enable', 'off')
            else
               set(app.CIUpperEditField, 'enable', 'off')
               set(app.CIUpperEditFieldLabel, 'enable', 'off')
                             
            end
            kernel_new.SS_SDT = str2double(app.SDThreshEditField.Value);
            kernel_new.SS_CI_Upper = str2double(app.CIUpperEditField.Value);
            
            %Self-Covar
                kernel_new.QC_Std = app.NumStdDevCheckBox_2.Value;
            if kernel_new.QC_Std
               set(app.SDThreshEditField_2, 'enable', 'on')
               set(app.SDThreshEditField_2Label, 'enable', 'on')
               
               set(app.ConfIntCheckBox_2, 'enable', 'off')
            else
               set(app.SDThreshEditField_2, 'enable', 'off')
               set(app.SDThreshEditField_2Label, 'enable', 'off')
            end
            
            kernel_new.QC_CI = app.ConfIntCheckBox_2.Value;
            if kernel_new.QC_CI
               set(app.CILowerEditField, 'enable', 'on')
               set(app.CILowerEditFieldLabel, 'enable', 'on')
               
               set(app.NumStdDevCheckBox_2, 'enable', 'off')
            else
               set(app.CILowerEditField, 'enable', 'off')
               set(app.CILowerEditFieldLabel, 'enable', 'off')
            end
            
       
             kernel_new.QC_SDT = str2double(app.SDThreshEditField_2.Value);
             kernel_new.QC_CI_Lower = str2double(app.CILowerEditField.Value);
            %RF type
            kernel_new.RF_Box = app.BoxCheckBox.Value;
            kernel_new.RF_AS = app.AllsignificantpixelCheckBox.Value;
            
            if app.AllsignificantpixelCheckBox.Value
                set(app.GaussianCheckBox,'enable', 'on')
            else
                set(app.GaussianCheckBox,'enable', 'off')
            end
            
            kernel_new.RF_Gaus = app.GaussianCheckBox.Value;
            kernel_new.RF_RBoxes = app.NumRingBoxesSpinner.Value;
            kernel_new.RF_SDGaus = app.NumSDGaussSpinner.Value;
            
            kernel_new.Num_STE_bins = app.NumSTEbinsSpinner.Value;
            kernel_new.plot_RF_overview = app.PlotRFinformationCheckBox.Value;
            %Add variables to settings 
            app.add_info.settings.kernel_new = kernel_new;
            savesettings(app);
            
            
            
            
        end

        % Button pushed function: FindRFsButton
        function FindRFsButtonPushed(app, event)
            %This button runs the RF analysis code written by Paul
            %
            app.selected_script = "RF_Ident_CL_4_GUI_Ready.m";
            run_script(app);
        end

        % Menu selected function: OnMenu
        function OnMenuSelected(app, event)
            app.add_info.settings.parpro = true;
            savesettings(app);
        end

        % Menu selected function: OFFMenu
        function OFFMenuSelected(app, event)
            app.add_info.settings.parpro = false;
            savesettings(app);
        end

        % Menu selected function: ImportdataMenu
        function ImportdataMenuSelected(app, event)
            app.StatusLamp.Color = 'red'; %Change lamp to red
            drawnow(); %Update GUI
            
            out = import_files_app; %Import dataset into gui
                                  
            
            %% inport dataset into gui
            Idx = strfind(app.data_info{:,1}, out.savename);
                            
            if Idx == 1
                error('File has already been loaded')
            end
            
            
            app.working_path(app.nr_r) = string(out.pathname); %updates the essential information in the app
            app.data_info{app.nr_r, 1} = out.savename;
            app.data_info{app.nr_r, 2} = out.all_clusters;
            app.data_info{app.nr_r, 3} = out.all_spikes;
            %                 app.savenames(app.nr_r) = data.savename;
            %                 app.spike_nr(app.nr_r) = data.all_spikes;
            %                 app.cluster_nr(app.nr_r) = data.all_clusters;
            app.UITable.Data = app.data_info;
            app.spike_info = out;
            
            if app.nr_r == 1 
                app.UITable.Enable = 'on'; %Update Table in the GUI (Table is switched off as long as
            %no recording is selected
            end
            %Update number of recordings
            app.nr_r = app.nr_r+1;

            
            
            
            
            app.StatusLamp.Color = 'green'; %Change lamp to red
            drawnow(); %Update GUI
        end

        % Button pushed function: PlotsinglecellButton_2
        function PlotsinglecellButton_2Pushed(app, event)
            cell_idx = app.CelltoplotSpinner_2.Value;
            app.add_info.kernel_raster_cell = cell_idx;
            app.selected_script = "noise_spiketrain_analysis.m";
            run_script(app)
        end

        % Value changed function: CelltoplotSpinner
        function CelltoplotSpinnerValueChanged(app, event)
            app.add_info.kernel_overview_cell = app.CelltoplotSpinner.Value;
            app.add_info.settings.kernel_new.RF_overview_panel = app.CellRFoverviewPanel;
            %Switch all plotting panels off (Before switching dynamically
            %on)
            set(app.NumStdDevCheckBox_SC_plot, 'enable', 'off')
            set(app.NumStdDevCheckBox_SS_plot, 'enable', 'off')
            set(app.ConfIntCheckBox_SC_plot, 'enable', 'off')
            set(app.ConfIntCheckBox_SS_plot, 'enable', 'off')
            set(app.BoxCheckBox_plot, 'enable', 'off')
            set(app.AllsignificantpixelCheckBox_plot, 'enable', 'off')
            set(app.GaussianCheckBox_plot, 'enable', 'off')
                        
            
            % Add panels to add_info to switch them on and off dynamically
            app.add_info.settings.kernel_new.panels.NumStdDevCheckBox_SC_plot = app.NumStdDevCheckBox_SC_plot;
            app.add_info.settings.kernel_new.panels.NumStdDevCheckBox_SS_plot = app.NumStdDevCheckBox_SS_plot;
            app.add_info.settings.kernel_new.panels.ConfIntCheckBox_SC_plot = app.ConfIntCheckBox_SC_plot;
            app.add_info.settings.kernel_new.panels.ConfIntCheckBox_SS_plot = app.ConfIntCheckBox_SS_plot;
            app.add_info.settings.kernel_new.panels.BoxCheckBox_plot = app.BoxCheckBox_plot;
            app.add_info.settings.kernel_new.panels.AllsignificantpixelCheckBox_plot = app.AllsignificantpixelCheckBox_plot;
            app.add_info.settings.kernel_new.panels.GaussianCheckBox_plot = app.GaussianCheckBox_plot;
            
            app.selected_script = "plot_RF_overview_SC.m";
            run_script(app)
        end

        % Button pushed function: PlotsinglecellButton
        function PlotsinglecellButtonPushed(app, event)
            %Collect status of the checkboxes in the plotting panel
            app.add_info.settings.kernel_new.singlecplot.SS_STA = app.NumStdDevCheckBox_SS_plot.Value;
            app.add_info.settings.kernel_new.singlecplot.SS_CI = app.ConfIntCheckBox_SS_plot.Value;
            app.add_info.settings.kernel_new.singlecplot.SC_STA = app.NumStdDevCheckBox_SC_plot.Value;
            app.add_info.settings.kernel_new.singlecplot.SC_CI = app.ConfIntCheckBox_SC_plot.Value;
            app.add_info.settings.kernel_new.singlecplot.Box = app.BoxCheckBox_plot.Value;
            app.add_info.settings.kernel_new.singlecplot.Allpix = app.AllsignificantpixelCheckBox_plot.Value;
            app.add_info.settings.kernel_new.singlecplot.Gaussian = app.GaussianCheckBox_plot.Value;
            app.add_info.settings.kernel_new.singlecplot.Heat_gray = app.GrayButton.Value;
            app.add_info.settings.kernel_new.singlecplot.Heat_colour = app.ColourButton.Value;
                       
            app.selected_script = "CL_Plotting_RF_SingleCell_v1.m";
            run_script(app)
        end
    end

    % Component initialization
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Create UIFigure and hide until all components are created
            app.UIFigure = uifigure('Visible', 'off');
            app.UIFigure.Position = [100 100 1043 748];
            app.UIFigure.Name = 'UI Figure';

            % Create FileMenu
            app.FileMenu = uimenu(app.UIFigure);
            app.FileMenu.Text = 'File';

            % Create LoadsortedspikesMenu
            app.LoadsortedspikesMenu = uimenu(app.FileMenu);
            app.LoadsortedspikesMenu.MenuSelectedFcn = createCallbackFcn(app, @LoadsortedspikesMenuSelected, true);
            app.LoadsortedspikesMenu.Text = 'Load sorted spikes';

            % Create LoaddatasetMenu
            app.LoaddatasetMenu = uimenu(app.FileMenu);
            app.LoaddatasetMenu.MenuSelectedFcn = createCallbackFcn(app, @LoaddatasetMenuSelected, true);
            app.LoaddatasetMenu.Text = 'Load dataset';

            % Create ExportdataMenu
            app.ExportdataMenu = uimenu(app.FileMenu);
            app.ExportdataMenu.Text = 'Export data';

            % Create ExportsubsetofspikesMenu
            app.ExportsubsetofspikesMenu = uimenu(app.ExportdataMenu);
            app.ExportsubsetofspikesMenu.MenuSelectedFcn = createCallbackFcn(app, @ExportsubsetofspikesMenuSelected, true);
            app.ExportsubsetofspikesMenu.Text = 'Export subset of spikes';

            % Create ExportanalysisfilesMenu
            app.ExportanalysisfilesMenu = uimenu(app.ExportdataMenu);
            app.ExportanalysisfilesMenu.Text = 'Export analysis files';

            % Create ExportfullrecordingMenu
            app.ExportfullrecordingMenu = uimenu(app.ExportdataMenu);
            app.ExportfullrecordingMenu.Text = 'Export full recording';

            % Create ImportdataMenu
            app.ImportdataMenu = uimenu(app.FileMenu);
            app.ImportdataMenu.MenuSelectedFcn = createCallbackFcn(app, @ImportdataMenuSelected, true);
            app.ImportdataMenu.Text = 'Import data';

            % Create EditMenu
            app.EditMenu = uimenu(app.UIFigure);
            app.EditMenu.Text = 'Edit';

            % Create CopypanelcontentMenu
            app.CopypanelcontentMenu = uimenu(app.EditMenu);
            app.CopypanelcontentMenu.Text = 'Copy panel content';

            % Create SettingsMenu
            app.SettingsMenu = uimenu(app.UIFigure);
            app.SettingsMenu.Text = 'Settings';

            % Create NoiseFileMenu
            app.NoiseFileMenu = uimenu(app.SettingsMenu);
            app.NoiseFileMenu.MenuSelectedFcn = createCallbackFcn(app, @NoiseFileMenuSelected, true);
            app.NoiseFileMenu.Text = 'Noise File';

            % Create KernelMenu
            app.KernelMenu = uimenu(app.SettingsMenu);
            app.KernelMenu.Text = 'Kernel';

            % Create MovieFrequencyMenu
            app.MovieFrequencyMenu = uimenu(app.KernelMenu);
            app.MovieFrequencyMenu.MenuSelectedFcn = createCallbackFcn(app, @MovieFrequencyMenuSelected, true);
            app.MovieFrequencyMenu.Text = 'Movie Frequency';

            % Create ScriptFolderMenu
            app.ScriptFolderMenu = uimenu(app.SettingsMenu);
            app.ScriptFolderMenu.Text = 'Script Folder';

            % Create AssignenewMenu
            app.AssignenewMenu = uimenu(app.ScriptFolderMenu);
            app.AssignenewMenu.MenuSelectedFcn = createCallbackFcn(app, @AssignenewMenuSelected, true);
            app.AssignenewMenu.Text = 'Assigne new';

            % Create ParallelprocessingMenu
            app.ParallelprocessingMenu = uimenu(app.SettingsMenu);
            app.ParallelprocessingMenu.Text = 'Parallel processing';

            % Create OnMenu
            app.OnMenu = uimenu(app.ParallelprocessingMenu);
            app.OnMenu.MenuSelectedFcn = createCallbackFcn(app, @OnMenuSelected, true);
            app.OnMenu.Text = 'On';

            % Create OFFMenu
            app.OFFMenu = uimenu(app.ParallelprocessingMenu);
            app.OFFMenu.MenuSelectedFcn = createCallbackFcn(app, @OFFMenuSelected, true);
            app.OFFMenu.Text = 'OFF';

            % Create WindowsMenu
            app.WindowsMenu = uimenu(app.UIFigure);
            app.WindowsMenu.Text = 'Windows';

            % Create CloseallMenu
            app.CloseallMenu = uimenu(app.WindowsMenu);
            app.CloseallMenu.MenuSelectedFcn = createCallbackFcn(app, @CloseallMenuSelected, true);
            app.CloseallMenu.Text = 'Close all';

            % Create TabGroup
            app.TabGroup = uitabgroup(app.UIFigure);
            app.TabGroup.SelectionChangedFcn = createCallbackFcn(app, @TabGroupSelectionChanged, true);
            app.TabGroup.Position = [2 40 1042 709];

            % Create MainAnalysisTab
            app.MainAnalysisTab = uitab(app.TabGroup);
            app.MainAnalysisTab.Title = 'Main Analysis';

            % Create Mainpagegrid
            app.Mainpagegrid = uigridlayout(app.MainAnalysisTab);
            app.Mainpagegrid.ColumnWidth = {'6x', '3x', '3x', '1x', '1x', '1x', 124};
            app.Mainpagegrid.RowHeight = {'2.69x', 22, '1.31x', 22, 22, 22, '14.31x', 24, '12.44x', 31, 22};
            app.Mainpagegrid.ColumnSpacing = 3.58823529411765;
            app.Mainpagegrid.RowSpacing = 2.83333333333333;
            app.Mainpagegrid.Padding = [3.58823529411765 2.83333333333333 3.58823529411765 2.83333333333333];

            % Create UITable
            app.UITable = uitable(app.Mainpagegrid);
            app.UITable.ColumnName = {'Filename'; 'Cluster Nr'; 'Spikes'; 'Column 4'};
            app.UITable.RowName = {};
            app.UITable.CellEditCallback = createCallbackFcn(app, @UITableCellEdit, true);
            app.UITable.CellSelectionCallback = createCallbackFcn(app, @UITableCellSelection, true);
            app.UITable.Layout.Row = [1 4];
            app.UITable.Layout.Column = [1 7];

            % Create RunButton
            app.RunButton = uibutton(app.Mainpagegrid, 'push');
            app.RunButton.ButtonPushedFcn = createCallbackFcn(app, @RunButtonPushed, true);
            app.RunButton.Layout.Row = 11;
            app.RunButton.Layout.Column = 6;
            app.RunButton.Text = 'Run';

            % Create RefreshListButton
            app.RefreshListButton = uibutton(app.Mainpagegrid, 'push');
            app.RefreshListButton.ButtonPushedFcn = createCallbackFcn(app, @RefreshListButtonPushed, true);
            app.RefreshListButton.Layout.Row = 11;
            app.RefreshListButton.Layout.Column = 4;
            app.RefreshListButton.Text = 'Refresh List';

            % Create ChangeFolderButton
            app.ChangeFolderButton = uibutton(app.Mainpagegrid, 'push');
            app.ChangeFolderButton.ButtonPushedFcn = createCallbackFcn(app, @ChangeFolderButtonPushed, true);
            app.ChangeFolderButton.Layout.Row = 11;
            app.ChangeFolderButton.Layout.Column = 5;
            app.ChangeFolderButton.Text = 'Change Folder';

            % Create Stim_Axes
            app.Stim_Axes = uiaxes(app.Mainpagegrid);
            title(app.Stim_Axes, 'Stim Trace')
            xlabel(app.Stim_Axes, 'X')
            ylabel(app.Stim_Axes, 'Y')
            app.Stim_Axes.PlotBoxAspectRatio = [1.73988439306358 1 1];
            app.Stim_Axes.Layout.Row = [9 10];
            app.Stim_Axes.Layout.Column = [2 3];

            % Create EditVariableButton
            app.EditVariableButton = uibutton(app.Mainpagegrid, 'push');
            app.EditVariableButton.ButtonPushedFcn = createCallbackFcn(app, @EditVariableButtonPushed, true);
            app.EditVariableButton.Layout.Row = 8;
            app.EditVariableButton.Layout.Column = 2;
            app.EditVariableButton.Text = 'Edit Variable';

            % Create RemoveselectedfileButton
            app.RemoveselectedfileButton = uibutton(app.Mainpagegrid, 'push');
            app.RemoveselectedfileButton.ButtonPushedFcn = createCallbackFcn(app, @RemoveselectedfileButtonPushed, true);
            app.RemoveselectedfileButton.Layout.Row = 5;
            app.RemoveselectedfileButton.Layout.Column = 7;
            app.RemoveselectedfileButton.Text = 'Remove selected file';

            % Create CreatenewstimulusfromzoomedButton
            app.CreatenewstimulusfromzoomedButton = uibutton(app.Mainpagegrid, 'push');
            app.CreatenewstimulusfromzoomedButton.ButtonPushedFcn = createCallbackFcn(app, @CreatenewstimulusfromzoomedButtonPushed, true);
            app.CreatenewstimulusfromzoomedButton.Layout.Row = 11;
            app.CreatenewstimulusfromzoomedButton.Layout.Column = [2 3];
            app.CreatenewstimulusfromzoomedButton.Text = 'Create new stimulus from zoomed';

            % Create DataListBox
            app.DataListBox = uilistbox(app.Mainpagegrid);
            app.DataListBox.ValueChangedFcn = createCallbackFcn(app, @DataListBoxValueChanged, true);
            app.DataListBox.Layout.Row = [6 7];
            app.DataListBox.Layout.Column = 2;

            % Create StimuliListBoxLabel
            app.StimuliListBoxLabel = uilabel(app.Mainpagegrid);
            app.StimuliListBoxLabel.HorizontalAlignment = 'center';
            app.StimuliListBoxLabel.Layout.Row = 5;
            app.StimuliListBoxLabel.Layout.Column = 3;
            app.StimuliListBoxLabel.Text = 'Stimuli';

            % Create StimuliListBox
            app.StimuliListBox = uilistbox(app.Mainpagegrid);
            app.StimuliListBox.ValueChangedFcn = createCallbackFcn(app, @StimuliListBoxValueChanged, true);
            app.StimuliListBox.Layout.Row = [6 7];
            app.StimuliListBox.Layout.Column = 3;

            % Create ScriptsListBoxLabel
            app.ScriptsListBoxLabel = uilabel(app.Mainpagegrid);
            app.ScriptsListBoxLabel.HorizontalAlignment = 'center';
            app.ScriptsListBoxLabel.Layout.Row = 5;
            app.ScriptsListBoxLabel.Layout.Column = 5;
            app.ScriptsListBoxLabel.Text = 'Scripts';

            % Create ScriptsListBox
            app.ScriptsListBox = uilistbox(app.Mainpagegrid);
            app.ScriptsListBox.ValueChangedFcn = createCallbackFcn(app, @ScriptsListBoxValueChanged, true);
            app.ScriptsListBox.Layout.Row = [6 10];
            app.ScriptsListBox.Layout.Column = [4 6];

            % Create AllspikesPanel
            app.AllspikesPanel = uipanel(app.Mainpagegrid);
            app.AllspikesPanel.TitlePosition = 'centertop';
            app.AllspikesPanel.Title = 'All spikes';
            app.AllspikesPanel.Layout.Row = [6 7];
            app.AllspikesPanel.Layout.Column = 1;

            % Create QualitycheckPanel
            app.QualitycheckPanel = uipanel(app.Mainpagegrid);
            app.QualitycheckPanel.TitlePosition = 'centertop';
            app.QualitycheckPanel.Title = 'Quality check';
            app.QualitycheckPanel.Layout.Row = 9;
            app.QualitycheckPanel.Layout.Column = 1;

            % Create DataLabel
            app.DataLabel = uilabel(app.Mainpagegrid);
            app.DataLabel.HorizontalAlignment = 'center';
            app.DataLabel.Layout.Row = 5;
            app.DataLabel.Layout.Column = 2;
            app.DataLabel.Text = 'Data';

            % Create KernelsTab
            app.KernelsTab = uitab(app.TabGroup);
            app.KernelsTab.Title = 'Kernels';

            % Create UsequalityindexCheckBox
            app.UsequalityindexCheckBox = uicheckbox(app.KernelsTab);
            app.UsequalityindexCheckBox.ValueChangedFcn = createCallbackFcn(app, @UsequalityindexCheckBoxValueChanged, true);
            app.UsequalityindexCheckBox.Text = 'Use quality index';
            app.UsequalityindexCheckBox.Position = [223 460 113 22];

            % Create ShowButton
            app.ShowButton = uibutton(app.KernelsTab, 'push');
            app.ShowButton.ButtonPushedFcn = createCallbackFcn(app, @ShowButtonPushed, true);
            app.ShowButton.Position = [10 181 89 22];
            app.ShowButton.Text = 'Show';

            % Create Kernel_Axes
            app.Kernel_Axes = uiaxes(app.KernelsTab);
            title(app.Kernel_Axes, 'UV Kernel')
            xlabel(app.Kernel_Axes, 'X')
            ylabel(app.Kernel_Axes, 'Y')
            app.Kernel_Axes.PlotBoxAspectRatio = [1.94615384615385 1 1];
            app.Kernel_Axes.Position = [385 352 300 185];

            % Create Kernel_Axes_2
            app.Kernel_Axes_2 = uiaxes(app.KernelsTab);
            title(app.Kernel_Axes_2, 'Blue Kernel')
            xlabel(app.Kernel_Axes_2, 'X')
            ylabel(app.Kernel_Axes_2, 'Y')
            app.Kernel_Axes_2.PlotBoxAspectRatio = [1.94615384615385 1 1];
            app.Kernel_Axes_2.Position = [684 352 300 185];

            % Create Kernel_Axes_4
            app.Kernel_Axes_4 = uiaxes(app.KernelsTab);
            title(app.Kernel_Axes_4, 'Red Kernel')
            xlabel(app.Kernel_Axes_4, 'X')
            ylabel(app.Kernel_Axes_4, 'Y')
            app.Kernel_Axes_4.PlotBoxAspectRatio = [1.94615384615385 1 1];
            app.Kernel_Axes_4.Position = [684 164 300 185];

            % Create Kernel_Axes_3
            app.Kernel_Axes_3 = uiaxes(app.KernelsTab);
            title(app.Kernel_Axes_3, 'Green Kernel')
            xlabel(app.Kernel_Axes_3, 'X')
            ylabel(app.Kernel_Axes_3, 'Y')
            app.Kernel_Axes_3.PlotBoxAspectRatio = [1.94615384615385 1 1];
            app.Kernel_Axes_3.Position = [384 164 300 185];

            % Create TimeSliderLabel
            app.TimeSliderLabel = uilabel(app.KernelsTab);
            app.TimeSliderLabel.HorizontalAlignment = 'right';
            app.TimeSliderLabel.Position = [715 31 32 22];
            app.TimeSliderLabel.Text = 'Time';

            % Create TimeSlider
            app.TimeSlider = uislider(app.KernelsTab);
            app.TimeSlider.ValueChangedFcn = createCallbackFcn(app, @TimeSliderValueChanged, true);
            app.TimeSlider.Position = [519 82 443 3];

            % Create PlayButton
            app.PlayButton = uibutton(app.KernelsTab, 'push');
            app.PlayButton.ButtonPushedFcn = createCallbackFcn(app, @PlayButtonPushed, true);
            app.PlayButton.Position = [408 72 100 22];
            app.PlayButton.Text = 'Play';

            % Create StopButton
            app.StopButton = uibutton(app.KernelsTab, 'push');
            app.StopButton.ButtonPushedFcn = createCallbackFcn(app, @StopButtonPushed, true);
            app.StopButton.Position = [408 50 100 22];
            app.StopButton.Text = 'Stop';

            % Create QualityCriteriaPanel
            app.QualityCriteriaPanel = uipanel(app.KernelsTab);
            app.QualityCriteriaPanel.Title = 'Quality Criteria';
            app.QualityCriteriaPanel.Position = [1 575 308 100];

            % Create ResetButton
            app.ResetButton = uibutton(app.QualityCriteriaPanel, 'push');
            app.ResetButton.Position = [221 10 86 22];
            app.ResetButton.Text = 'Reset';

            % Create CalculateButton
            app.CalculateButton = uibutton(app.QualityCriteriaPanel, 'push');
            app.CalculateButton.ButtonPushedFcn = createCallbackFcn(app, @CalculateButtonPushed, true);
            app.CalculateButton.Position = [222 43 85 22];
            app.CalculateButton.Text = 'Calculate';

            % Create StdThresholdSliderLabel
            app.StdThresholdSliderLabel = uilabel(app.QualityCriteriaPanel);
            app.StdThresholdSliderLabel.HorizontalAlignment = 'right';
            app.StdThresholdSliderLabel.Position = [69 52 80 22];
            app.StdThresholdSliderLabel.Text = 'Std Threshold';

            % Create StdThresholdSlider
            app.StdThresholdSlider = uislider(app.QualityCriteriaPanel);
            app.StdThresholdSlider.Limits = [1 4];
            app.StdThresholdSlider.MajorTicks = [1 2 3 4];
            app.StdThresholdSlider.ValueChangedFcn = createCallbackFcn(app, @StdThresholdSliderValueChanged, true);
            app.StdThresholdSlider.MinorTicks = [];
            app.StdThresholdSlider.Position = [66 40 87 3];
            app.StdThresholdSlider.Value = 2;

            % Create KernelsListBoxLabel
            app.KernelsListBoxLabel = uilabel(app.KernelsTab);
            app.KernelsListBoxLabel.HorizontalAlignment = 'right';
            app.KernelsListBoxLabel.Position = [87 554 46 22];
            app.KernelsListBoxLabel.Text = 'Kernels';

            % Create KernelsListBox
            app.KernelsListBox = uilistbox(app.KernelsTab);
            app.KernelsListBox.ValueChangedFcn = createCallbackFcn(app, @KernelsListBoxValueChanged, true);
            app.KernelsListBox.Position = [10 202 201 343];

            % Create ShowtimekernelsCheckBox
            app.ShowtimekernelsCheckBox = uicheckbox(app.KernelsTab);
            app.ShowtimekernelsCheckBox.ValueChangedFcn = createCallbackFcn(app, @ShowtimekernelsCheckBoxValueChanged, true);
            app.ShowtimekernelsCheckBox.Text = 'Show time kernels';
            app.ShowtimekernelsCheckBox.Position = [224 431 120 22];

            % Create NrneighboursSliderLabel
            app.NrneighboursSliderLabel = uilabel(app.KernelsTab);
            app.NrneighboursSliderLabel.HorizontalAlignment = 'right';
            app.NrneighboursSliderLabel.Position = [233 408 81 22];
            app.NrneighboursSliderLabel.Text = 'Nr neighbours';

            % Create NrneighboursSlider
            app.NrneighboursSlider = uislider(app.KernelsTab);
            app.NrneighboursSlider.Limits = [1 4];
            app.NrneighboursSlider.MajorTicks = [1 2 3 4];
            app.NrneighboursSlider.ValueChangedFcn = createCallbackFcn(app, @NrneighboursSliderValueChanged, true);
            app.NrneighboursSlider.MinorTicks = [];
            app.NrneighboursSlider.Position = [233 404 92 3];
            app.NrneighboursSlider.Value = 1;

            % Create LinkAxesButton
            app.LinkAxesButton = uibutton(app.KernelsTab, 'push');
            app.LinkAxesButton.ButtonPushedFcn = createCallbackFcn(app, @LinkAxesButtonPushed, true);
            app.LinkAxesButton.Position = [885 143 100 22];
            app.LinkAxesButton.Text = 'Link Axes';

            % Create RejectCellQIButton
            app.RejectCellQIButton = uibutton(app.KernelsTab, 'push');
            app.RejectCellQIButton.ButtonPushedFcn = createCallbackFcn(app, @RejectCellQIButtonPushed, true);
            app.RejectCellQIButton.Position = [225 348 100 22];
            app.RejectCellQIButton.Text = 'Reject Cell QI';

            % Create KernelsnewTab
            app.KernelsnewTab = uitab(app.TabGroup);
            app.KernelsnewTab.Title = 'Kernels new';

            % Create Kernel_new_grid
            app.Kernel_new_grid = uigridlayout(app.KernelsnewTab);
            app.Kernel_new_grid.ColumnWidth = {264, 50, 268, '1x', '1x'};
            app.Kernel_new_grid.RowHeight = {26, 66, 101, 164, 89, 97, 19, 64, 22};
            app.Kernel_new_grid.RowSpacing = 3.5;
            app.Kernel_new_grid.Padding = [10 3.5 10 3.5];

            % Create RFIdentificationMethodPanel
            app.RFIdentificationMethodPanel = uipanel(app.Kernel_new_grid);
            app.RFIdentificationMethodPanel.Title = 'RF Identification Method';
            app.RFIdentificationMethodPanel.Layout.Row = 2;
            app.RFIdentificationMethodPanel.Layout.Column = [1 2];

            % Create STASDCheckBox
            app.STASDCheckBox = uicheckbox(app.RFIdentificationMethodPanel);
            app.STASDCheckBox.ValueChangedFcn = createCallbackFcn(app, @NumStdDevCheckBoxValueChanged, true);
            app.STASDCheckBox.Text = 'STA-SD';
            app.STASDCheckBox.Position = [13 20 65 22];

            % Create SelfCovarianceCheckBox
            app.SelfCovarianceCheckBox = uicheckbox(app.RFIdentificationMethodPanel);
            app.SelfCovarianceCheckBox.ValueChangedFcn = createCallbackFcn(app, @NumStdDevCheckBoxValueChanged, true);
            app.SelfCovarianceCheckBox.Text = 'Self Covariance';
            app.SelfCovarianceCheckBox.Position = [13 -1 107 22];

            % Create RFOptionsPanel
            app.RFOptionsPanel = uipanel(app.Kernel_new_grid);
            app.RFOptionsPanel.Title = 'RF Options';
            app.RFOptionsPanel.Layout.Row = [5 7];
            app.RFOptionsPanel.Layout.Column = [1 2];

            % Create RFTypePanel
            app.RFTypePanel = uipanel(app.RFOptionsPanel);
            app.RFTypePanel.Title = 'RF Type';
            app.RFTypePanel.Position = [13 95 203 86];

            % Create BoxCheckBox
            app.BoxCheckBox = uicheckbox(app.RFTypePanel);
            app.BoxCheckBox.ValueChangedFcn = createCallbackFcn(app, @NumStdDevCheckBoxValueChanged, true);
            app.BoxCheckBox.Text = 'Box';
            app.BoxCheckBox.Position = [13 40 43 22];

            % Create AllsignificantpixelCheckBox
            app.AllsignificantpixelCheckBox = uicheckbox(app.RFTypePanel);
            app.AllsignificantpixelCheckBox.ValueChangedFcn = createCallbackFcn(app, @NumStdDevCheckBoxValueChanged, true);
            app.AllsignificantpixelCheckBox.Text = 'All significant pixel';
            app.AllsignificantpixelCheckBox.Position = [13 21 120 22];

            % Create GaussianCheckBox
            app.GaussianCheckBox = uicheckbox(app.RFTypePanel);
            app.GaussianCheckBox.ValueChangedFcn = createCallbackFcn(app, @NumStdDevCheckBoxValueChanged, true);
            app.GaussianCheckBox.Text = 'Gaussian';
            app.GaussianCheckBox.Position = [13 1 73 22];

            % Create NumSDGaussSpinnerLabel
            app.NumSDGaussSpinnerLabel = uilabel(app.RFOptionsPanel);
            app.NumSDGaussSpinnerLabel.HorizontalAlignment = 'right';
            app.NumSDGaussSpinnerLabel.Position = [16 31 89 22];
            app.NumSDGaussSpinnerLabel.Text = 'Num SD Gauss';

            % Create NumSDGaussSpinner
            app.NumSDGaussSpinner = uispinner(app.RFOptionsPanel);
            app.NumSDGaussSpinner.ValueChangedFcn = createCallbackFcn(app, @NumStdDevCheckBoxValueChanged, true);
            app.NumSDGaussSpinner.Position = [120 31 100 22];
            app.NumSDGaussSpinner.Value = 2;

            % Create NumRingBoxesSpinnerLabel
            app.NumRingBoxesSpinnerLabel = uilabel(app.RFOptionsPanel);
            app.NumRingBoxesSpinnerLabel.HorizontalAlignment = 'right';
            app.NumRingBoxesSpinnerLabel.Position = [15 63 96 22];
            app.NumRingBoxesSpinnerLabel.Text = 'Num Ring Boxes';

            % Create NumRingBoxesSpinner
            app.NumRingBoxesSpinner = uispinner(app.RFOptionsPanel);
            app.NumRingBoxesSpinner.ValueChangedFcn = createCallbackFcn(app, @NumStdDevCheckBoxValueChanged, true);
            app.NumRingBoxesSpinner.Position = [126 63 100 22];
            app.NumRingBoxesSpinner.Value = 1;

            % Create QualityControlThresholdPanel
            app.QualityControlThresholdPanel = uipanel(app.Kernel_new_grid);
            app.QualityControlThresholdPanel.Title = 'Quality Control Threshold';
            app.QualityControlThresholdPanel.Layout.Row = [3 4];
            app.QualityControlThresholdPanel.Layout.Column = [1 2];

            % Create STASDPanel
            app.STASDPanel = uipanel(app.QualityControlThresholdPanel);
            app.STASDPanel.Title = 'STA-SD';
            app.STASDPanel.Position = [12 140 251 100];

            % Create QCTypePanel
            app.QCTypePanel = uipanel(app.STASDPanel);
            app.QCTypePanel.Title = 'QC Type';
            app.QCTypePanel.Position = [13 0 120 75];

            % Create NumStdDevCheckBox
            app.NumStdDevCheckBox = uicheckbox(app.QCTypePanel);
            app.NumStdDevCheckBox.ValueChangedFcn = createCallbackFcn(app, @NumStdDevCheckBoxValueChanged, true);
            app.NumStdDevCheckBox.Text = 'Num. Std. Dev.';
            app.NumStdDevCheckBox.Position = [2 26 102 22];

            % Create ConfIntCheckBox
            app.ConfIntCheckBox = uicheckbox(app.QCTypePanel);
            app.ConfIntCheckBox.ValueChangedFcn = createCallbackFcn(app, @NumStdDevCheckBoxValueChanged, true);
            app.ConfIntCheckBox.Text = 'Conf. Int.';
            app.ConfIntCheckBox.Position = [2 5 71 22];

            % Create SDThreshEditFieldLabel
            app.SDThreshEditFieldLabel = uilabel(app.STASDPanel);
            app.SDThreshEditFieldLabel.HorizontalAlignment = 'right';
            app.SDThreshEditFieldLabel.Position = [136 47 66 22];
            app.SDThreshEditFieldLabel.Text = 'SD Thresh.';

            % Create SDThreshEditField
            app.SDThreshEditField = uieditfield(app.STASDPanel, 'text');
            app.SDThreshEditField.ValueChangedFcn = createCallbackFcn(app, @NumStdDevCheckBoxValueChanged, true);
            app.SDThreshEditField.Position = [216 47 21 22];
            app.SDThreshEditField.Value = '0';

            % Create CIUpperEditFieldLabel
            app.CIUpperEditFieldLabel = uilabel(app.STASDPanel);
            app.CIUpperEditFieldLabel.HorizontalAlignment = 'right';
            app.CIUpperEditFieldLabel.Position = [135 13 68 22];
            app.CIUpperEditFieldLabel.Text = 'CI Upper %';

            % Create CIUpperEditField
            app.CIUpperEditField = uieditfield(app.STASDPanel, 'text');
            app.CIUpperEditField.ValueChangedFcn = createCallbackFcn(app, @NumStdDevCheckBoxValueChanged, true);
            app.CIUpperEditField.Position = [216 13 21 22];
            app.CIUpperEditField.Value = '0';

            % Create SelfCovariancePanel
            app.SelfCovariancePanel = uipanel(app.QualityControlThresholdPanel);
            app.SelfCovariancePanel.Title = 'Self Covariance';
            app.SelfCovariancePanel.Position = [12 16 251 117];

            % Create QCTypePanel_2
            app.QCTypePanel_2 = uipanel(app.SelfCovariancePanel);
            app.QCTypePanel_2.Title = 'QC Type';
            app.QCTypePanel_2.Position = [13 16 120 75];

            % Create NumStdDevCheckBox_2
            app.NumStdDevCheckBox_2 = uicheckbox(app.QCTypePanel_2);
            app.NumStdDevCheckBox_2.ValueChangedFcn = createCallbackFcn(app, @NumStdDevCheckBoxValueChanged, true);
            app.NumStdDevCheckBox_2.Text = 'Num. Std. Dev';
            app.NumStdDevCheckBox_2.Position = [3 31 100 22];

            % Create ConfIntCheckBox_2
            app.ConfIntCheckBox_2 = uicheckbox(app.QCTypePanel_2);
            app.ConfIntCheckBox_2.ValueChangedFcn = createCallbackFcn(app, @NumStdDevCheckBoxValueChanged, true);
            app.ConfIntCheckBox_2.Text = 'Conf. Int.';
            app.ConfIntCheckBox_2.Position = [3 10 71 22];

            % Create SDThreshEditField_2Label
            app.SDThreshEditField_2Label = uilabel(app.SelfCovariancePanel);
            app.SDThreshEditField_2Label.HorizontalAlignment = 'right';
            app.SDThreshEditField_2Label.Position = [136 58 63 22];
            app.SDThreshEditField_2Label.Text = 'SD Thresh';

            % Create SDThreshEditField_2
            app.SDThreshEditField_2 = uieditfield(app.SelfCovariancePanel, 'text');
            app.SDThreshEditField_2.ValueChangedFcn = createCallbackFcn(app, @NumStdDevCheckBoxValueChanged, true);
            app.SDThreshEditField_2.Position = [216 58 21 22];
            app.SDThreshEditField_2.Value = '0';

            % Create CILowerEditFieldLabel
            app.CILowerEditFieldLabel = uilabel(app.SelfCovariancePanel);
            app.CILowerEditFieldLabel.HorizontalAlignment = 'right';
            app.CILowerEditFieldLabel.Position = [137 26 68 22];
            app.CILowerEditFieldLabel.Text = 'CI Lower %';

            % Create CILowerEditField
            app.CILowerEditField = uieditfield(app.SelfCovariancePanel, 'text');
            app.CILowerEditField.ValueChangedFcn = createCallbackFcn(app, @NumStdDevCheckBoxValueChanged, true);
            app.CILowerEditField.Position = [216 26 21 22];
            app.CILowerEditField.Value = '0';

            % Create StimulusFrequencyEditFieldLabel
            app.StimulusFrequencyEditFieldLabel = uilabel(app.Kernel_new_grid);
            app.StimulusFrequencyEditFieldLabel.HorizontalAlignment = 'right';
            app.StimulusFrequencyEditFieldLabel.Layout.Row = 1;
            app.StimulusFrequencyEditFieldLabel.Layout.Column = 1;
            app.StimulusFrequencyEditFieldLabel.Text = 'Stimulus Frequency';

            % Create StimulusFrequencyEditField
            app.StimulusFrequencyEditField = uieditfield(app.Kernel_new_grid, 'text');
            app.StimulusFrequencyEditField.ValueChangedFcn = createCallbackFcn(app, @NumStdDevCheckBoxValueChanged, true);
            app.StimulusFrequencyEditField.Layout.Row = 1;
            app.StimulusFrequencyEditField.Layout.Column = 2;
            app.StimulusFrequencyEditField.Value = '20';

            % Create HzLabel
            app.HzLabel = uilabel(app.Kernel_new_grid);
            app.HzLabel.Layout.Row = 1;
            app.HzLabel.Layout.Column = 3;
            app.HzLabel.Text = 'Hz';

            % Create FindRFsButton
            app.FindRFsButton = uibutton(app.Kernel_new_grid, 'push');
            app.FindRFsButton.ButtonPushedFcn = createCallbackFcn(app, @FindRFsButtonPushed, true);
            app.FindRFsButton.Layout.Row = 9;
            app.FindRFsButton.Layout.Column = [1 2];
            app.FindRFsButton.Text = 'Find RFs';

            % Create PlottingResultsPanel
            app.PlottingResultsPanel = uipanel(app.Kernel_new_grid);
            app.PlottingResultsPanel.Title = 'Plotting Results';
            app.PlottingResultsPanel.Layout.Row = [4 6];
            app.PlottingResultsPanel.Layout.Column = 3;

            % Create PlotsinglecellButton
            app.PlotsinglecellButton = uibutton(app.PlottingResultsPanel, 'push');
            app.PlotsinglecellButton.ButtonPushedFcn = createCallbackFcn(app, @PlotsinglecellButtonPushed, true);
            app.PlotsinglecellButton.Position = [11 14 100 22];
            app.PlotsinglecellButton.Text = 'Plot single cell';

            % Create CelltoplotSpinnerLabel
            app.CelltoplotSpinnerLabel = uilabel(app.PlottingResultsPanel);
            app.CelltoplotSpinnerLabel.HorizontalAlignment = 'right';
            app.CelltoplotSpinnerLabel.Position = [3 303 62 22];
            app.CelltoplotSpinnerLabel.Text = 'Cell to plot';

            % Create CelltoplotSpinner
            app.CelltoplotSpinner = uispinner(app.PlottingResultsPanel);
            app.CelltoplotSpinner.ValueChangedFcn = createCallbackFcn(app, @CelltoplotSpinnerValueChanged, true);
            app.CelltoplotSpinner.Position = [80 303 100 22];

            % Create STASDLabel
            app.STASDLabel = uilabel(app.PlottingResultsPanel);
            app.STASDLabel.Position = [31 269 49 22];
            app.STASDLabel.Text = 'STA-SD';

            % Create SelfCovarianceLabel
            app.SelfCovarianceLabel = uilabel(app.PlottingResultsPanel);
            app.SelfCovarianceLabel.Position = [149 269 90 22];
            app.SelfCovarianceLabel.Text = 'Self Covariance';

            % Create NumStdDevCheckBox_SC_plot
            app.NumStdDevCheckBox_SC_plot = uicheckbox(app.PlottingResultsPanel);
            app.NumStdDevCheckBox_SC_plot.Text = 'Num. Std. Dev.';
            app.NumStdDevCheckBox_SC_plot.Position = [143 240 102 22];

            % Create NumStdDevCheckBox_SS_plot
            app.NumStdDevCheckBox_SS_plot = uicheckbox(app.PlottingResultsPanel);
            app.NumStdDevCheckBox_SS_plot.Text = 'Num. Std. Dev.';
            app.NumStdDevCheckBox_SS_plot.Position = [5 240 102 22];

            % Create ConfIntCheckBox_SC_plot
            app.ConfIntCheckBox_SC_plot = uicheckbox(app.PlottingResultsPanel);
            app.ConfIntCheckBox_SC_plot.Text = 'Conf. Int.';
            app.ConfIntCheckBox_SC_plot.Position = [143 210 71 22];

            % Create ConfIntCheckBox_SS_plot
            app.ConfIntCheckBox_SS_plot = uicheckbox(app.PlottingResultsPanel);
            app.ConfIntCheckBox_SS_plot.Text = 'Conf. Int.';
            app.ConfIntCheckBox_SS_plot.Position = [5 209 71 22];

            % Create RFTypePanel_2
            app.RFTypePanel_2 = uipanel(app.PlottingResultsPanel);
            app.RFTypePanel_2.Title = 'RF Type';
            app.RFTypePanel_2.Position = [3 106 140 86];

            % Create BoxCheckBox_plot
            app.BoxCheckBox_plot = uicheckbox(app.RFTypePanel_2);
            app.BoxCheckBox_plot.Text = 'Box';
            app.BoxCheckBox_plot.Position = [13 40 43 22];

            % Create AllsignificantpixelCheckBox_plot
            app.AllsignificantpixelCheckBox_plot = uicheckbox(app.RFTypePanel_2);
            app.AllsignificantpixelCheckBox_plot.Text = 'All significant pixel';
            app.AllsignificantpixelCheckBox_plot.Position = [13 21 120 22];

            % Create GaussianCheckBox_plot
            app.GaussianCheckBox_plot = uicheckbox(app.RFTypePanel_2);
            app.GaussianCheckBox_plot.Text = 'Gaussian';
            app.GaussianCheckBox_plot.Position = [13 1 73 22];

            % Create HeatmapButtonGroup
            app.HeatmapButtonGroup = uibuttongroup(app.PlottingResultsPanel);
            app.HeatmapButtonGroup.Title = 'Heatmap';
            app.HeatmapButtonGroup.Position = [3 45 117 63];

            % Create GrayButton
            app.GrayButton = uiradiobutton(app.HeatmapButtonGroup);
            app.GrayButton.Text = 'Gray';
            app.GrayButton.Position = [11 17 58 22];
            app.GrayButton.Value = true;

            % Create ColourButton
            app.ColourButton = uiradiobutton(app.HeatmapButtonGroup);
            app.ColourButton.Text = 'Colour';
            app.ColourButton.Position = [61 17 65 22];

            % Create PlotstatisticsforallcellsButton
            app.PlotstatisticsforallcellsButton = uibutton(app.PlottingResultsPanel, 'push');
            app.PlotstatisticsforallcellsButton.Position = [110 13 145 23];
            app.PlotstatisticsforallcellsButton.Text = 'Plot statistics for all cells';

            % Create GeneraloptionsPanel
            app.GeneraloptionsPanel = uipanel(app.Kernel_new_grid);
            app.GeneraloptionsPanel.Title = 'General options';
            app.GeneraloptionsPanel.Layout.Row = [7 8];
            app.GeneraloptionsPanel.Layout.Column = [1 2];

            % Create NumSTEbinsSpinnerLabel
            app.NumSTEbinsSpinnerLabel = uilabel(app.GeneraloptionsPanel);
            app.NumSTEbinsSpinnerLabel.HorizontalAlignment = 'right';
            app.NumSTEbinsSpinnerLabel.Position = [1 34 83 22];
            app.NumSTEbinsSpinnerLabel.Text = 'Num STE bins';

            % Create NumSTEbinsSpinner
            app.NumSTEbinsSpinner = uispinner(app.GeneraloptionsPanel);
            app.NumSTEbinsSpinner.ValueChangedFcn = createCallbackFcn(app, @NumStdDevCheckBoxValueChanged, true);
            app.NumSTEbinsSpinner.Position = [99 34 100 22];
            app.NumSTEbinsSpinner.Value = 20;

            % Create PlotRFinformationCheckBox
            app.PlotRFinformationCheckBox = uicheckbox(app.GeneraloptionsPanel);
            app.PlotRFinformationCheckBox.ValueChangedFcn = createCallbackFcn(app, @NumStdDevCheckBoxValueChanged, true);
            app.PlotRFinformationCheckBox.Text = 'Plot RF information';
            app.PlotRFinformationCheckBox.Position = [7 4 125 22];

            % Create PlotrasterplotPanel
            app.PlotrasterplotPanel = uipanel(app.Kernel_new_grid);
            app.PlotrasterplotPanel.Title = 'Plot raster plot';
            app.PlotrasterplotPanel.Layout.Row = [7 8];
            app.PlotrasterplotPanel.Layout.Column = 3;

            % Create CelltoplotSpinner_2Label
            app.CelltoplotSpinner_2Label = uilabel(app.PlotrasterplotPanel);
            app.CelltoplotSpinner_2Label.HorizontalAlignment = 'right';
            app.CelltoplotSpinner_2Label.Position = [4 33 62 22];
            app.CelltoplotSpinner_2Label.Text = 'Cell to plot';

            % Create CelltoplotSpinner_2
            app.CelltoplotSpinner_2 = uispinner(app.PlotrasterplotPanel);
            app.CelltoplotSpinner_2.Position = [81 33 100 22];

            % Create PlotsinglecellButton_2
            app.PlotsinglecellButton_2 = uibutton(app.PlotrasterplotPanel, 'push');
            app.PlotsinglecellButton_2.ButtonPushedFcn = createCallbackFcn(app, @PlotsinglecellButton_2Pushed, true);
            app.PlotsinglecellButton_2.Position = [80 5 100 22];
            app.PlotsinglecellButton_2.Text = 'Plot single cell';

            % Create CellRFoverviewPanel
            app.CellRFoverviewPanel = uipanel(app.Kernel_new_grid);
            app.CellRFoverviewPanel.Title = 'Cell RF overview';
            app.CellRFoverviewPanel.Layout.Row = [4 5];
            app.CellRFoverviewPanel.Layout.Column = 4;

            % Create RFcellsoverviewPanel
            app.RFcellsoverviewPanel = uipanel(app.Kernel_new_grid);
            app.RFcellsoverviewPanel.Title = 'RF cells overview';
            app.RFcellsoverviewPanel.Layout.Row = [2 3];
            app.RFcellsoverviewPanel.Layout.Column = [3 4];

            % Create GratingsTab
            app.GratingsTab = uitab(app.TabGroup);
            app.GratingsTab.Title = 'Gratings';

            % Create GridLayout
            app.GridLayout = uigridlayout(app.GratingsTab);
            app.GridLayout.ColumnWidth = {'1x', '4x'};
            app.GridLayout.RowHeight = {89, 22, 24, 142, 22, '1x', '1x', '1x'};

            % Create PolarPlotsPanel
            app.PolarPlotsPanel = uipanel(app.GridLayout);
            app.PolarPlotsPanel.Title = 'Polar Plots';
            app.PolarPlotsPanel.BackgroundColor = [1 1 1];
            app.PolarPlotsPanel.Layout.Row = [1 6];
            app.PolarPlotsPanel.Layout.Column = 2;

            % Create Panel2
            app.Panel2 = uipanel(app.GridLayout);
            app.Panel2.Title = 'Panel2';
            app.Panel2.Layout.Row = [7 8];
            app.Panel2.Layout.Column = 2;

            % Create SignificancecriteriaPanel
            app.SignificancecriteriaPanel = uipanel(app.GridLayout);
            app.SignificancecriteriaPanel.Title = 'Significance criteria';
            app.SignificancecriteriaPanel.Layout.Row = 1;
            app.SignificancecriteriaPanel.Layout.Column = 1;

            % Create RTestCheckBox
            app.RTestCheckBox = uicheckbox(app.SignificancecriteriaPanel);
            app.RTestCheckBox.ValueChangedFcn = createCallbackFcn(app, @RTestCheckBoxValueChanged, true);
            app.RTestCheckBox.Text = 'RTest';
            app.RTestCheckBox.Position = [5 38 52 22];

            % Create OTestCheckBox
            app.OTestCheckBox = uicheckbox(app.SignificancecriteriaPanel);
            app.OTestCheckBox.ValueChangedFcn = createCallbackFcn(app, @OTestCheckBoxValueChanged, true);
            app.OTestCheckBox.Text = 'OTest';
            app.OTestCheckBox.Position = [5 9 53 22];

            % Create RTestTimeWindowCheckBox
            app.RTestTimeWindowCheckBox = uicheckbox(app.SignificancecriteriaPanel);
            app.RTestTimeWindowCheckBox.Text = 'RTest Time Window';
            app.RTestTimeWindowCheckBox.Position = [82 39 128 22];

            % Create OTestTimeWindowCheckBox
            app.OTestTimeWindowCheckBox = uicheckbox(app.SignificancecriteriaPanel);
            app.OTestTimeWindowCheckBox.Text = 'OTest Time Window';
            app.OTestTimeWindowCheckBox.Position = [82 9 129 22];

            % Create ShowButton_2
            app.ShowButton_2 = uibutton(app.GridLayout, 'push');
            app.ShowButton_2.ButtonPushedFcn = createCallbackFcn(app, @ShowButton_2Pushed, true);
            app.ShowButton_2.Layout.Row = 5;
            app.ShowButton_2.Layout.Column = 1;
            app.ShowButton_2.Text = 'Show';

            % Create GratingsListBox
            app.GratingsListBox = uilistbox(app.GridLayout);
            app.GratingsListBox.Multiselect = 'on';
            app.GratingsListBox.ValueChangedFcn = createCallbackFcn(app, @GratingsListBoxValueChanged, true);
            app.GratingsListBox.Layout.Row = 4;
            app.GratingsListBox.Layout.Column = 1;
            app.GratingsListBox.Value = {'Item 1'};

            % Create Grating_cells_qcPanel
            app.Grating_cells_qcPanel = uipanel(app.GridLayout);
            app.Grating_cells_qcPanel.Title = 'Grating_cells_qc';
            app.Grating_cells_qcPanel.Layout.Row = 6;
            app.Grating_cells_qcPanel.Layout.Column = 1;

            % Create Grating_qc_summaryPanel
            app.Grating_qc_summaryPanel = uipanel(app.GridLayout);
            app.Grating_qc_summaryPanel.Title = 'Grating_qc_summary';
            app.Grating_qc_summaryPanel.Layout.Row = 7;
            app.Grating_qc_summaryPanel.Layout.Column = 1;

            % Create Grating_frequencyPanel
            app.Grating_frequencyPanel = uipanel(app.GridLayout);
            app.Grating_frequencyPanel.Title = 'Grating_frequency';
            app.Grating_frequencyPanel.Layout.Row = 8;
            app.Grating_frequencyPanel.Layout.Column = 1;

            % Create FFFTab
            app.FFFTab = uitab(app.TabGroup);
            app.FFFTab.Title = 'FFF';

            % Create GridLayout2
            app.GridLayout2 = uigridlayout(app.FFFTab);
            app.GridLayout2.ColumnWidth = {'1x', '2x', '2x'};
            app.GridLayout2.RowHeight = {22, '2x', '1x'};
            app.GridLayout2.RowSpacing = 7;
            app.GridLayout2.Padding = [10 7 10 7];

            % Create PlotsPanel
            app.PlotsPanel = uipanel(app.GridLayout2);
            app.PlotsPanel.Title = 'Plots';
            app.PlotsPanel.Layout.Row = 2;
            app.PlotsPanel.Layout.Column = [2 3];

            % Create CellListBoxLabel
            app.CellListBoxLabel = uilabel(app.GridLayout2);
            app.CellListBoxLabel.HorizontalAlignment = 'center';
            app.CellListBoxLabel.Layout.Row = 1;
            app.CellListBoxLabel.Layout.Column = 1;
            app.CellListBoxLabel.Text = 'Cell';

            % Create CellListBox
            app.CellListBox = uilistbox(app.GridLayout2);
            app.CellListBox.ValueChangedFcn = createCallbackFcn(app, @CellListBoxValueChanged, true);
            app.CellListBox.Layout.Row = 2;
            app.CellListBox.Layout.Column = 1;

            % Create QualityOverviewPanel
            app.QualityOverviewPanel = uipanel(app.GridLayout2);
            app.QualityOverviewPanel.Title = 'Quality Overview';
            app.QualityOverviewPanel.Layout.Row = 3;
            app.QualityOverviewPanel.Layout.Column = 1;

            % Create AlltracesoverviewPanel
            app.AlltracesoverviewPanel = uipanel(app.GridLayout2);
            app.AlltracesoverviewPanel.Title = 'All traces overview';
            app.AlltracesoverviewPanel.Layout.Row = 3;
            app.AlltracesoverviewPanel.Layout.Column = 2;

            % Create QualityoverviewdetailedPanel
            app.QualityoverviewdetailedPanel = uipanel(app.GridLayout2);
            app.QualityoverviewdetailedPanel.Title = 'Quality overview detailed';
            app.QualityoverviewdetailedPanel.Layout.Row = 3;
            app.QualityoverviewdetailedPanel.Layout.Column = 3;

            % Create ClusteringTab
            app.ClusteringTab = uitab(app.TabGroup);
            app.ClusteringTab.Title = 'Clustering';

            % Create StatusLampLabel
            app.StatusLampLabel = uilabel(app.UIFigure);
            app.StatusLampLabel.HorizontalAlignment = 'right';
            app.StatusLampLabel.Position = [930 9 40 22];
            app.StatusLampLabel.Text = 'Status';

            % Create StatusLamp
            app.StatusLamp = uilamp(app.UIFigure);
            app.StatusLamp.Position = [985 9 20 20];

            % Create SelectedstimulusEditFieldLabel
            app.SelectedstimulusEditFieldLabel = uilabel(app.UIFigure);
            app.SelectedstimulusEditFieldLabel.HorizontalAlignment = 'right';
            app.SelectedstimulusEditFieldLabel.Position = [707 8 100 22];
            app.SelectedstimulusEditFieldLabel.Text = 'Selected stimulus';

            % Create SelectedstimulusEditField
            app.SelectedstimulusEditField = uieditfield(app.UIFigure, 'text');
            app.SelectedstimulusEditField.Position = [822 8 100 22];

            % Show the figure after all components are created
            app.UIFigure.Visible = 'on';
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = MEA_analysis_exported

            % Create UIFigure and components
            createComponents(app)

            % Register the app with App Designer
            registerApp(app, app.UIFigure)

            % Execute the startup function
            runStartupFcn(app, @startupFcn)

            if nargout == 0
                clear app
            end
        end

        % Code that executes before app deletion
        function delete(app)

            % Delete UIFigure when app is deleted
            delete(app.UIFigure)
        end
    end
end