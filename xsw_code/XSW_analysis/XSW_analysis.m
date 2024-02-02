% Chris Judd 21/09/18
function XSW_analysis
% Code for the analysis of XSW measurements.
% Currently set up for analysing DITP/TIPB on Ag(111) looking at I 3d and
% C 1s spectra for 111, -111 and 200 planes.
% Diamond experiments 17/7/18 - 26/7/18.

% Required codes:
% - i09_export_XSW.m
% - XPS_fitter.m
% - XSW_fitter.m                                                           %%%%%%%%%% Editable %%%%%%%%%%
% - Base_c_c1s.txt and Base_c_I3d.txt
%       (or Base_c.txt for different experiments)
% - fpfpp folder containing element data for XSW in specified location.
%       (See XSW_fitter).

% Codes required by other functions
% - fitpeakn.m
% - Doniach_Sunjic.m
% - fconv.m
% - Gaussian.m
% - step_fcn.m
% - q_param.m


% IMPORTANT!!!
% Code currently requires data to be stored in specifically named file
% folders to work.
% (It does this to extract data required for analysis without the user
% needing to input anything)
% Required foldername sequence: (Ignore "[]")
% .../Diamond data\[Molecule name]\[]\[Orbital]\[Reflecting plane]
% E.g. : \Diamond Data\DITP\Pre-Anneal\I3d\111

% Editing for different experiments
% Important sections that require editing for each different set of
% experimental conditions are marked:
% " %%%%%%%%%% Editable %%%%%%%%%% "
% at the end of the line.
% These may be values that need changing (e.g. Bragg energy values,
% reflecting planes etc) or a current part of the code relying of
% foldernames (see above).
% There may also be editable sections in "required codes" fucntions (listed
% above) not in this .m file.

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Load Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Loads in a series of XSW measurements for the same bragg reflection. The
% average of all of these measurements is then used for fitting data (shown
% in graphs on the right).
% Initial fitting parameters and number of peaks can be set here as well as
% loading previosly saved fitting parameters. All energies are measured
% relative to the energy of the first peak.

% %%%%%%%%%%%%%%%%%%%%%%%%%%%% Fitting Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fitting program uses a convolution of a guassian and a doniach sunjic
% with a step function, to fit each peak.
% Up to 6 peaks can be fitted to any data set.

% The fitting program will display the averaged data (red circles) and the
% fitted curves for each photon energy on the rocking curve. Once finished,
% the fit at each photon energy can be viewed by using the left and right
% arrows.

% Fitting is optimised by varying each peak position and widths, adjusting
% the range of values they can take and setting whether they can vary or
% not (tick boxes). Clicking a row on the parameter list displays that
% parameters variation with photon energy. Clicking teh average box sets
% the inital value of that parameter to the current average. Ticking the
% 'Varies' tick box off sets that parameter to be constant for all photon
% energies at the initial value.

% Fitting parameters can be saved using the save button and will be saved
% as a text file in the current directory as 'c_data.txt'. This will be
% loaded as the default file next time this dataset is selected.

% A good method for optimising the fit is to go find good initial positions
% for the peaks and allow everything to vary within a reasonable range.
% Then find a parameter that does not vary much with photon energy and fix
% that at its average value, then refit. Repeat until all variables are
% sensible and then go back and allow individual parameters to vary to
% optimise their positions.

% All parameters except brg0, brg1 and intensity1,2,3... should be
% approximately constant with photon energy and fitting should be optimised
% to allow this. Once this is done press XSW to begin XSW fitting.

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% XSW %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fit button fits an XSW curve to the intensity data (varying with photon
% energy) for the selected peak and saves it as a .mat, .png and .svg file
% if 'Save data ?' is selected.

% Main output graph shows the fitted curve in red and data with errorbars
% in green. Coherent fraction and position are displayed.

% Smaller colormap figure shows grid of all possible coherent fractions and
% positions with red being large error and blue meaning small error.

% Drop down box allows user to select which peak to fit. Sum calculates fit
% from sum of several peaks (currently up to 1,3,4 and 6).


clear; clc; close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Load data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Loads in all user selected XSW data files and extracts data from them.
[filenames,pathname] = uigetfile({'*.nxs','nexus'},'MultiSelect','on');
%[data_square,i0,hv,nixswr,angles,energy,polar] = i09_export_XSW(filenames,pathname);
[data_square,i0,hv,nixswr,angles,energy,polar] = i09_export_XPS1(filenames,pathname);

% Change FE= to whatever you get from the fermi energy file

FE=5.04;
energy= energy-FE;

% Averages the data if more than one dataset is loaded.
data_square4d = data_square;
if ndims(data_square) == 4
    hv = hv(1:size(data_square,2));
    nixswr = nixswr(1:size(data_square,2));
    data_square = mean(data_square4d,4);
end

% Plots all loaded and averaged data, summing over photon energy and
% binding energy (top and bottom graphs respectively. Also prpeares data
% for exporting to fitting function.
figure(2)
set(gcf,'units','normalized','outerposition',[0 0.04 1 0.96])
if ndims(data_square) > 2
    subplot(2,3,[2,3]); hold all; xlabel('Binding Energy /eV'); ylabel('Intensity'); title('Sum over hv')
    subplot(2,3,[5,6]); hold all; xlabel('hv /eV'); ylabel('Intensity'); title('Sum over binding energies');
    for n = 1:size(data_square4d,4) + 1
        % For all loaded data.
        if n <= size(data_square4d,4)
            for nn = 1:length(hv)
                % Prepares data for fitting function.
                data_export(nn,:) = squeeze(sum(data_square4d(:,nn,:,n),1))';
            end
            subplot(2,3,[5,6]); plot(hv,mean(data_export,2),'.')
            subplot(2,3,[2,3]); plot(energy,mean(data_export,1),'.')
            % For Averaged data.
        else
            for nn = 1:length(hv)
                % Prepares data for fitting function.
                data_export(nn,:) = squeeze(sum(data_square(:,nn,:),1))';
                % For XPS data
                %data_export = data_square;
            end
            subplot(2,3,[5,6]); plot(hv,mean(data_export,2))
            subplot(2,3,[2,3]); plot(energy,mean(data_export,1))
        end
    end
else
    
    data_export = mean(data_square,1);
    subplot(2,3,[2,3,5,6]); 
    hold all
    plot(energy,data_square,'.')
    plot(energy,data_export,'-')
    
    
    
end
% Loads in initial fitting parameters.
% c_data is loaded if it already exists as previous fitting parameters.
% Otherwise relevant bace c_data is loaded.
if exist([pathname,'c_data.txt'],'file') == 2                              %%%%%%%%%% Editable %%%%%%%%%%
    c_data = dlmread([pathname,'c_data.txt']);
elseif contains(pathname,'I3d')
    c_data = dlmread('Base_c_I3d.txt');
elseif contains(pathname,'C1s')
    c_data = dlmread('Base_c_C1s.txt');
else
    c_data = dlmread('Base_c.txt');
end

% Prepares List of parameters.
titles = {  'bgr0','bgr1',  'asym1',  'step_h',  'lorentz_width', 'intensity1',   'energy1',  'width1',  'intensity2',  'energy2',  'width2',  'intensity3',  'energy3',   'width3','intensity4',   'energy4',  'width4',  'intensity5',  'energy5',  'width5',  'intensity6',  'energy6',   'width6','asym2','bgr3','bgr4'};
set(gcf,'units','normalized')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% UIs and callback functions %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Parameter list table.
c_table = uitable('Data',[num2cell(c_data(:,1:3))';num2cell(logical(c_data(:,4)))']',...
    'rowname',titles',...
    'columnname',{'Value','Upper Bound','Lower Bound','Variable?'},...
    'units','normalized',...
    'position',[0.05,0.05,0.25,0.67],...
    'fontsize',12,...
    'columnEditable',logical([1,1,1,1]),...
    'columnformat',{[],[],[],'logical'});

% Load previous parameter (c_data) button.
load_c_data_button = uicontrol('Style','pushbutton',...
    'units','normalized',...
    'position',[0.2,0.75,0.1,0.08],...
    'string',sprintf('Load preset parameters'),...
    'fontsize',12,...
    'callback',@load_c_data);

% Runs initial fitting program and advances to next stage.
fit_button = uicontrol('Style','pushbutton',...
    'units','normalized',...
    'position',[0.05,0.75,0.1,0.08],...
    'string','Fit',...
    'fontsize',30,...
    'callback',@run_fit_first);

% Dropdown box for selecting number of peaks.
peak_no_select = uicontrol('Style','popupmenu',...
    'units','normalized',...
    'position',[0.2,0.82,0.1,0.05],...
    'string',[1;2;3;4;5;6],...
    'fontsize',16,...
    'Callback',@peak_no_listen);

% Dsplays currently slected number of peaks.
peak_no_text = uicontrol('style','text',...
    'units','normalized',...
    'position',[0.2,0.87,0.1,0.05],...
    'string', sprintf('Peaks = %i',peak_no_select.Value),...
    'fontsize',16);

    function load_c_data(source,event)
        % Callback function for loading previous c_data parameters.
        [loadfile, loadpath] = uigetfile('.txt');
        c_data = dlmread([loadpath, loadfile]);
        set(c_table,'Data', [num2cell(c_data(:,1:3))';num2cell(logical(c_data(:,4)))']')
    end

    function peak_no_listen(source,event)
        % Callback function for changing the peak number text.
        set(peak_no_text,'string', sprintf('Peaks = %i',peak_no_select.Value))
    end

    function run_fit_first(source,event)
        % Callback for running the fitting program. Advances to next step
        % in fitting.
        c_data = [cell2mat(c_table.Data(:,1:3)),cell2mat(c_table.Data(:,4))];
        XPS_fit_main(data_export,hv,energy,c_data,peak_no_select.Value,nixswr,polar);
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Fitting Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function XPS_fit_main(data_square,hv,energy,c_data,no_peaks,nixswr,polar)
% Function for optimising the fitting parameters for XSW data prior to
% final fitting of XSW curve.
figure(3)
% Runs fitting code.
% Applies a least square non-linear minimisation algorithm to fit input
% data (data_square) with fitting parameters c_data and number of peaks
% (no_peaks).
% An extra figure is produced showing the fitted data for eahc photon
% energy. Top plot shows the residuals from the fit, bottom shows the fit
% itself. After fitting is applied to all photon energies, figures and be
% scrolled through using the arrow keys.

[c_out,titles,err] = XPS_fitter(data_square,hv,energy,c_data,no_peaks);

% Adding peak integrals, averaged over all hv
peak_int = peak_integrate(hv,energy,c_data,c_out,no_peaks);
for n = 1:no_peaks
    leg_string{n} = ['peak ',num2str(n), ': ', num2str(mean(peak_int(:,n),1),'%0.1f')];
end
legend(leg_string)
figure(2)
clf
pathname = cd;

% Plot first c_data parameter (bgr0) from fitted function varying with
% photon energy.
subplot(1,3,2:3)
plot_c_out = plot(hv,c_out(:,1));
xlabel('hv /eV')
ylabel(titles(1))
title(titles(1))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% UIs and callback functions %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Text showing the file name of the data currently processing.
file_text = uicontrol('style','text',...
    'units','normalized',...
    'position',[0.05,0.95,0.3,0.03],...
    'string', strrep(pathname(strfind(pathname,'Data')+5:end),'\',' '),... %%%%%%%%%% Editable %%%%%%%%%%
    'fontsize',16);

% Interactive table displaying parameter name, initial value, maximum and
% minimum values and whether it is set to vary or be fixed. Selecting any
% of the values (except name) will display a graph showinghow the parameter
% varies with photon energy.
c_table = uitable('Data',[num2cell(c_data(:,1:3))';num2cell(logical(c_data(:,4)))']',...
    'rowname',titles',...
    'columnname',{'Value','Upper Bound','Lower Bound','Variable?'},...
    'units','normalized',...
    'position',[0.05,0.05,0.24,0.67],...
    'fontsize',12,...
    'columnEditable',logical([1,1,1,1]),...
    'columnformat',{[],[],[],'logical'},...
    'cellselectioncallback',@plot_new,...
    'celleditcallback',@table_update);

% Table showing average values for each parameter. Selecting one will set
% the initial value to the average.
avg_table = uitable('Data',num2cell(mean(c_out,1))',...
    'columnname',{'Average'},...
    'rowname',{},...
    'units','normalized',...
    'position',[0.3,0.05,0.06,0.67],...
    'fontsize',12,...
    'columnEditable',false,...
    'columnformat',{[]},...
    'cellselectioncallback',@avg_val);

% Reruns the fitting program.
fit_button = uicontrol('Style','pushbutton',...
    'units','normalized',...
    'position',[0.05,0.75,0.1,0.08],...
    'string','Refit',...
    'fontsize',30,...
    'callback',@run_fit);

% Saves all current parameter (c_data) values in c_data.txt in current
% directory.
save_params_button = uicontrol('Style','pushbutton',...
    'units','normalized',...
    'position',[0.05,0.85,0.1,0.08],...
    'string','Save Parameters',...
    'fontsize',18,...
    'callback',@save_params);

% Selects number of peaks to fit.
peak_no_select = uicontrol('Style','popupmenu',...
    'units','normalized',...
    'position',[0.2,0.74,0.1,0.05],...
    'string',[1;2;3;4;5;6],...
    'value',no_peaks,...
    'fontsize',16,...
    'Callback',@peak_no_listen);

% Displays current number of peaks fitting.
peak_no_text = uicontrol('style','text',...
    'units','normalized',...
    'position',[0.2,0.79,0.1,0.03],...
    'string', sprintf('Peaks = %i',peak_no_select.Value),...
    'fontsize',16);

% Advances code to fitting XSW data based on current fit.
% Should only do this if all values are fixed except bgr0, bgr1 and
% relevant intensities.
XSW_button = uicontrol('Style','pushbutton',...
    'units','normalized',...
    'position',[0.2,0.85,0.1,0.08],...
    'string','XSW',...
    'fontsize',30,...
    'callback',@run_XSW);

    function plot_new(source,event)
        % Callback function to plot selected parameter varying with photon
        % energy.
        if isempty(event.Indices) == 0
            a = event.Indices(1);
            subplot(1,3,2:3)
            
            set(plot_c_out,'XData',hv,'Ydata',c_out(:,a))
            xlabel('hv /eV')
            ylabel(titles(a))
            title(titles(a))
        end
    end


    function run_fit(source,event)
        % Callback funtion for 'refit' button.
        % Extracts c_data from parameter table and reruns parent function.
        c_data = [cell2mat(c_table.Data(:,1:3)),cell2mat(c_table.Data(:,4))];
        XPS_fit_main(data_square,hv,energy,c_data,no_peaks,nixswr,polar)
    end

    function peak_no_listen(source,event)
        % Callback function for peak number selection.
        % Adjusts peak number text.
        no_peaks = peak_no_select.Value;
        set(peak_no_text,'string', sprintf('Peaks = %i',peak_no_select.Value))
    end

    function save_params(source,event)
        % Callback function for save params button.
        % Saves parameter table data as c_data.txt in current directory.
        c_data = [cell2mat(c_table.Data(:,1:3)),cell2mat(c_table.Data(:,4))];
        dlmwrite([cd,'\c_data.txt'],c_data)
    end

    function avg_val(source,event)
        % Callback function for average table.
        % Changes selected parameter initial value to the selected average.
        a = event.Indices(1);
        Previous_data = c_data(a,1);
        c_data(a,1) = source.Data{a,1};
        c_table.Data = [num2cell(c_data(:,1:3))';num2cell(logical(c_data(:,4)))']';
        
        events.Indices = event.Indices;
        events.PreviousData = Previous_data;
        table_update(c_table,events)
    end

    function table_update(source,event)
        % Callback function for when data is changed in the parameter
        % table.
        % Adjusts c_data variable and checks that all values obey maximum
        % and minimum constraints.
        c_data = [cell2mat(c_table.Data(:,1:3)),cell2mat(c_table.Data(:,4))];
        a = event.Indices(1);
        b = event.Indices(2);
        if c_data(a,1) > c_data(a,2) || c_data(a,1) < c_data(a,3)
            c_data(a,b) = event.PreviousData;
            source.Data = [num2cell(c_data(:,1:3))';num2cell(logical(c_data(:,4)))']';
            warning('Value must be between upper and lower bounds')
        elseif c_data(a,2) < c_data(a,3)
            c_data(a,b) = event.PreviousData;
            source.Data = [num2cell(c_data(:,1:3))';num2cell(logical(c_data(:,4)))']';
            warning('Upper bound must be greater than lower bound')
        elseif c_data(a,3) > c_data(a,2)
            c_data(a,b) = event.PreviousData;
            source.Data = [num2cell(c_data(:,1:3))';num2cell(logical(c_data(:,4)))']';
            warning('Lower bound must be less than upper bound')
        end
    end

    function run_XSW(source,event)
        % Callback for XSW button.
        % Advances code to next function used for XSW fitting.
        XSW_fit_main(hv,c_out,nixswr,err,polar)
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% XSW %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function XSW_fit_main(hv,c_out,nixswr,err,polar)
% Fits XSW curve based on fitting parameters found previously.
figure(1)
set(gcf,'units','normalized','outerposition',[0 0.04 1 0.96])
clf
no_peaks = size(err,1);
peak_no = 1;

% Button to start XSW fit for currently selected peak.
XSW_button = uicontrol('Style','pushbutton',...
    'units','normalized',...
    'position',[0.05,0.8,0.1,0.08],...
    'string','Fit',...
    'fontsize',30,...
    'callback',@run_XSW);

% Text displaying selected peak.
peak_no_text = uicontrol('style','text',...
    'units','normalized',...
    'position',[0.05,0.75,0.1,0.03],...
    'string', 'Peak to fit',...
    'fontsize',16);

% Drop down box for selecting which peak to fit.
peak_no_select = uicontrol('Style','popupmenu',...
    'units','normalized',...
    'position',[0.05,0.7,0.1,0.05],...
    'string',{1:no_peaks,'Sum'},...
    'value',1,...
    'fontsize',16,...
    'Callback',@peak_no_listen);

% Checkbos for determining whether images of final data are saved.
save_check = uicontrol('Style','checkbox',...
    'units','normalized',...
    'position',[0.05,0.65,0.1,0.05],...
    'string','   Save data ?',...
    'fontsize',16,...
    'Value',1);

% Sets the filename for outputs to be based on current path and peak
% number.
pathname = cd;
plane_str = pathname(find(pathname=='\',1,'last')+1:end);                  %%%%%%%%%% Editable %%%%%%%%%%
filename = ['XSW_',strrep(pathname(strfind(pathname,'Data')+5:end),'\','_'),'_',num2str(peak_no)];

% Sets plane and bragg energies based on file location.
if strcmp(plane_str,'111')                                                 %%%%%%%%%% Editable %%%%%%%%%%
    plane = [1,1,1];
    bragg = 2631.2;
elseif strcmp(plane_str,'-111')
    plane = [-1,1,1];
    bragg = 2631.2;     %%%maybe this is wrong? ask chris%%%%%
elseif strcmp(plane_str,'200')
    plane = [2,0,0];
    bragg = 3035;
end

    function peak_no_listen(source,event)
        % Callback function for changing peak text.
        peak_no = peak_no_select.Value;
        filename(end) = num2str(peak_no);
    end

    function run_XSW(source,event)
        % Callback function for running XSW fit.
        % Deletes prvious figures.
        delete(findobj(gcf,'type','axes'))
        set(source,'string','Fitting...')
        pause(0.1)
        
        % If 'all' is not selected in dropdown.
        if peak_no <= no_peaks
            % Extracts error for particular peak
            err_in = err(peak_no,:);
            % Runs XSW fitting program.
            fpfpp_dir = 'C:\Users\pcxef1\OneDrive - The University of Nottingham\Diamond Data\NIXSW analysis guide\XSW codes\XSW_analysis\fpfpp'; %%%%%%%%%% Editable %%%%%%%%%%
            [q,err_fh,err_ph] = XSW_fitter(fpfpp_dir,filename,[hv',c_out(:,3*(peak_no+1)),nixswr',err_in'],15,1,bragg,plane,-polar,save_check.Value); %%%%%%%%%% Editable %%%%%%%%%%
            
        else                                                               %%%%%%%%%% Editable %%%%%%%%%%
            % If 'Sum' is selected in dropdown.
            % Runs XSW fitting based on sum of several peaks.
            % Sums all peaks apart from excluded peaks (peak_exc)
            
            % Adjusts filename
            filename = ['XSW_',strrep(pathname(strfind(pathname,'Data')+5:end),'\','_'),'_a'];
            % Excludes these peaks from the summation
            peak_exc = [0];                                                %%%%%%%%%% Editable %%%%%%%%%%
            % Calculates sum of all peaks and their errors for all peaks
            % not excluded.
            c_out_all = 0;
            err_in_all = 0;
            for na = 1:no_peaks
                if na ~= peak_exc
                    c_out_all = c_out_all+c_out(:,3*(na+1));
                    err_in_all = err_in_all + err(na,:).^2;
                end
            end
            err_in_all = sqrt(err_in_all);
            
            % Runs XSW fitting program.
            fpfpp_dir = 'C:\Users\pcxef1\OneDrive - The University of Nottingham\Diamond Data\NIXSW analysis guide\XSW codes\XSW_analysis\fpfpp'; %%%%%%%%%% Editable %%%%%%%%%%
            [q,err_fh,err_ph] = XSW_fitter(fpfpp_dir,filename,[hv',c_out_all,nixswr',err_in_all'],15,1,bragg,plane,-polar,save_check.Value); %%%%%%%%%% Editable %%%%%%%%%%
        end
        set(source,'string','Fit')
    end
end