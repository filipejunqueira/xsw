%% Loads in data and plots
clear; 
% Select XSW data files to read
[filenames,pathname] = uigetfile({'*.nxs','nexus'},'MultiSelect','on');
% Reads XSW data files (data square) and important parameters
[data_square,i0,hv,nixswr,angles,energy,polar] = i09_export_XSW(filenames,pathname);
% data_square format: [detector energies (eV),photon energy (hv), binding energy (energy), Repeats (optional)]
% i0 = [?]
% hv = X-ray photon energy
% nixswr = rocking curve intensity measurement
% angles = [?]
% energy = measured binding energy of detected photo-electrons
% polar = polar angle of detector


% Calacutes average of data if more than one data set loaded
data_square4d = data_square;
if ndims(data_square) == 4
    hv = hv(1:size(data_square,2));
    nixswr = nixswr(1:size(data_square,2));
    data_square = mean(data_square4d,4);    
end

% Plots all data and compares to average
figure(2); clf; hold all; xlabel('Binding Energy /eV'); ylabel('Intensity'); title('Sum over hv')
figure(1); clf; hold all; xlabel('hv /eV'); ylabel('Intensity'); title('Sum over binding energies'); 
for n = 1:size(data_square4d,4) + 1
    % Multiple datasets loaded
    if n <= size(data_square4d,4)
        for nn = 1:length(hv)
            data_export(nn,:) = squeeze(sum(data_square4d(:,nn,:,n),1))';
        end
        figure(1); plot(hv,mean(data_export,2),'.')
        figure(2); plot(energy,mean(data_export,1),'.')
    % Single dataset loaded    
    else
        for nn = 1:length(hv)
            data_export(nn,:) = squeeze(sum(data_square(:,nn,:)))';
        end
        figure(1); plot(hv,mean(data_export,2))
        figure(2); plot(energy,mean(data_export,1))
    end
end

% Plot 3D data
% figure(3)
% clf
% [mesh_hv,mesh_energy] = meshgrid(hv,energy);
% surface(mesh_hv,mesh_energy,data_export','EdgeColor','none')
% xlabel('hv /eV')
% ylabel('binding energy /eV')
% zlabel('Intensity')

% Rocking Curve plot
figure(4)
plot(hv,nixswr)

%% XPS Fitting - Using pre-prepared c_data file
c_data = dlmread([cd,'\c_data.txt']);
no_peaks = 4;
[c_out,titles,err] = XPS_fitter(data_export,hv,energy,c_data,no_peaks);
%% XPS fitting using base xps file - generally don't need to use this
[c_out,titles,err] = XPS_Base(data_export,hv,energy);

%dlmwrite([pathname,strrep(pathname(strfind(pathname,'Data'):end-1),'\','_'),'.txt'],c_out);


% need to check normaliasing each data point to i0

%     bgr0         =    [600,      1e6,   0.0005       1];   
%     bgr1         =    [-4e-5, 2e5,    -2e5            1]; 
%     
%     asym1 =           [0.0489    0.4     1e-5    1];
%    
%     lorentz_width =   [0.9778   1.5    0.1    1];
%     
%     step_h =          [0.001    0.9     0.00   1];    
%    
%     energy1         = [620,  621, 619,    1];
%     intensity1      = [2000     1e7     0       1];
%     width1          = [0.5951    0.9     0.2    1];
%     
%     energy2         = [2.5,   3.5,  1.5   1];
%     intensity2      = [250    1e6     1e-8    1];
%     width2          = [5.4448    8  0.2    1];

%% Display c_out final values - generally don't need to use this
figure(4)
for a =1:size(titles,2)
    disp(['n = ',num2str(a),', ',titles{a},' Average = ',num2str(mean(c_out(:,a)))])
    plot(hv,c_out(:,a))
    xlabel('hv /eV')
    ylabel(titles(a))
    title(titles(a))
    pause
end
%% Run XSW fitting programXSW
fpfpp_location = 'C:\Users\ppxcj\Documents\Experiments\XSW\XSW_codes\David\fpfpp';

[q,err_fh,err_ph] = XSW_Ag_111_I3d(fpfpp_location,'I3d_test1',[hv',c_out(:,6),nixswr',err'],15,1);
%[q,err_fh,err_ph] = XSW_Ag_111_C1s(fpfpp_location,'C1s_test1',[hv',c_out(:,9),nixswr',err_1',err_2'],15,1);

%% Manual XSW fitting 
% IMPORTANT NOTE: This takes a long time to start, BE PATIENT it is running

addpath('C:\Users\ppxcj\Documents\Experiments\XSW\XSW_codes\My_code\Bodges - For individual_XSW_codes')

% Ensure bragg energy and reflecting plane is correct otherwise you get
% errors
bragg = 2629;%3035
plane = [1,1,1];

f = 0.5;
p = 0.5;
peak_no = 2;
name = [strrep(pathname(strfind(pathname,'Diamond Data')+13:end),'\',' '),num2str(peak_no)];
[f,p] = XSW_bodger(hv,c_out(:,3*(peak_no+1)),nixswr,err(peak_no,:),bragg,-polar,f,p,plane,name);
%% Cycles through all possible XSW curves (General not necessary)
% IMPORTANT NOTE: This takes a long time to start, BE PATIENT it is running

subplot('Position',[0.32 0.1 0.63 0.6]);
cla
clc
fa = linspace(0.0,1,10);
pa = linspace(0.01,1.0,10);
q_scale = 0;
for n1 = 1:numel(fa)
    for n2 = 1:numel(pa)
    q = [1,fa(n1),pa(n2)];

    pause(0.01)
    q_scale = XSW_Bodge('C:\Users\ppxcj\Documents\Experiments\XSW\XSW_codes\My_code\fpfpp',[hv',c_out(:,3*(peak_no+1)),nixswr',err(peak_no,:)'],15,q,bragg,-polar,plane,q_scale);
    fclose all;
    end
end
disp('done')
%% Integrate peak areas (for use after peak fitting is finalised
clf
for n = 1:size(hv,2)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Taken from XPS_fitter function
    % All code between '%' lines used to fit peaks to XPS data
    cs = c_data(:,1)';
    flags = logical(c_data(:,4))';
    
    
    energy_hold = energy;
    energy_hold(find(energy_hold(1:end-1)-energy(2:end)==0)) = [];
    e_range = 2*(max(energy)-min(energy));
    estep = abs(min(energy_hold(1:end-1)-energy_hold(2:end)));
    new_e = min(energy)-e_range:estep:max(energy)+e_range;
    
    for bb = 1:length(energy)
        abs_newe = abs(new_e-energy(bb));
        index_e(bb) = find(abs_newe==min(abs_newe),1);
    end
    
    
    %[lineshape,di] = fit_peakn(c_out(n,flags),cs,flags,energy,1,6,new_e,index_e,no_peaks,3);
    
    c = c_out(n,flags);
    ein = energy;
    background = 1;
    di = 6;
    plotter = 0;
    
    colors = 'gbmcrk';
    cs(flags) = c;
    zz = 1;
    bgr0 =          cs(1);
    bgr1 =          cs(2);
    asym1 =         cs(3);
    step_h =        cs(4);
    lorentz_width = cs(5);
    for kk = 1:no_peaks
        intensity1 =    cs(di);di = di+1;
        if kk == 1
            energy1 =       cs(di);di = di+1;
            e1 = energy1(1);
            bck = (bgr0(1)+(ein(:)-e1(1))*bgr1(1)).*background(:);
            lineshape = bck(:);
        else
            energy2 =       cs(di);di = di+1;
            e1 = energy1(1)+energy2(1);
        end
        width1 =        cs(di);di = di+1;
        peak1_G = Gaussian(new_e,[intensity1(1),e1,width1(1)]);
        peak1_DS = Doniach_Sunjic(-new_e,[intensity1(1),asym1(1),-e1,lorentz_width(1)]);
        peak1 = fconv(peak1_DS,peak1_G);
        peak1 = peak1./max(peak1)*intensity1(1);
        peak1 = peak1(1:2:end);
        peak1 = peak1(index_e);
        step1 = step_fcn(ein,[intensity1(1),e1,width1(1)]);
        step = step_h(1)*step1;
        lineshape = lineshape + peak1(:)+step;
        if plotter > 0
            %figure(plotter)
            obj = (findobj(gca,'type','line'));
            if length(obj) < no_peaks+2
                plot(ein,peak1(:)+bck(:),colors(zz),'Linewidth',2)
                set(gca,'Xdir','reverse')
            else
                set(obj(no_peaks-kk+3),'Xdata',ein,'Ydata',peak1(:)+bck(:));
            end
            hold on
            zz = zz+1;
            if zz>length(colors)
                zz = 1;
            end
        end
        peak_all(kk,:) = peak1;
    end
    
    clf
    plot(energy,peak_all)
    pause
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Integrate area under fitted peaks
    
    for nn = 1:no_peaks
        int(n,nn) = -trapz(energy,peak_all(nn,:));
    end
end



% Plots variation of peak intensity with photon energy (hv)
plot(hv,int')
legend(['peak 1';'peak 2';'peak 3';'peak 4'])
mean(int)
max(int)
disp(mean(int)./min(mean(int)))