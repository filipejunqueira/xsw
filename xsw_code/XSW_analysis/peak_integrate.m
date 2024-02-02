function peak_int = peak_integrate(hv,energy,c_data,c_out,no_peaks)


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
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Integrate area under fitted peaks
    
    for nn = 1:no_peaks
        peak_int(n,nn) = -trapz(energy,peak_all(nn,:));
    end
end



