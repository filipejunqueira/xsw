function [lineshape,di] = fit_peakn(c,cs,flags,ein,background,di,new_e,index_e,no_peaks,plotter)
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
    end
end