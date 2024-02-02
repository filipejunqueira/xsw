function [c_out,titles,err] = XPS_forced(data_square,hv,energy,c_data,no_peaks)
    
    PhD_data_square = data_square;

    
    %                   1       2        3         4              5    %         6              7           8           9             10         11          12             13          14        15              16         17          18            19         20            21           22           23     24       
    titles =      {  'bgr0','bgr1',  'asym1',  'step_h',  'lorentz_width', 'intensity1',   'energy1',  'width1',  'intensity2',  'energy2',  'width2',  'intensity3',  'energy3',   'width3','intensity4',   'energy4',  'width4',  'intensity5',  'energy5',  'width5',  'intensity6',  'energy6',   'width6','asym2','bgr3','bgr4'};
        
    cs = c_data(:,1)';
    ubs = c_data(:,2)';
    lbs = c_data(:,3)';
    flags = logical(c_data(:,4))';
    
       
    ub = ubs(flags);
    lb = lbs(flags);
    c = cs(flags);
    energy_hold = energy;
    energy_hold(find(energy_hold(1:end-1)-energy(2:end)==0)) = [];
    e_range = 2*(max(energy)-min(energy));
    estep = abs(min(energy_hold(1:end-1)-energy_hold(2:end)));
    new_e = min(energy)-e_range:estep:max(energy)+e_range; 


    for bb = 1:length(energy)
        abs_newe = abs(new_e-energy(bb));
        index_e(bb) = find(abs_newe==min(abs_newe),1);
    end

    c_out = zeros(length(hv),length(cs));
    %PhD_data_square = (sum(PhD_data_square));
    err = zeros(no_peaks,length(hv));
    for nn = 1%:length(hv)
        c_out(nn,:) = cs;
        if nn == 1
            c = cs(flags);
            
            %c = [bgr0(1),intensity1(1)];%,intensity4(1),energy4(1),width4(1),intensity5(1),energy5(1),width5(1),asym1(1),asym2(1)];
            %c = [bgr0(1),bgr1(1),intensity1(1),intensity2(1),intensity3(1)];%,intensity4(1),energy4(1),width4(1),intensity5(1),energy5(1),width5(1),asym1(1),asym2(1)];
             test_fit(c,cs,flags,energy,squeeze((PhD_data_square(nn,:))),new_e,index_e,0);
            %pause(0.1)
        end
        [c_hold,RESNORM,RESIDUAL,EXITFLAG,OUTPUT,LAMBDA,jacobian] = lsqnonlin(@(c)  test_fit(c,cs,flags,energy,squeeze((PhD_data_square(nn,:))),new_e,index_e,0),c,lb,ub);
        c_out(nn,flags) = c_hold;
        %sum(sqrt(jacobian(:,3).^2))
        %whos jacobian
        %jacobian(:,3)
        %pause
        lb_jacob = lb;
        ub_jacob = ub;
        c_jacob = c_out(nn,:);
        
        res =  test_fit(c_hold,cs,flags,energy,squeeze((PhD_data_square(nn,:))),new_e,index_e,3);
        chisq = sum(res.^2);
        %(length(iin(nn,:))-length(c))
        
        
        for en = 1:no_peaks
            jacob = sqrt(sum(jacobian(:,en+nnz(flags)-no_peaks).^2));
            err(en,nn) = sqrt(chisq/(length(PhD_data_square(nn,:))-length(c_hold))/jacob*2);
        end       

        %[c_out(nn,3) err_3(nn)]
        %pause
        c = c_out(nn,flags);
        %c(2:4)
        %nn
        pause(0.01)
    end
  
    
    function [f] = test_fit(c,cs,flags,ein,iin,new_e,index_e,plotter)
        if plotter > 0
            figure(plotter)
            subplot(4,3,[5,6,8,9,11,12])
            hold off

            if nn == 1
                cla
            end
        end
        di = 6;
        %no_peaks = 1;
        [lineshape,di] = fit_peakn(c,cs,flags,ein,1,di,new_e,index_e,no_peaks,plotter);

        %di
        f = iin(:)-lineshape(:);
        if plotter > 0

            %figure(plotter)
            subplot(4,3,[5,6,8,9,11,12])
            set(gca,'Xdir','reverse')
            obj = (findobj(gca,'type','line'));
            if length(obj) < no_peaks + 2
                hold on            
                plot(ein,iin(:),'ro')
                plot(ein,lineshape,'k','Linewidth',2)
                hold off
            else
                set(obj(2),'Xdata',ein,'Ydata',iin(:));            
                set(obj(1),'Xdata',ein,'Ydata',lineshape); 
            end
            
            subplot(4,3,2:3)
            set(gca,'Xdir','reverse')
            obj = (findobj(gca,'type','line'));
            if isempty(obj)
                plot(ein,f)
                hold on
                plot(ein,zeros(size(ein)),'r--')
                hold off
                set(gca,'xtick',[])
                ylim([min(f),max(f)]*1.2)
            else
                set(obj(2),'Xdata',ein,'Ydata',f); 
            end
        end
        

    end   

end


