function [c_out,titles,err_1,err_2] = XPS_Base(data_square,hv,energy)

    PhD_data_square = data_square;

  %background_mult = [10000   1000000 1000];
   % energy1         = [400, 405,        395];
   bgr0             = [3,      1e6,   0.0005       1];
   bgr1             = [0.0367, 2e5,    -2e5            1];
   
   asym1            = [0.040826    0.4     1e-5    1];
   
   lorentz_width    = [0.4   1.5    0.1    1];
   
   step_h           = [0.02   0.9     0.00  1];
  
   energy1          = [285.6,  285, 286,    1];
   intensity1       = [10     1e7     0       1];
   width1           = [0.48412   0.9     0.2   1];
   
   energy2          = [2.5,   3.5,  1.5   0];
   intensity2       = [250    1e6     1e-8    0];
   width2           = [0.48765   8  0.2    0];
   
    
    energy3         = [0.5,    0.7,   0.4    0];
    intensity3      = [70     1000     0       0];
    width3          = [0.75    1   0.2    0];    
    
    energy4         = [3.6,     6,   3     0];
    intensity4      = [50     1e6     0       0];
    width4          = [3    8    0.1  0]; 
    
    energy5         = [11.3691,      11.5,   11.1         0];
    intensity5      = [700     1e6     0      0];
    width5          = [5.6991     9           0.1        0]; 
    
    energy6         = [8.9397,      10,     8    0];
    intensity6      = [500      1e6     0         0];
    width6          = [0.6330        0.9    0.1    0];     
    
     asym2 =           [0.45     0.5    0.3    0];
     
    bgr3         =    [0,      1,   0.000000005       0];   
    bgr4         =    [0, 2e-0,    -2e-0            0];      
    

    cc=1;
    %                   1       2        3         4              5    %         6              7           8           9             10         11          12             13          14        15              16         17          18            19         20            21           22           23     24       
    titles =      {  'bgr0','bgr1',  'asym1',  'step_h',  'lorentz_width', 'intensity1',   'energy1',  'width1',  'intensity2',  'energy2',  'width2',  'intensity3',  'energy3',   'width3','intensity4',   'energy4',  'width4',  'intensity5',  'energy5',  'width5',  'intensity6',  'energy6',   'width6','asym2','bgr3','bgr4'};
        
    cs =            [bgr0(cc),bgr1(cc),asym1(cc),step_h(cc),lorentz_width(cc),intensity1(cc),energy1(cc),width1(cc),intensity2(cc),energy2(cc),width2(cc),intensity3(cc),energy3(cc),width3(cc),intensity4(cc),energy4(cc),width4(cc),intensity5(cc),energy5(cc),width5(cc),intensity6(cc),energy6(cc),width6(cc),asym2(cc),bgr3(cc),bgr4(cc)];
    cc=2;
    ubs =           [bgr0(cc),bgr1(cc),asym1(cc),step_h(cc),lorentz_width(cc),intensity1(cc),energy1(cc),width1(cc),intensity2(cc),energy2(cc),width2(cc),intensity3(cc),energy3(cc),width3(cc),intensity4(cc),energy4(cc),width4(cc),intensity5(cc),energy5(cc),width5(cc),intensity6(cc),energy6(cc),width6(cc),asym2(cc),bgr3(cc),bgr4(cc)];
    cc=3;
    lbs =           [bgr0(cc),bgr1(cc),asym1(cc),step_h(cc),lorentz_width(cc),intensity1(cc),energy1(cc),width1(cc),intensity2(cc),energy2(cc),width2(cc),intensity3(cc),energy3(cc),width3(cc),intensity4(cc),energy4(cc),width4(cc),intensity5(cc),energy5(cc),width5(cc),intensity6(cc),energy6(cc),width6(cc),asym2(cc),bgr3(cc),bgr4(cc)];
    cc=4;
    flags = logical([bgr0(cc),bgr1(cc),asym1(cc),step_h(cc),lorentz_width(cc),intensity1(cc),energy1(cc),width1(cc),intensity2(cc),energy2(cc),width2(cc),intensity3(cc),energy3(cc),width3(cc),intensity4(cc),energy4(cc),width4(cc),intensity5(cc),energy5(cc),width5(cc),intensity6(cc),energy6(cc),width6(cc),asym2(cc),bgr3(cc),bgr4(cc)]);
    %flags
    ub = ubs(flags);
    lb = lbs(flags);
    c = cs(flags);
    energy_hold = energy;
    energy_hold(find(energy_hold(1:end-1)-energy(2:end)==0)) = [];
    e_range = 2*(max(energy)-min(energy));
    estep = abs(min(energy_hold(1:end-1)-energy_hold(2:end)));
    new_e = min(energy)-e_range:estep:max(energy)+e_range; 
%     figure(5)
%     plot(energy)
%     hold on
%     plot(new_e)
%     whos new_e
%     hold off
%     pause
    for bb = 1:length(energy)
        abs_newe = abs(new_e-energy(bb));
        index_e(bb) = find(abs_newe==min(abs_newe),1);
    end
    

    find(ub<lb)
    c_out = zeros(length(hv),length(cs));
    %PhD_data_square = (sum(PhD_data_square));
    for nn = 1:length(hv)
        c_out(nn,:) = cs;
        if nn == 1
            c = cs(flags);
            
            %c = [bgr0(1),intensity1(1)];%,intensity4(1),energy4(1),width4(1),intensity5(1),energy5(1),width5(1),asym1(1),asym2(1)];
            %c = [bgr0(1),bgr1(1),intensity1(1),intensity2(1),intensity3(1)];%,intensity4(1),energy4(1),width4(1),intensity5(1),energy5(1),width5(1),asym1(1),asym2(1)];
             test_fit(c,cs,flags,energy,squeeze((PhD_data_square(nn,:))),new_e,index_e,1);
            pause(0.1)
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
        jacob = sqrt(sum(jacobian(:,3).^2));
        err_1(nn) = sqrt(chisq/(length(PhD_data_square(nn,:))-length(c_hold))/jacob*2);
        
        
%         jacob = sqrt(sum(jacobian(:,4).^2));
%         err_2(nn) = sqrt(chisq/(length(PhD_data_square(nn,:))-length(c_hold))/jacob*2);
        
%         jacob = sqrt(sum(jacobian(:,5).^2));
%         err_3(nn) = sqrt(chisq/(length(PhD_data_square(nn,:))-length(c_hold))/jacob*2);      
         
%          jacob = sqrt(sum(jacobian(:,6).^2));
%          err_4(nn) = sqrt(chisq/(length(PhD_data_square(nn,:))-length(c_hold))/jacob*2);  
         
%          jacob = sqrt(sum(jacobian(:,7).^2));
%          err_5(nn) = sqrt(chisq/(length(PhD_data_square(nn,:))-length(c_hold))/jacob*2);  
%          
%          jacob = sqrt(sum(jacobian(:,8).^2));
%          err_6(nn) = sqrt(chisq/(length(PhD_data_square(nn,:))-length(c_hold))/jacob*2);  
         
        %[c_out(nn,3) err_3(nn)]
        %pause
        c = c_out(nn,flags);
        %c(2:4)
        %nn
        pause(0.2)
    end
    
    
    function [f] = test_fit(c,cs,flags,ein,iin,new_e,index_e,plotter)
        if plotter > 0
            cla
            figure(plotter)
            hold off
            set(gca,'Xdir','reverse')

            %pause
        end
        di = 6;
        no_peaks = 1;
        [lineshape,di] = fit_peakn(c,cs,flags,ein,1,di,new_e,index_e,no_peaks,plotter);

        %di
        if plotter > 0
            figure(plotter)
            %size(ein)
            %size(iin)
            hold on
            plot(ein,iin(:),'ro')         
            plot(ein,lineshape,'k','Linewidth',2)
            
            hold off
        end
        f = iin(:)-lineshape(:);

    end

    
end

