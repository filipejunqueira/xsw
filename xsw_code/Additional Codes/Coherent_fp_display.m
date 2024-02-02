close all;clc;clear

% Coherent fraction (f_val) and position (p_val) values and errors
% Example - Post anneal TIPB
f_val = [0.43,0.05,0.17]; 
f_err = [0.06,0.03,0.03];
p_val = [0.3,0.84,0.85]; 
p_err = [0.04,0.05,0.03];

% Plots errorbars in polar coordinates for each f and p value
for n = 1:size(f_val,2)
    subplot(1,size(f_val,2),n)
    
    h(n,:) = polarerror(p_val(n)*2*pi,f_val(n),p_err(n)*2*pi,f_err(n),'r');
    hold all

    text(135/180*pi,1.5,sprintf('f = %0.2f+/-%0.2f \np = %0.2f+/-%0.2f',f_val(n),f_err(n),p_val(n),p_err(n)))
    
end


% Formats polar plots
AX = findobj(gcf,'type','polaraxes');
numax = length(AX); %number of axes
for n = 1:numax
    set(AX(n),'fontsize',16)
    subplot(1,3,n)    
    rlim([0,1])
    thetaticks([0 90 180 270])
    thetaticklabels({'0.00','0.25','0.50','0.75'})
    rticks([0 0.25 0.50 0.75 1.00])
    rticklabels({'0','0.25','0.50','0.75','1.00'})
    AX(n).RAxisLocation = 90;
end
