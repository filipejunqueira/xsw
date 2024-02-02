function [f,p] = XSW_bodger(hv,c_out,nixswr,err_1,bragg,polar,f,p,plane,name)
clf
%f = 0.5; p = 0.5;
q = [1,f,p];
q_scale = 0;

f_text = uicontrol('style','edit',...
    'units','normalized',...
    'position',[0.05,0.8,0.2,0.05],...
    'fontsize',30,...
    'string',sprintf('f = %0.2f',f),...
    'callback',@f_text_func);


f_slider = uicontrol('style','slider',...
    'units','normalized',...
    'position',[0.05,0.75,0.2,0.05],...
    'min',0,'max',1,'value',f,...
    'sliderstep',[0.01,0.1],...
    'callback',@f_func);




p_text = uicontrol('style','edit',...
    'units','normalized',...
    'position',[0.05,0.65,0.2,0.05],...
    'fontsize',30,...
    'string',sprintf('p = %0.2f',p),...
    'callback',@p_text_func);


        
p_slider = uicontrol('style','slider',...
    'units','normalized',...
    'position',[0.05,0.6,0.2,0.05],...
    'min',0,'max',1,'value',p,...
    'sliderstep',[0.01,0.1],...
    'callback',@p_func);

i_text = uicontrol('style','text',...
    'units','normalized',...
    'position',[0.05,0.5,0.2,0.05],...
    'fontsize',20,...
    'string',name);



errorbutton = uicontrol('style','togglebutton',...
    'units','normalized',...
    'position',[0.05,0.87,0.095,0.1],...
    'string','Toggle error',...
    'value',1,...
    'fontsize',20,...
    'callback',@err_func);

linebutton = uicontrol('style','togglebutton',...
    'units','normalized',...
    'position',[0.155,0.87,0.095,0.1],...
    'string','Toggle line',...
    'value',1,...
    'fontsize',20,...
    'callback',@line_func);

    pause(0.1)
    q_scale = XSW_Bodge('C:\Users\ppxcj\Documents\Experiments\XSW\XSW_codes\My_code\fpfpp',[hv',c_out,nixswr',err_1'],15,q,bragg,polar,plane,q_scale);


function f_text_func(source,event)
    f = str2num(f_text.String(end-3:end));
    set(f_slider,'Value',f)
    q = [1,f,p];
    q_scale = XSW_Bodge('C:\Users\ppxcj\Documents\Experiments\XSW\XSW_codes\My_code\fpfpp',[hv',c_out,nixswr',err_1'],15,q,bragg,polar,plane,q_scale);
end


function p_text_func(source,event)
    p = str2num(p_text.String(end-3:end));
    set(p_slider,'Value',p)
    q = [1,f,p];
    q_scale = XSW_Bodge('C:\Users\ppxcj\Documents\Experiments\XSW\XSW_codes\My_code\fpfpp',[hv',c_out,nixswr',err_1'],15,q,bragg,polar,plane,q_scale);
end



function f_func(source,event);
    f = f_slider.Value;
    set(f_text,'string',sprintf('f = %0.2f',f))
    q = [1,f,p];
    q_scale = XSW_Bodge('C:\Users\ppxcj\Documents\Experiments\XSW\XSW_codes\My_code\fpfpp',[hv',c_out,nixswr',err_1'],15,q,bragg,polar,plane,q_scale);
end

function p_func(source,event)
    p = p_slider.Value;
    set(p_text,'string',sprintf('p = %0.2f',p))
    q = [1,f,p];
    q_scale = XSW_Bodge('C:\Users\ppxcj\Documents\Experiments\XSW\XSW_codes\My_code\fpfpp',[hv',c_out,nixswr',err_1'],15,q,bragg,polar,plane,q_scale);
end


% for f = 0.1:0.1:1
%     for p = 0.1:0.1:1
%         cla
%         q = [1,f,p];
%         XSW_Bodge('C:\Users\ppxcj\Documents\Experiments\XSW\XSW_codes\My_code\fpfpp',[hv',c_out,nixswr',err_1'],15,q);
%         set(fp,'string',sprintf('f = %0.1f q = %0.1f',f,p))
%         pause(0.2)
%         
%     end
% end
function err_func(source,event)
    ebars = findobj(gcf,'type','errorbar');
    if errorbutton.Value == 0
        set(ebars,'visible','off')
    else
        set(ebars,'visible','on')
    end
end

function line_func(source,event)
    line_no = findobj(gcf,'type','line');
    if linebutton.Value == 1
        set(line_no(2),'linestyle','none')
    else
        set(line_no(2),'linestyle','-')
    end
end


end