function f=Doniach_Sunjic(x,c)
%f = Doniach_Sunjic(x,c) where c needs 4 parameters [amplitude,asymmetry,
%centre of peak,width
%f = c(1)*cos(pi*c(2)/2+(1-c(2))*atan((x-c(3))/c(4)))./((c(4).^2+(x-c(3)).^2).^((1-c(2))/2));
x=x;
%c(3) = -c(3);
f = cos(pi*c(2)/2+(1-c(2))*atan((x-c(3))/c(4)))./((c(4).^2+(x-c(3)).^2).^((1-c(2))/2));
f = f./max(f)*c(1);
