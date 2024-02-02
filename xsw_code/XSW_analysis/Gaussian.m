function f=Gaussian(x,c)
%Gaussian(x,c): creates a Gaussian function f(x) with an amplitude c(1), its 
%centre at c(2) and the a width of c(3)
f = c(1)*exp(-((x-c(2)).^2)/(2*c(3)^2));
