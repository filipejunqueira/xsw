function [out] = step_fcn(x,c)
%[out] = step_fcn(x,c) - c(1) - amplitude, c(2) - position, c(3) - width


out =c(1) * (erf((x-c(2))/c(3))+1)/2;