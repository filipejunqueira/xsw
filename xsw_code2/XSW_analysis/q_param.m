function [beta,gamma,delta,Eb] = q_param(datadir,Z,n,l,ms)
% q_param(datadir,Z,n,l,ms)
cd_old = cd;
cd(datadir)
fid = fopen('q_param.txt');
if isstring(l)~=1

    switch l
        case 1
            l='s';
        case 2

            l='p';
        case 3
            l='d';
        case 4
            l='f';
    end
   
end

switch l
    case 's'
        js = '';
    case 'p'
    	if ms < 0
            js = '1/2';
        else
            js = '3/2';
        end
    case 'd'
    	if ms < 0
            js = '3/2';
        else
            js = '5/2';
        end
    case 'f'
    	if ms < 0
            js = '5/2';
        else
            js = '7/2';
        end        
end
name = [num2str(Z) ' ' num2str(n) l js];
line = [name(1:end-1) 'Z'];
nn = 0;
%line ~= name
while strcmp(line,name)~=1
    nn = nn +1;
    line = fgetl(fid);
    while length(line) < length(name)
        line = fgetl(fid);
        nn = nn+1;
    end
end
Eb = str2num(fgetl(fid));
Erow = str2num(fgetl(fid));
length(Erow);
sigma = [Erow;str2num(fgetl(fid))];
beta = [Erow;str2num(fgetl(fid))];
gamma = [Erow;str2num(fgetl(fid))];
delta = [Erow;str2num(fgetl(fid))];

cd(cd_old);
