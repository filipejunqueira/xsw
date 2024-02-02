function [data_square,i0,hv,nixswr,angles,energy,polar] = i09_export_XPS1(filenames,pathname)
angles = [];
polar = [];

%% XPS Data
% Reads and plots XPS data
%% %%%%%%%%%%%%%%%%%%%%%%% Read Data%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear data_square_hold
%[filenames,pathname] = uigetfile({'*.nxs','nexus'},'MultiSelect','on');
%[filename pathname] = uigetfile({'*.nxs','nexus'},'MultiSelect','on');
%filenames = [['i09-11268.nxs']; ['i09-11269.nxs']; ['i09-11270.nxs']];
cd(pathname)
data = [];
hv = [];
i0 = [];
nixswr = [];
energy = [];
if iscell(filenames)~= 1
    filenames = {filenames};
end%filenames = {filename};
for aa = length(filenames):-1:1
    filename = filenames{aa};
%cd(pathname)
hold_info = h5info(filename);
try 
%name1 = hold_info.Groups.Name;
    name2 = hold_info.Groups.Groups(1).Name;
    %name2 = hold_info.Groups.Groups(2).Name;       %N1s
    name3 = hold_info.Groups.Groups(1).Datasets(1).Name;
    name4 = hold_info.Groups.Groups(2).Name;
    name5 = hold_info.Groups.Groups(3).Name;
    name6 = hold_info.Groups.Groups(4).Name;
    name7 = hold_info.Groups.Groups(6).Name;
    name8 = hold_info.Groups.Groups(7).Name;
catch
    name9 = 'no'
end
%name10 = hold_info.Groups.Groups(8).Name;


%name3 = hold_info.Groups.Groups(2).Datasets(2).Name;
whole_name = [name2 '/' name3];
%hv = [hv   h5read(filename,[name2 '/excitation_energy'])];


% Carbon
data{aa} =  h5read(filename,[name2 '/image_data']);
energy = h5read(filename,[name2 '/energies']);

% Nitrogen
% data{aa} =  h5read(filename,[name4 '/image_data']);%
% energy = h5read(filename,[name4 '/energies']);

position{aa} = [h5read(filename,[name2 '/smpmx']);h5read(filename,[name2 '/smpmy']);h5read(filename,[name2 '/smpmz']);];
disp(position{aa})
% data3{aa} =  h5read(filename,[name5 '/image_data']);
% energy3 = h5read(filename,[name5 '/energies']);
% 
% data4{aa} =  h5read(filename,[name6 '/image_data']);
% energy4 = h5read(filename,[name6 '/energies']);
% 
% data5{aa} =  h5read(filename,[name7 '/image_data']);
% energy5 = h5read(filename,[name7 '/energies']);
% 
% data6{aa} =  h5read(filename,[name8 '/image_data']);
% energy6 = h5read(filename,[name8 '/energies']);



%energy = [energy; h5read(filename,[name2 '/energies'])];

%hv = [hv   h5read(filename,[name2 '/jenergy'])];
%nixswr = [nixswr  h5read(filename,[name2 '/nixswr'])];
%i_drain = h5read(filename,[name2 '/bfmiamp19']);
%i0 = [i0  h5read(filename,['/entry1/hm3amp20/hm3amp20'])];
%i0 = [i0  h5read(filename,[name2 '/sm5iamp8'])];
end
data_square = [];
%whos data_square
data_square_1 = [];
data_square_2 = [];
data_square_3 = [];
for aa = length(filenames):-1:1
    data_size = size(data{aa});
    data_vert = [];
    data_square_hold = [];
    data_square_hold_1 = zeros(data_size(end),1);
    data_square_hold_2 = zeros(data_size(end),1);
    data_square_hold_3 = zeros(data_size(end),1);
    data_hold = data{aa};
%     data_hold2 = data2{aa};
%     data_hold3 = data3{aa};
    %for nn = 1:data_size(end)
    %    for oo = 1:data_size(1)
    %        for zz = 1:data_size(2)
            %data_square_hold(oo) = sum(data_hold(:,oo));
            %data_square_hold(oo) = data_hold(oo,1);
            %data_square_hold(oo) = sum(data_hold(:,1,1,oo));
            %data_square_hold_1(oo) = sum(data_hold(oo,1:100));
            %data_square_hold_2(oo) = sum(data_hold(oo,101:200));
            %data_square_hold_3(oo) = sum(data_hold(oo,201:300));
            %data_square_hold(oo,nn) = sum(data_hold(oo,:,1,nn));
    %           data_square_hold(zz,oo,nn) = data_hold(oo,zz,1,nn);
    %        end
    %    end
    %end
    %whos data_square
    %size(data_square)
    %if isempty(data_square) == 1
    %    data_square = data_square_hold;
    %else
    %    size_data_square = size(data_square);
    %    size_data_hold = size(data_square_hold);
    %    data_square(:,:,size_data_square(3):size_data_square(3)+size_data_hold(3)-1) = data_square_hold;
    %end
    data_square(aa,:) = shiftdim(data_hold,1);
%     data_square2 = shiftdim(data_hold2,1);
%     data_square3 = shiftdim(data_hold3,1);
    %size(data_square)
    %data_square_1 = [data_square_1 data_square_hold_1];
    %data_square_2 = [data_square_2 data_square_hold_2];
    %data_square_3 = [data_square_3 data_square_hold_3];
end


%energy(269) = [];
%energy(180) = [];
%energy(91) = [];

%energy(59) = [];
%energy(40) = [];
%energy(21) = [];
%energy(2) = [];
%data_square(269,:) = [];
%data_square(180,:) = [];
%data_square(91,:) = [];
%data_square(59,:) = [];
%data_square(40,:) = [];
%data_square(21,:) = [];
%data_square(2,:) = [];
%for nn = 1:data_size(end)
%   data_vert = [data_vert;energy data_square(:,nn)];
%end

%1:30 80:end


%% %%%%%%%%%%%%%%%%%%%%%%% Plotting %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% figure(2); clf; hold all; xlabel('Binding Energy /eV'); ylabel('Intensity'); title('Sum over hv')
% for n = 1:size(data_square,1) + 1
%     if n <= size(data_square,1)
%         plot(energy,data_square(n,:),'.')
%     else
%         plot(energy,mean(data_square,1))
%     end
% end

%%
% hv = 435;
% [c_out,titles,err_1,err_2] = XPS_C1s_DITP(data_square,hv,energy);
% for a =1:14
%     disp(['n = ',num2str(a),', ',titles{a},' Average = ',num2str(mean(c_out(:,a)))])
% end

hv = 435;
data_export = data_square;





