function [data_square,i0,hv,nixswr,angles,energy,polar] = i09_export_XSW(filenames,pathname)
% Read in XSW data and metadata from [filenames,pathname] (allows for multiselect)
% Outputs:
% data_square format: [detector energies (eV),photon energy (hv), binding energy (energy), Repeats (optional)]
% i0 = [?]
% hv = X-ray photon energy
% nixswr = rocking curve intensity measurement
% angles = [?]
% energy = measured binding energy of detected photo-electrons
% polar = polar angle of detector

%%
if isempty(filenames)
    [filenames,pathname] = uigetfile({'*.nxs','nexus'},'MultiSelect','on');
    %1[filename pathname] = uigetfile({'*.nxs','nexus'},'MultiSelect','on');
end

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
%name1 = hold_info.Groups.Name;
name2 = hold_info.Groups.Groups(1).Name;
name3 = hold_info.Groups.Groups(1).Datasets(1).Name;
%name3 = hold_info.Groups.Groups(2).Datasets(2).Name;
whole_name = [name2 '/' name3];
nixswr = [nixswr  h5read(filename,[name2 '/nixswr'])'];
data{aa} =  h5read(filename,[name2 '/image_data']);
energy = h5read(filename,[name2 '/energies']);
angles = h5read(filename,[name2 '/angles']);
%energy = [energy; h5read(filename,[name2 '/energies'])];
hv = [hv   h5read(filename,[name2 '/excitation_energy'])'];
%hv = [hv   h5read(filename,[name2 '/jenergy'])];
polar = mean(h5read(filename,[name2 '/smpmpolar']));

%i_drain = h5read(filename,[name2 '/bfmiamp19']);
try
i0 = [i0  h5read(filename,['/entry1/hm3amp20/hm3amp20'])];
catch
    i0 = [i0  h5read(filename,[name2 '/smpmiamp39'])];
end
%i0 = [i0  h5read(filename,[name2 '/hm3iamp20'])];
end
data_square = [];
%whos data_square
data_square_1 = [];
data_square_2 = [];
data_square_3 = [];
for aa = length(filenames):-1:1
    data_size = size(data{aa});
    data_vert = [];
    data_square_hold = {};
    data_square_hold_1 = zeros(data_size(end),1);
    data_square_hold_2 = zeros(data_size(end),1);
    data_square_hold_3 = zeros(data_size(end),1);
    data_hold = data{aa};
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
    %data_square(:,:,aa) = data_hold;
    data_square_hold = shiftdim(data_hold,1);
    data_square(:,:,:,aa) = data_square_hold;

    %size(data_square)
    %data_square_1 = [data_square_1 data_square_hold_1];
    %data_square_2 = [data_square_2 data_square_hold_2];
    %data_square_3 = [data_square_3 data_square_hold_3];
end
%license('inuse')

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

