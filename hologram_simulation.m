% function hologram_simulation(filename, resolution, wavelength, refractiveindex)
%HOLOGRAM_SIMULATION Simulate Hologram using Rayleigh-Sommerfeld kernel
%   Input:
%   filename: a .dat file specifing the location and size of droplets
%   resolution: the resolution of the x-y plane, in m/px
%   wavelength: the wavelength of the incident beam
%   refractiveindex: the refractive index of the medium
%   Output:
%   a simulated hologram

%% Defalut values
clear
clc
close all
resolution = 4.03e-6;
wavelength = 635e-9;        
refractiveindex = 1;
%% Define parameters
dx = resolution;        % The Resolution in meter
wl = wavelength;        % The Wavelength
nm = refractiveindex;   % surrounding medium refractive indices

%% Read in data file for objects
Nx = 1024; Ny = Nx;

% define centers in px
%! fix this with data file read in
obj.centers = [512, 512;
    482, 150;
    532, 150];
obj.z = [40e-3;
    10e-3;
    10e-3];% temp in m
obj.radii = [40;
    20;
    20];

%% sort objects by z depth
obj.all = [obj.z, obj.centers, obj.radii];
obj.sorted = sortrows(obj.all);

[obj.unique_z, id] = unique(obj.sorted(:,1));

mask = zeros(Ny, Nx, length(obj.unique_z));
for i = 1:length(obj.unique_z)
    if id(i) ~= size(obj.sorted, 1) && length(obj.unique_z) ~= 1
        mask(:,:,i) = generate_circularmask([Ny, Nx], obj.sorted(id(i):id(i+1)-1, 2:3), obj.sorted(id(i):id(i+1)-1, 4));    
    else
        mask(:,:,i) = generate_circularmask([Ny, Nx], obj.sorted(id(i), 2:3), obj.sorted(id(i), 4));
    end
%     figure, imshow(mask(:,:,i))    
end
mask = 1-double(mask);
%% Propogate the plane wave
E = mask(:, :, end); % initial value
i = size(mask, 3);
z = obj.unique_z(i)-obj.unique_z(i-1);

M = propogate(E, wl, Nx, Ny, dx, z, 'shift');
E = M.*mask(:, :, i-1);

for i = size(mask, 3)-1:-1:2
    z = obj.unique_z(i)-obj.unique_z(i-1);
    M = propogate(E, wl, Nx, Ny, dx, z, 'unshift');
    E = M.*mask(:, :, i-1);
end
z = obj.unique_z(1);
M = propogate(E, wl, Nx, Ny, dx, z, 'unshift');

hologram = rescale(abs(M));

figure, imshow(prod(mask, 3))

figure,imshow(hologram);
axis equal
