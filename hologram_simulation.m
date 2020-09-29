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

x =((1:Nx)-(Nx+1)/2)/(Nx*dx);
y =((1:Ny)-(Ny+1)/2)/(Ny*dx);
[X,Y] = meshgrid(x,y);

% define centers in px
%! fix this with data file read in
obj.centers = [512, 512;
    482, 150;
    532, 150];
obj.z = [20e-3;
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
    mask = double(imcomplement(mask));
%     figure, imshow(mask(:,:,i))    
end

%% Propogate plane wave
M = mask(:, :, end);
for i = size(mask, 3):-1:1         
    mask_fft = fftshift(fft2(M));
    if i > 1
        z = obj.unique_z(i)-obj.unique_z(i-1);
    elseif i == 1
        z = obj.unique_z(i);
    end
        
    d = 1-(wl.*X).^2-(wl.*Y).^2;
    n = exp(-1i*2*pi*z/wl*d(d>0).^0.5);
    n = reshape(n, [Ny, Nx]);
        
    M1 = ifft2(mask_fft.*n(1:Ny,1:Nx));
    if i ~= 1
        M = M1.*mask(:, :, i-1);
    elseif i == 1
        M = M1;
    end
end
        
hologram = rescale(abs(M));

figure, imshow(sum(mask, 3))

figure,imagesc(hologram);
axis equal

% figure,plot(hologram(512,:));
% hold on
% plot(hologram(150,:));
% ylim([0, 0.1]);
    % imwrite(hologram, 'sim_hologram1.tif')
    
% end
% end

