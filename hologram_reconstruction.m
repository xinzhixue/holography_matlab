function hologram_reconstruction(filename, resolution, wavelength, ...
    refractiveindex, z_depth)
%HOLOGRAM_RECONSTRUCTION Reconstruct Hologram using Rayleigh-Sommerfeld
%   Input:
%   filename: a .dat file specifing the location and size of droplets
%   resolution: the resolution of the x-y plane, in m/px
%   wavelength: the wavelength of the incident beam
%   refractiveindex: the refractive index of the medium
%   z_depth: an 1D array representing the reconstruction volume, relative 
%   to the focal plane, in mm
%   Output:
%   A folder "Rec/" with reconstructed holograms under current directory


tic;
%%%% parameter
dx = resolution;         %% The Resolution in meter
wl = wavelength;        %% The Wavelength
Nm = refractiveindex;                 %% surrounding medium refractive indices
zSeq = z_depth/Nm*10^(-3);     %% reconstruction depth, this distance is measured in air
k=2*pi*Nm/wl;              %% wave number

I_raw = imread(filename);
I1 = I_raw;
% I1 = imcrop(I_raw, [1, 400, 1499, 1499]);
% Add mean intensity padding around original image
% I_big = mean(I1(:))*ones(size(I1, 1)+20, size(I1, 2)+20);
% I1 = I_big;

%     % Add mean padding around original image
% I_big = padarray(I1,[500 500],mean(I1(:)));
% I1 = I_big;

Nx=size(I1,2);             %% Image size
Ny=size(I1,1);
x =((1:Nx)-(Nx+1)/2)/(Nx*dx);
y =((1:Ny)-(Ny+1)/2)/(Ny*dx);
[X,Y] = meshgrid(x,y);
%%%%%%%%%%%%%%%%%%%%%%%%%% Reconstruction using fast fourier transform
%% Reconstruction using Kirchhoff-Fresnel
% I1=double(I1);
% I_fft = fftshift(fft2(I1));
% for i=1:size(zSeq,2)
%     z = zSeq(i);
%     n = (WL*1i)*exp(-1i*pi*WL*z*(X.^2+Y.^2));
%     M1 = I_fft.*n;
%     M2 = ifft2(M1);
%     M(:,:,i) = abs(M2);
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Reconstruction using Rayleigh-Sommerfeld
I1=double(I1);
I = fftshift(fft2(I1));
for m=1:size(zSeq,2)
    z = zSeq(m);
    
    d = 1-(wl.*X).^2-(wl.*Y).^2;
    n = zeros(Ny, Nx);
    n = exp(-1i*2*pi*z/wl*d(d>0).^0.5);
    n = reshape(n, [Ny, Nx]);
    
    M1 = I.*n(1:Ny,1:Nx);
    M2 = ifft2(M1);
    M(:,:,m) = abs(M2);
end
%% Writing the file on PC
if isdir('Rec')
    rmdir('Rec', 's');
end
direc = 'Rec/';
mkdir(direc);
M=(M-min(M(:)))/(max(M(:))-min(M(:)));
% M=1-M;
% MInten=0.86; %% mean intensity of the reconstruction volume
for n=1:size(zSeq,2)
    I_rec=M(:,:,n);
    %     Diff=0.86-mean(I(:));
    %     I=I+Diff; %% rescale to the mean intensity, so keep each plane have same mean intensity
    imwrite(I_rec,sprintf([direc,'Rec_%3.3d.tif'],n),'tif','compression','none');
end
toc;



end

