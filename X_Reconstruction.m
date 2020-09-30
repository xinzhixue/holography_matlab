% Reconstruct Hologram using Kirchhoff-Fresnel or Rayleigh-Sommerfeld

% =========================================================================
% Revise Log (only recent):
% 2015/04/14--Add image cropping part to test square vs. rectangular
% image reconstruction difference on edge. Theoretically it should be the
% same; Test result shows similar ringing effect on the image edge.
% 2015/04/23--Add mean intensity padding around original image.
% 2015/08/12--Modify mean intensity padding to zero padding, failed to
% reduce ringing effect, change to mean intensity padding
% =========================================================================
clear
clc

tic;
%%%% parameter
dx=4.03e-6;         %% The Resolution in meter
WL=635e-9;        %% The Wavelength
Nm=1;                 %% surrounding medium refractive indices
zSeq=(18:0.05:22)/Nm*10^(-3);     %% reconstruction depth, this distance is measured in air
k=2*pi*Nm/WL;              %% wave number

I_raw = imread('sim_hologram1.tif');
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
%% Reconstruction using Kirchhoff-Fresnel
% I1=double(I1);
% I_fft = fftshift(fft2(I1));
% for i=1:size(zSeq,2)
%     z = zSeq(i);
%     for r = 1:Nx
%         n(r,1:Ny) = (WL*j)*exp(-j*pi*WL*z*(x(r)^2 + y.^2));
%     end
%     M1 = I_fft.*n(1:Ny,1:Nx);
%     M2 = ifft2(M1);
%     M(:,:,i) = abs(M2);
% end
%% Reconstruction using Rayleigh-Sommerfeld
I1=double(I1);
I = fftshift(fft2(I1));
for m=1:size(zSeq,2)
    z = zSeq(m);
    for r = 1:Nx
        d(1:Ny)=1- (WL*x(r))^2 - (WL*y).^2;
        id=find(d>0);
        n(1:Ny, r)=0;
        n(id, r)=exp(-1i*2*pi*z/WL*d(id).^0.5);
    end
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


