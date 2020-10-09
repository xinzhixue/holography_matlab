function [M] = propogate(E, wl, Nx, Ny, dx, z, opt)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

x =((1:Nx)-(Nx+1)/2)/(Nx*dx);
y =((1:Ny)-(Ny+1)/2)/(Ny*dx);
[X,Y] = meshgrid(x,y);

d = 1-(wl.*X).^2-(wl.*Y).^2;
n = exp(-1i*2*pi*z/wl*d(d>0).^0.5);
n = reshape(n, [Ny, Nx]);

if strcmp(opt, 'shift')
    E_fft = fftshift(fft2(E));
elseif strcmp(opt, 'unshift')
    E_fft = (fft2(E));
else
    error('select: shift or unshift');
end

M = ifft2(E_fft.*n(1:Ny,1:Nx));
end

