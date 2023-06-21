function [ C ] = simpleGaussianPlume( s, m, p )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
% Output value of concentration at points given by a matrix produced by a
% source defined by s

ang = atan2((p.x_matrix-s.x),(p.y_matrix-s.y));
ang = ang-pi/2;
ang_diff = (s.phi)+ang;
dist = ((p.x_matrix-s.x).^2+(p.y_matrix-s.y).^2).^0.5;
xr = dist .* cos(ang_diff);
yr = dist .* sin(ang_diff);
zr = p.z_matrix; % Assume
% sy = s.ci.*xr.^p.b;
% sz = s.cii.*xr.^p.d + p.f;
sy = s.ci.*xr./((1+0.0001.*xr).^0.5);
sz = s.cii.*xr./((1+0.0003.*xr).^0.5);
sy(sy<0)=0;
sz(sz<0)=0;
sy(imag(sy) ~= 0) = 0;
sz(imag(sz) ~= 0) = 0;
C = (s.Q./(s.u.*sy.*sz*2*pi)).*exp(-(yr.^2)./(2*sy.^2)).*(exp(-(zr-s.z).^2./(2*sz.^2))+exp(-(zr+s.z).^2./(2*sz.^2)));
C(isnan(C))=0;




end

