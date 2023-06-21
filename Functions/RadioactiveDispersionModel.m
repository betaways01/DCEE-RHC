function [ C ] = RadioactiveDispersionModel( s, p )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

s.D=s.ci;
s.t = s.cii;
epsilon = 1;
module_dist = sqrt(((s.x-p.x_matrix)).^2 + ((s.y-p.y_matrix)).^2);
%module_dist = sqrt(((s.x-p.x_matrix)).^2 + ((s.y-p.y_matrix)).^2 + ((s.z-p.z_matrix)).^2);

C = s.u.*s.t.*s.D.*(s.Q./(4*pi.*(module_dist.^2+epsilon))+s.phi);

end

