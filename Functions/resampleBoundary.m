function [theta] = resampleBoundary(theta,Wp,N,domain)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
j=0;
for i=1:N
        if inpolygon(theta.x(i),theta.y(i),domain(1:2),domain(3:4))
            j=j+1;
            theta.x(j)=theta.x(i);
            theta.y(j)=theta.y(i);
            %theta.z(j)=theta.z(i);
            theta.Q(j)=theta.Q(i);
            %theta.u(j)=theta.u(i);
            %theta.phi(j)=theta.phi(i);
            %theta.ci(j)=theta.ci(i);
            %theta.cii(j)=theta.cii(i);
            Wp(j)=Wp(i);
        end
end
j=j+1;
theta.x(j:end)=[];
theta.y(j:end)=[];
%theta.z(j:end)=[];
theta.Q(j:end)=[];
%theta.u(j:end)=[];
%theta.phi(j:end)=[];
%theta.ci(j:end)=[];
%theta.cii(j:end)=[];
Wp(j:end)=[];

indx=resampleStratified(Wp,N);

theta.x = theta.x(indx);
theta.y = theta.y(indx);
%theta.z = theta.z(indx);
theta.Q = theta.Q(indx);
%theta.u = theta.u(indx);
%theta.phi = theta.phi(indx);
%theta.ci = theta.ci(indx);
%theta.cii = theta.cii(indx);

end

