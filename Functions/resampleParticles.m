function [ theta dk hopt ] = resampleParticles( theta,Wpnorm,N )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
%   indx = resampleStratified(Wpnorm);
    indx = resampleSystematic(Wpnorm);

theta.x = theta.x(indx);
theta.y = theta.y(indx);
%theta.z = theta.z(indx);
theta.Q = theta.Q(indx);
%theta.u = theta.u(indx);
%theta.phi = theta.phi(indx);
%theta.ci = theta.ci(indx);
%theta.cii = theta.cii(indx);

Covx = cov(theta.x);
Covy = cov(theta.y);
%Covz = cov(theta.z);
CovQ = cov(theta.Q);
%Covu = cov(theta.u);
%Covphi = cov(theta.phi);
%Covci = cov(theta.ci);
%Covcii = cov(theta.cii);


        dk.x = cholcov(Covx);%sqrt(CovXxp(1));
        dk.y = cholcov(Covy);%sqrt(CovXyp(1));
        %dk.z = cholcov(Covz);%sqrt(CovXyp(1));
        dk.Q = cholcov(CovQ);%sqrt(CovXqp(1));
        %dk.u = cholcov(Covu);%sqrt(CovXqp(1));
        %dk.phi = cholcov(Covphi);%sqrt(CovXqp(1));
        %dk.ci = cholcov(Covci);
        %dk.cii = cholcov(Covcii);
        
        mm=3; % The dimension of parameter space
        A=(4/(mm+2))^(1/(mm+4));%InstantaneousGaussian
        % Need to make many changes if dimension of parameter space changes
        cx = 4*pi/3;
%         A = (8/cx*(m+4)*((2*sqrt(pi))^m))^(1/(m+4));
        hopt=A*(N^(-1/(mm+4)));
        


end

