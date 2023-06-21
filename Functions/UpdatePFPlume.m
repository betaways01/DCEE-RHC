function [ theta, Wpnorm  ] = UpdatePFPlume( D_k_store, theta, Wpnorm, pos,P_k_store, m,N,PF_Memory,domain )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

% Estimated from particle filter
    % C = simpleGaussianPlume(theta,m,pos);
    C = RadioactiveDispersionModel(theta,pos);
    % C(C<m.thresh)=0;
    
    % Wp = Likelihood_Plain_Gaussian(C, D, Wpnorm);
    Wp = Likelihood_Like_Yee(C, D_k_store(end), Wpnorm,m);
    keep=inpolygon(theta.x,theta.y,[domain(1),domain(2), domain(2), domain(1)],[domain(3),domain(3),domain(4),domain(4)]);
    Wp(~keep)=0;
    Wpnorm = normaliseWeight(theta,Wp,N);
    Neff = 1/sum(Wpnorm.^2);
    
    

%     xmean = sum(Wpnorm.*theta.x);ymean = sum(Wpnorm.*theta.y);zmean = sum(Wpnorm.*theta.z);
%     qmean = sum(Wpnorm.*theta.Q);umean = sum(Wpnorm.*theta.u);phimean = sum(Wpnorm.*theta.phi);
%     ci = sum(Wpnorm.*theta.ci);cii = sum(Wpnorm.*theta.cii);
    
    % Resample step
    if Neff < 0.9*N% for low release rate or sensing area
%         resample = resample +1; 
%         indx = resampleStratified(Wpnorm);
        indx = resampleSystematic(Wpnorm);
        [theta, dk, hopt] = resampleParticles(theta,Wpnorm,N);
        theta = mcmcResampleStep_Memory(theta,Wpnorm,N, D_k_store,m, dk, hopt,pos,P_k_store,PF_Memory);
        Wpnorm = ones(N,1)/N;
    end

end

