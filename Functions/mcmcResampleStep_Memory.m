function [ theta ] = mcmcResampleStep_Memory( theta,Wpnorm,N, D_k_store,m, dk, hopt,pos,P_k_store,PF_Memory )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

keep = [];
for zz = 1:1%5
    keep = [];
    n.x = theta.x + (hopt.*dk.x.*randn(N,1));
    n.y = theta.y + (hopt.*dk.y.*randn(N,1));
    n.z = theta.z ;%+ (hopt.*dk.z.*randn(N,1));
    n.Q = theta.Q + (hopt.*dk.Q.*randn(N,1));
    n.u = theta.u ;%+ (hopt.*dk.u.*randn(N,1));
    n.phi = theta.phi ;%+ (hopt.*dk.phi.*randn(N,1));
    n.ci = theta.ci ;%+ (hopt.*dk.ci.*randn(N,1));
    n.cii = theta.cii ;%+ (hopt.*dk.cii.*randn(N,1));
    for loops = 1:2
        pri = prior(n);
        numPri = length(pri);
        n.x(pri) = theta.x(pri) + (hopt.*dk.x.*randn(numPri,1));
        n.y(pri) = theta.y(pri) + (hopt.*dk.y.*randn(numPri,1));
        n.z(pri) = theta.z(pri) ;%+ (hopt.*dk.z.*randn(numPri,1));
        n.Q(pri) = theta.Q(pri) + (hopt.*dk.Q.*randn(numPri,1));
        n.u(pri) = theta.u(pri) ;%+ (hopt.*dk.u.*randn(numPri,1));
        n.phi(pri) = theta.phi(pri) ;%+ (hopt.*dk.phi.*randn(numPri,1));
        n.ci(pri) = theta.ci(pri) ;%+ (hopt.*dk.ci.*randn(numPri,1));
        n.cii(pri) = theta.cii(pri) ;%+ (hopt.*dk.cii.*randn(numPri,1));
    end
    pri = prior(n);
    n.x(pri) = theta.x(pri);
    n.y(pri) = theta.y(pri);
    n.z(pri) = theta.z(pri);
    n.Q(pri) = theta.Q(pri);
    n.u(pri) = theta.u(pri);
    n.phi(pri) = theta.phi(pri);
    n.ci(pri) = theta.ci(pri);
    n.cii(pri) = theta.cii(pri);
    
    
    nWpnorm = Wpnorm;
    r = size(D_k_store);
    if PF_Memory==1 || length(P_k_store(1,:))==1
        pos.X = P_k_store(end,1);
        pos.Y = P_k_store(end,2);
        D = D_k_store(end);
        
        nC = RadioactiveDispersionModel(n,pos);
        
        
        % nWp = Likelihood_Plain_Gaussian(nC,D,Wpnorm);
        nWp = Likelihood_Like_Yee(nC, D, Wpnorm, m);
        
        nWpnorm = normaliseWeight(n,nWp,N);
    else
        for mem = 1:PF_Memory
            % Estimeated from particle filter
            % nC = simpleGaussianPlume(n,m,pos);
            ind = ceil(rand*r(1));
            pos.X = P_k_store(ind,1);
            pos.Y = P_k_store(ind,2);
            D = D_k_store(ind);
            
            nC = RadioactiveDispersionModel(n,pos);
            
            
            % nWp = Likelihood_Plain_Gaussian(nC,D,Wpnorm);
            nWp = Likelihood_Like_Yee(nC, D, nWpnorm, m);
            
            nWpnorm = normaliseWeight(n,nWp,N);
        end
    end
    alpha = nWpnorm./(Wpnorm);%(Wpnorm(indx))
    %         alpha = nWpnorm/(1/N);
    
    mcrand = rand(N,1);
    keep = find(alpha>mcrand);
    % keep = 1:N;
    notkeep = find(alpha<rand);
    
    theta.x(keep) = n.x(keep);
    theta.y(keep) = n.y(keep);
    theta.z(keep) = n.z(keep);
    theta.Q(keep) = n.Q(keep);
    theta.u(keep) = n.u(keep);
    theta.phi(keep) = n.phi(keep);
    theta.ci(keep) = n.ci(keep);
    theta.cii(keep) = n.cii(keep);
    
    
    Wpnorm = ones(N,1)/N;
    
end
end

