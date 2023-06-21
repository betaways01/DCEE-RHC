function [ Wp ] = Likelihood_Like_Yee(C, D, Wpnorm, m)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

PH0 = 0.3;
PH1 = 1-PH0;
ProbBackground= 1;

NDsigma = m.thresh;%1e-6;
NDsigma = 1e-4+C;

% sigma = 5e-5+5e-3*exp(-dist/10);%0.06;%
sigma = 1e-4+C;
% sigma = m.thresh+C;

if D<=m.thresh
    Wp = Wpnorm .* ((PH0*ProbBackground) + (PH1*(1/2)*(1+erf((m.thresh-C)./(NDsigma.*sqrt(2))))));
    Wp(isnan(Wp))=0;
else
    Wp = Wpnorm .* 1./sigma.*sqrt(2*pi).*exp((-15*(C - D).^(2))./(sigma.^2));
    Wp(isnan(Wp))=0;
end    

% sigma = sigmaFun(C);
% 
% if D==0
%     prob= (1/2)*(1+erf((m.thresh-C)./(sigma.*sqrt(2))));
%    Wp = Wpnorm .* prob;
%    length(find(prob==1));
% else
%     prob = 1./sigma.*sqrt(2*pi).*exp((-(C - D).^2)./(sigma.^2));
%     length(find(prob<=0));
% Wp = Wpnorm .* prob;
% Wp(isnan(Wp))=0;
% end   

end

