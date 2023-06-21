function [ Wpnorm ] = normaliseWeight(theta,Wp, N)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
    pri = prior(theta);
    if length(pri) >= 1
        Wp(pri')=0;
    end
    Wp(isnan(Wp))=0;
    
    Wpnorm = Wp./sum(Wp);
     if sum(Wpnorm)==0
         Wpnorm = ones(N,1)/N;
    end

end

