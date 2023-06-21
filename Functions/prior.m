function [px] = prior_uncertain_met(s)

% if x(1)>=-10&&x(1)<=60&&x(2)>=-20&&x(2)<=30&&x(3)>=500&&x(3)<=20000
%     px=1/70/50/20000;
% else
%     px=0;
% end

[m n]=size(s.x);%n=3
tempx=zeros(m,1);
tempx(s.Q>0&s.u>=0&s.ci>0&s.cii>0&s.z>0)=4;
px=find(tempx~=4);

end

