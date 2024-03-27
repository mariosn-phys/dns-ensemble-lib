function [advU,advW] = advect_mean(up,vp,wp,Fn)

global gamma DYF
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

    gn=repmat(gamma,[1 Fn]);

    if Fn==1
    advU=DYF*mean(mean(up.*vp,3),2)+gn;
    advW=DYF*mean(mean(wp.*vp,3),2);
    else
    advU=DYF*permute(mean(mean(up.*vp,3),2),[1 4 2 3])+gn;
    advW=DYF*permute(mean(mean(wp.*vp,3),2),[1 4 2 3]);
    end

end