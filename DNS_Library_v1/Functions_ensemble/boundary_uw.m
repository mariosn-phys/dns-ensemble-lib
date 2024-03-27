function [u1,w1] = boundary_uw(u1,w1,U1,W1,Fn)

global b1 b2 NX MZ 
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

    u1=repmat(permute(U1,[1 3 4 2]),[1,NX,MZ,1])+u1;u1=[ones(1,NX,MZ,Fn)*b1;u1(2:end-1,:,:,:);ones(1,NX,MZ,Fn)*b2];
    w1=repmat(permute(W1,[1 3 4 2]),[1,NX,MZ,1])+w1;w1=[zeros(1,NX,MZ,Fn);w1(2:end-1,:,:,:);zeros(1,NX,MZ,Fn)];

end