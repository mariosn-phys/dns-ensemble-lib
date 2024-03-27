function [u0,v0,w0,g0] = make_uw_ens_Lyap(vi,gi,UP1,WP1,b1,b2,N_ens,K_ens)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

for ik=1:K_ens
for ij=1:N_ens(ik)
[u0(:,:,:,ij,ik),v0(:,:,:,ij,ik),w0(:,:,:,ij,ik),g0(:,:,:,ij,ik)] = make_uw(vi(:,:,:,ij,ik),gi(:,:,:,ij,ik),UP1,WP1,b1,b2);
end
end

end

