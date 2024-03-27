function [up,vp,wp,gp] = adj_norm_k(up,vp,wp,gp,scf)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

%        Ef_p=Ener_k(up,vp,wp);

%        scf = scf/ sqrt(Ef_p);

        up=scf*up;
        vp=scf*vp;
        wp=scf*wp;
        gp=scf*gp;


end

