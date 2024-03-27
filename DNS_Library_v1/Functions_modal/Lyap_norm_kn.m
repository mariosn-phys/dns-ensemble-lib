function [up,vp,wp,gp] = Lyap_norm_kn(up,vp,wp,gp,scf,Ef_p)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

%        Ef_p=Ener_k(up,vp,wp);

        up=scf*up/sqrt(Ef_p);
        vp=scf*vp/sqrt(Ef_p);
        wp=scf*wp/sqrt(Ef_p);
        gp=scf*gp/sqrt(Ef_p);


end

