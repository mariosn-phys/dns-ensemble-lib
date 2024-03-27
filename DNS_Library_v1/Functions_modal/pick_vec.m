function [uLo,vLo,wLo,gLo] = pick_vec(uL,vL,wL,gL,Ln,Lk)

        uLo=uL(:,:,:,Ln,Lk);
        vLo=vL(:,:,:,Ln,Lk);
        wLo=wL(:,:,:,Ln,Lk);
        gLo=gL(:,:,:,Ln,Lk);

end

