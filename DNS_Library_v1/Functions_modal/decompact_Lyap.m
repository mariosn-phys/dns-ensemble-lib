function [uLn,vLn,wLn,gLn]=decompact_Lyap(uLnk,vLnk,wLnk,gLnk)

global kLyap NX

uLn_hat=repmat(permute(zeros(size(uLnk)),[1 5 2 3 4]),[1 NX 1 1 1]);
vLn_hat=uLn_hat;
wLn_hat=uLn_hat;
gLn_hat=uLn_hat;

Lkiter=0;
for Lk=kLyap
    Lkiter=Lkiter+1;
    
    uLn_hat(:,Lk+1,:,:,Lk)=ipermute(uLnk(:,:,:,Lkiter),[1 3 4 5 2]);
    vLn_hat(:,Lk+1,:,:,Lk)=ipermute(vLnk(:,:,:,Lkiter),[1 3 4 5 2]);
    wLn_hat(:,Lk+1,:,:,Lk)=ipermute(wLnk(:,:,:,Lkiter),[1 3 4 5 2]);
    gLn_hat(:,Lk+1,:,:,Lk)=ipermute(gLnk(:,:,:,Lkiter),[1 3 4 5 2]);

    uLn_hat(:,NX-Lk+1,:,:,Lk)=conj(uLn_hat(:,Lk+1,:,:,Lkiter));
    vLn_hat(:,NX-Lk+1,:,:,Lk)=conj(vLn_hat(:,Lk+1,:,:,Lkiter));
    wLn_hat(:,NX-Lk+1,:,:,Lk)=conj(wLn_hat(:,Lk+1,:,:,Lkiter));
    gLn_hat(:,NX-Lk+1,:,:,Lk)=conj(gLn_hat(:,Lk+1,:,:,Lkiter));
    
end

uLn=real(ifft(uLn_hat,[],2));
vLn=real(ifft(vLn_hat,[],2));
wLn=real(ifft(wLn_hat,[],2));
gLn=real(ifft(gLn_hat,[],2));


end
