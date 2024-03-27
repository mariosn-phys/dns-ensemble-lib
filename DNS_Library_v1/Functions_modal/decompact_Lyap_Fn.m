function [uLn,vLn,wLn,gLn]=decompact_Lyap_Fn(uLnk,vLnk,wLnk,gLnk,Lkn)

global kLyap NLyap NX

fnn = (Lkn == kLyap');

NL = max(sum(fnn,2));
KL = length(kLyap);

uLn_hat=repmat(permute(0*uLnk(:,:,1),[1 5 2 3 4]),[1 NX 1 NL KL]);
vLn_hat=uLn_hat;
wLn_hat=uLn_hat;
gLn_hat=uLn_hat;

Lkiter=0;
for Lk=kLyap
    Lkiter=Lkiter+1;
    
    Ln = find( fnn(Lkiter,:) == 1);
    
    uLn_hat(:,Lk+1,:,:,Lkiter)=ipermute(uLnk(:,:,Ln(1):Ln(end)),[1 3 4 5 2]);
    vLn_hat(:,Lk+1,:,:,Lkiter)=ipermute(vLnk(:,:,Ln(1):Ln(end)),[1 3 4 5 2]);
    wLn_hat(:,Lk+1,:,:,Lkiter)=ipermute(wLnk(:,:,Ln(1):Ln(end)),[1 3 4 5 2]);
    gLn_hat(:,Lk+1,:,:,Lkiter)=ipermute(gLnk(:,:,Ln(1):Ln(end)),[1 3 4 5 2]);

    uLn_hat(:,NX-Lk+1,:,:,Lkiter)=conj(uLn_hat(:,Lk+1,:,:,Lkiter));
    vLn_hat(:,NX-Lk+1,:,:,Lkiter)=conj(vLn_hat(:,Lk+1,:,:,Lkiter));
    wLn_hat(:,NX-Lk+1,:,:,Lkiter)=conj(wLn_hat(:,Lk+1,:,:,Lkiter));
    gLn_hat(:,NX-Lk+1,:,:,Lkiter)=conj(gLn_hat(:,Lk+1,:,:,Lkiter));
    
end

uLn=real(ifft(uLn_hat,[],2));
vLn=real(ifft(vLn_hat,[],2));
wLn=real(ifft(wLn_hat,[],2));
gLn=real(ifft(gLn_hat,[],2));


end
