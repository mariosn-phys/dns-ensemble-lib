function [uLnk,vLnk,wLnk,gLnk]=compact_Lyap(uLn,vLn,wLn,gLn)

global kLyap

uLn_hat=fft(uLn,[],2);
vLn_hat=fft(vLn,[],2);
wLn_hat=fft(wLn,[],2);
gLn_hat=fft(gLn,[],2);

Lkiter=0;
for Lk=kLyap
    Lkiter=Lkiter+1;

    uLnk(:,:,:,Lkiter)=permute(uLn_hat(:,Lk+1,:,:,Lkiter),[1 3 4 5 2]);
    vLnk(:,:,:,Lkiter)=permute(vLn_hat(:,Lk+1,:,:,Lkiter),[1 3 4 5 2]);
    wLnk(:,:,:,Lkiter)=permute(wLn_hat(:,Lk+1,:,:,Lkiter),[1 3 4 5 2]);
    gLnk(:,:,:,Lkiter)=permute(gLn_hat(:,Lk+1,:,:,Lkiter),[1 3 4 5 2]);
    
end

end
