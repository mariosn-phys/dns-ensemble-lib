function [uLnkn,vLnkn,wLnkn,gLnkn,Lkn]=compact_Lyap_Fn(uLn,vLn,wLn,gLn)

global kLyap NLyap

uLn_hat=fft(uLn,[],2);
vLn_hat=fft(vLn,[],2);
wLn_hat=fft(wLn,[],2);
gLn_hat=fft(gLn,[],2);

Lkn=[];

uLnkn=[];
vLnkn=[];
wLnkn=[];
gLnkn=[];

Lkiter=0;
MKiter=0;
for Lk=kLyap

    Lkiter=Lkiter+1;
    
%    Lkn = [Lkn kLyap(Lkiter)*ones(1,NLyap(Lkiter)/2)];
    Lkn = [Lkn kLyap(Lkiter)*ones(1,NLyap(Lkiter))];

    uLnkn = cat(3,uLnkn,permute(uLn_hat(:,Lk+1,:,1:NLyap(Lkiter),Lkiter),[1 3 4 5 2]));
    vLnkn = cat(3,vLnkn,permute(vLn_hat(:,Lk+1,:,1:NLyap(Lkiter),Lkiter),[1 3 4 5 2]));
    wLnkn = cat(3,wLnkn,permute(wLn_hat(:,Lk+1,:,1:NLyap(Lkiter),Lkiter),[1 3 4 5 2]));
    gLnkn = cat(3,gLnkn,permute(gLn_hat(:,Lk+1,:,1:NLyap(Lkiter),Lkiter),[1 3 4 5 2]));
    
%     uLnk(:,:,Lkiter)=permute(uLn_hat(:,Lk+1,:,:,Lkiter),[1 3 4 5 2]);
%     vLnk(:,:,Lkiter)=permute(vLn_hat(:,Lk+1,:,:,Lkiter),[1 3 4 5 2]);
%     wLnk(:,:,Lkiter)=permute(wLn_hat(:,Lk+1,:,:,Lkiter),[1 3 4 5 2]);
%     gLnk(:,:,Lkiter)=permute(gLn_hat(:,Lk+1,:,:,Lkiter),[1 3 4 5 2]);
     
end

end
