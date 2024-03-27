function [Ptot] = project_Lyap(u1,v1,w1,u2,v2,w2,Lk)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

global M k dz dy NX L B kLyap yE

%    E1=Ener(u1,v1,w1);
%    E2=Ener(u2,v2,w2);

    uLhat=fft(u1,[],2);
    vLhat=fft(v1,[],2);
    wLhat=fft(w1,[],2);

    uL1hat=fft(u2,[],2);
    vL1hat=fft(v2,[],2);
    wL1hat=fft(w2,[],2);


    kk=Lk;%length(M)/2
    

    PK=trapz(yE,sum(conj(uLhat(:,kk,:)).*uL1hat(:,kk,:)+conj(vLhat(:,kk,:)).*vL1hat(:,kk,:)+conj(wLhat(:,kk,:)).*wL1hat(:,kk,:),3),1)/L*dz/B/NX/NX;
    %PK=sum(sum(abs(uLhat(:,kk,:)).^2+abs(vLhat(:,kk,:)).^2+abs(wLhat(:,kk,:)).^2,3),1)*dy*dz/sqrt(E1*E2)
    



Ptot=real(sum(PK));

% kk=2:length(M)/2
% P2=sum(sum(conj(uLhat(:,kk,:)).*uL1hat(:,kk,:)+conj(vLhat(:,kk,:)).*vL1hat(:,kk,:)+conj(wLhat(:,kk,:)).*wL1hat(:,kk,:),3),1)*dy/L*dz/B/NX/NX


%P01=InnerP(u1,v1,w1,u2,v2,w2)/sqrt(E1*E2);


