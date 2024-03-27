function [Ptot] = project_adj_k(u1,v1,w1,u2,v2,w2)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

global M k dz dy NX L B kLyap yE

%    E1=Ener(u1,v1,w1);
%    E2=Ener(u2,v2,w2);

%for ix=0:.1:2*pi

%    cf=exp(1i*ix);
    
    uLhat=u1;
    vLhat=v1;
    wLhat=w1;

    uL1hat=u2;
    vL1hat=v2;
    wL1hat=w2;


%for kk=1:1%length(M)/2
    

    PK=trapz(yE,sum(conj(uLhat).*uL1hat+conj(vLhat).*vL1hat+conj(wLhat).*wL1hat,2),1)/L*dz/B/NX/NX;
    %PK=sum(sum(abs(uLhat(:,kk,:)).^2+abs(vLhat(:,kk,:)).^2+abs(wLhat(:,kk,:)).^2,3),1)*dy*dz/sqrt(E1*E2)
    
%end
Ptot=sum(PK);

%end

% kk=2:length(M)/2
% P2=sum(sum(conj(uLhat(:,kk,:)).*uL1hat(:,kk,:)+conj(vLhat(:,kk,:)).*vL1hat(:,kk,:)+conj(wLhat(:,kk,:)).*wL1hat(:,kk,:),3),1)*dy/L*dz/B/NX/NX


%P01=InnerP(u1,v1,w1,u2,v2,w2)/sqrt(E1*E2);


