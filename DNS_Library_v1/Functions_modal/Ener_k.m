function [Ptot]=Ener_k(u1,v1,w1,kk)

global M k dz dy NX L B kLyap yE

%    E1=Ener(u1,v1,w1);
%    E2=Ener(u2,v2,w2);

    uLhat=u1;
    vLhat=v1;
    wLhat=w1;


%for kk=2:kLyap+1%length(M)/2
    

    PK=trapz(yE,sum(conj(uLhat).*uLhat+conj(vLhat).*vLhat+conj(wLhat).*wLhat,2),1)*dz/L/B/NX/NX;
    %PK=sum(sum(abs(uLhat(:,kk,:)).^2+abs(vLhat(:,kk,:)).^2+abs(wLhat(:,kk,:)).^2,3),1)*dy*dz/sqrt(E1*E2)
    
%end


Ptot=real(PK);