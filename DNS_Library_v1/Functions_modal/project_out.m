function [uL1,vL1,wL1,gL1] = project_out(u1,v1,w1,g1,u2,v2,w2,g2,P)
%PROJECT_OUT Summary of this function goes here
%   Detailed explanation goes here
global M k dz dy NX L B kLyap



%     uLhat=fft(u1,[],2);
%     vLhat=fft(v1,[],2);
%     wLhat=fft(w1,[],2);
%     gLhat=fft(g1,[],2);
% 
%     uL1hat=fft(u2,[],2);
%     vL1hat=fft(v2,[],2);
%     wL1hat=fft(w2,[],2);
%     gL1hat=fft(g2,[],2);

    uL1=u2-P.*u1;
    vL1=v2-P.*v1;
    wL1=w2-P.*w1;
    gL1=g2-P.*g1;

%   %for kk=[1 2:kLyap+1 length(M)-kLyap+1:length(M)]%:length(M)
% for kk=1:length(M)
% 
%     if kk<=length(M)/2
%         P1=P;
%     elseif kk>=length(M)/2+1
%         P1=conj(P);
%     end
%     
%     %kp=k(kk);
%     uL2hat(:,kk,:)=uL1hat(:,kk,:)-P1*uLhat(:,kk,:);
%     vL2hat(:,kk,:)=vL1hat(:,kk,:)-P1*vLhat(:,kk,:);
%     wL2hat(:,kk,:)=wL1hat(:,kk,:)-P1*wLhat(:,kk,:);
%     gL2hat(:,kk,:)=gL1hat(:,kk,:)-P1*gLhat(:,kk,:);
% 
%     %PK(kk)=sum(sum(conj(uLhat(:,kk,:)).*uL1hat(:,kk,:)+conj(vLhat(:,kk,:)).*vL1hat(:,kk,:)+conj(wLhat(:,kk,:)).*wL1hat(:,kk,:),3),1)*dy/L*dz/B/sqrt(E1*E2)/NX/NX
%     %PK=sum(sum(abs(uLhat(:,kk,:)).^2+abs(vLhat(:,kk,:)).^2+abs(wLhat(:,kk,:)).^2,3),1)*dy*dz/sqrt(E1*E2)
% end
  
    
% for kk=[1 2:kLyap+1 length(M)-kLyap+1:length(M)]%:length(M)
% 
%     if kk<=length(M)/2
%         P1=P;
%     elseif kk>=length(M)/2+1
%         P1=conj(P);
%     end
%     
%     %kp=k(kk);
%     uL1hat(:,kk,:)=uL1hat(:,kk,:)-P1*uLhat(:,kk,:);
%     vL1hat(:,kk,:)=vL1hat(:,kk,:)-P1*vLhat(:,kk,:);
%     wL1hat(:,kk,:)=wL1hat(:,kk,:)-P1*wLhat(:,kk,:);
%     gL1hat(:,kk,:)=gL1hat(:,kk,:)-P1*gLhat(:,kk,:);
% 
%     %PK(kk)=sum(sum(conj(uLhat(:,kk,:)).*uL1hat(:,kk,:)+conj(vLhat(:,kk,:)).*vL1hat(:,kk,:)+conj(wLhat(:,kk,:)).*wL1hat(:,kk,:),3),1)*dy/L*dz/B/sqrt(E1*E2)/NX/NX
%     %PK=sum(sum(abs(uLhat(:,kk,:)).^2+abs(vLhat(:,kk,:)).^2+abs(wLhat(:,kk,:)).^2,3),1)*dy*dz/sqrt(E1*E2)
% end
% 
% 
% % norm(uL2hat(:)-uL1hat(:))
% % norm(vL2hat(:)-vL1hat(:))
% % norm(wL2hat(:)-wL1hat(:))
% % norm(gL2hat(:)-gL1hat(:))
% 
% %uL1=real(ifft(uL1hat,[],2));
% %vL1=real(ifft(vL1hat,[],2));
% %wL1=real(ifft(wL1hat,[],2));
% %gL1=real(ifft(gL1hat,[],2));
% 
% uL1=ifft(uL1hat,[],2);
% vL1=ifft(vL1hat,[],2);
% wL1=ifft(wL1hat,[],2);
% gL1=ifft(gL1hat,[],2);
% 
% % norm(uL1(:)-uL2(:))
% % norm(vL1(:)-vL2(:))
% % norm(wL1(:)-wL2(:))
% % norm(gL1(:)-gL2(:))
% 
% 
% 
% %[up,vp,wp,gp] = Lyap_norm(uL1,vp,wp,gp,scf)

end

