function [gi,vi,ui,wi] = solv_vg_lyap_RK_kron_v2_k(g0,v0,advg1,advv1,step,h)
% Solve C_N

global N NX MZ M P k l D2 I ICvRK ICgRK dt CgRK CvRK DYF NLyap kLyap DYkronm kkmm llmm msolv Nmsolv

global ICvCvRK1 ICvIDELvRK1 ICgRK1 ICgCgRK1
global ICvCvRK2 ICvIDELvRK2 ICgRK2 ICgCgRK2

    NL=sum(NLyap)/2;

    mrange=[msolv MZ-flip(msolv(2:end))+2];

    ghat=fft(g0,[],2);g1v=reshape(ghat(2:end-1,mrange,:),[N*Nmsolv*NL,1]);
    vvhat=fft(v0,[],2);v1v=reshape(vvhat(2:end-1,mrange,:),[N*Nmsolv*NL,1]);
    ag1hat=fft(advg1,[],2)*(-dt)*h;ag1v=reshape(ag1hat(2:end-1,mrange,:),[N*Nmsolv*NL,1]);
    av1hat=fft(advv1,[],2)*(-dt)*h;av1v=reshape(av1hat(2:end-1,mrange,:),[N*Nmsolv*NL,1]);
    
    
        vi=0 * g0;
        gi=vi;
        ui=vi;
        wi=vi;

%sind=((kL-1)*N*(2*MZ/3+1)+1):(kL*N*(2*MZ/3+1));  
        
if step==1

g1v=ICgRK1*ag1v+ICgCgRK1*g1v;
v1v=ICvIDELvRK1*av1v+ICvCvRK1*v1v;

elseif step==2
    
g1v=ICgRK2*ag1v+ICgCgRK2*g1v;
v1v=ICvIDELvRK2*av1v+ICvCvRK2*v1v;
    
end


dv1v=DYkronm*v1v;
u1v=-llmm*g1v+kkmm*dv1v;
w1v=kkmm*g1v+llmm*dv1v;
    

%Rebuild
gi(2:end-1,mrange,:)=reshape(g1v,[N Nmsolv NL]);
vi(2:end-1,mrange,:)=reshape(v1v,[N Nmsolv NL]);
ui(2:end-1,mrange,:)=reshape(u1v,[N Nmsolv NL]);
wi(2:end-1,mrange,:)=reshape(w1v,[N Nmsolv NL]);
               
end

%  old code
%         vi=zeros(N+2,MZ,NL);
%         gi=zeros(N+2,MZ,NL);
%         ui=zeros(N+2,MZ,NL);
%         wi=zeros(N+2,MZ,NL);
% 
% %sind=((kL-1)*N*(2*MZ/3+1)+1):(kL*N*(2*MZ/3+1));  
%         
% if step==1
% 
% gi(2:end-1,[1:MZ/3+1 MZ-MZ/3+1:MZ],:)=reshape(ICgRK1*ag1v+ICgCgRK1*g1v,[N 2*MZ/3+1 NL]);
% vi(2:end-1,[1:MZ/3+1 MZ-MZ/3+1:MZ],:)=reshape(ICvIDELvRK1*av1v+ICvCvRK1*v1v,[N 2*MZ/3+1 NL]);
% 
% elseif step==2
%     
% gi(2:end-1,[1:MZ/3+1 MZ-MZ/3+1:MZ],:)=reshape(ICgRK2*ag1v+ICgCgRK2*g1v,[N 2*MZ/3+1 NL]);
% vi(2:end-1,[1:MZ/3+1 MZ-MZ/3+1:MZ],:)=reshape(ICvIDELvRK2*av1v+ICvCvRK2*v1v,[N 2*MZ/3+1 NL]);
%     
% end
% 
%                
% Miter=0;
% Lkiter=0;
% for Lk=kLyap
%     Lkiter=Lkiter+1;
%     
%         for Ln=1:2:NLyap(Lkiter)
%         Miter=Miter+1;   
%         
%         dydvihat=DYF*vi(:,:,Miter);        
%         
%         kp=k(Lk+1);
%        
%         for jj=[1:MZ/3+1 MZ-MZ/3+1:MZ]
%             
%         lp=l(jj);       
%         if (kp==0) && (lp==0)                     
%         else
%         ui(:,jj,Miter)=(-1i*lp*gi(:,jj,Miter)+1i*kp*dydvihat(:,jj))/(kp^2+lp^2);
%         wi(:,jj,Miter)=(1i*kp*gi(:,jj,Miter)+1i*lp*dydvihat(:,jj))/(kp^2+lp^2);
%         end
% 
%         end     
%         
%         end
%         
% end
