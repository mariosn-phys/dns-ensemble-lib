function [gi,vi,ui,wi] = solv_vgn_f_kron_zx_RK3(g0,v0,advg1,advv1,step,h,Fn)
% Solve C_N

global N NX MZ k l dt DYkron llm kkm ksolv msolv Nksolv Nmsolv
global ICvkron1 ICvDvkron1 ICgkron1 ICggkron1
global ICvkron2 ICvDvkron2 ICgkron2 ICggkron2
%global ICvkron3 ICvDvkron3 ICgkron3 ICggkron3

% ghat=fft2_cube(g0);g1v=reshape(permute(ghat(2:end-1,1:NX/3+1,[1:MZ/3+1 MZ-MZ/3+1:MZ]),[1 3 2]),[N*(NX/3+1)*(2*MZ/3+1),1]);
% vvhat=fft2_cube(v0);v1v=reshape(permute(vvhat(2:end-1,1:NX/3+1,[1:MZ/3+1 MZ-MZ/3+1:MZ]),[1 3 2]),[N*(NX/3+1)*(2*MZ/3+1),1]);
% ag1hat=fft2_cube(advg1)*(-dt)*h;ag1v=reshape(permute(ag1hat(2:end-1,1:NX/3+1,[1:MZ/3+1 MZ-MZ/3+1:MZ]),[1 3 2]),[N*(NX/3+1)*(2*MZ/3+1),1]);
% av1hat=fft2_cube(advv1)*(-dt)*h;av1v=reshape(permute(av1hat(2:end-1,1:NX/3+1,[1:MZ/3+1 MZ-MZ/3+1:MZ]),[1 3 2]),[N*(NX/3+1)*(2*MZ/3+1),1]);

krange=NX-flip(ksolv(2:end))+2;
mrange=[msolv MZ-flip(msolv(2:end))+2];

ghat=fft2_cuben(g0,Fn);g1v=reshape(permute(ghat(2:end-1,ksolv,mrange,:),[1 3 2 4]),[N*(Nksolv)*(Nmsolv),Fn]);
vvhat=fft2_cuben(v0,Fn);v1v=reshape(permute(vvhat(2:end-1,ksolv,mrange,:),[1 3 2 4]),[N*(Nksolv)*(Nmsolv),Fn]);
ag1hat=fft2_cuben(advg1,Fn)*(-dt)*h;ag1v=reshape(permute(ag1hat(2:end-1,ksolv,mrange,:),[1 3 2 4]),[N*(Nksolv)*(Nmsolv),Fn]);
av1hat=fft2_cuben(advv1,Fn)*(-dt)*h;av1v=reshape(permute(av1hat(2:end-1,ksolv,mrange,:),[1 3 2 4]),[N*(Nksolv)*(Nmsolv),Fn]);

vi=0*v0;
gi=vi;
ui=vi;
wi=vi;


% Calculate

if step==1
g1v=ICgkron1*ag1v+ICggkron1*g1v;
v1v=ICvDvkron1*av1v+ICvkron1*v1v;
elseif step==2 | step==3
g1v=ICgkron2*ag1v+ICggkron2*g1v;
v1v=ICvDvkron2*av1v+ICvkron2*v1v;   
end

tic;
dv1v=DYkron*v1v;
u1v=-llm*g1v+kkm*dv1v;
w1v=kkm*g1v+llm*dv1v;
toc

tic;
%dv1v=DYkron*v1v;
u1c=llmm*g1v+DYkronk*v1v;
w1c=kkm*g1v+DYkronl*v1v;
toc

for it=1:100
gputimeit(@() DYkron*v1v)
gputimeit(@() DYkronk*v1v)
end

tic;
%dv1v=DYkron*v1v;
u1c=-llm*g1v+DYkronk*v1v;
w1c=kkm*g1v+DYkronl*v1v;
toc

%Rebuild
gi(2:end-1,ksolv,mrange,:)=ipermute(reshape(g1v,[N Nmsolv Nksolv Fn]),[1 3 2 4]);
vi(2:end-1,ksolv,mrange,:)=ipermute(reshape(v1v,[N Nmsolv Nksolv Fn]),[1 3 2 4]);
ui(2:end-1,ksolv,mrange,:)=ipermute(reshape(u1v,[N Nmsolv Nksolv Fn]),[1 3 2 4]);
wi(2:end-1,ksolv,mrange,:)=ipermute(reshape(w1v,[N Nmsolv Nksolv Fn]),[1 3 2 4]);

gi(2:end-1,krange,mrange(2:end),:)=conj(flip(flip(gi(2:end-1,ksolv(2:end),mrange(2:end),:),3),2));
vi(2:end-1,krange,mrange(2:end),:)=conj(flip(flip(vi(2:end-1,ksolv(2:end),mrange(2:end),:),3),2));
ui(2:end-1,krange,mrange(2:end),:)=conj(flip(flip(ui(2:end-1,ksolv(2:end),mrange(2:end),:),3),2));
wi(2:end-1,krange,mrange(2:end),:)=conj(flip(flip(wi(2:end-1,ksolv(2:end),mrange(2:end),:),3),2));

gi(2:end-1,krange,1,:)=conj(flip(gi(2:end-1,ksolv(2:end),1,:),2));
vi(2:end-1,krange,1,:)=conj(flip(vi(2:end-1,ksolv(2:end),1,:),2));
ui(2:end-1,krange,1,:)=conj(flip(ui(2:end-1,ksolv(2:end),1,:),2));
wi(2:end-1,krange,1,:)=conj(flip(wi(2:end-1,ksolv(2:end),1,:),2));


end

% old code, NX/3+1 MZ/3+1
% 
% vi=zeros(N+2,NX,MZ);
% gi=zeros(N+2,NX,MZ);
% ui=zeros(N+2,NX,MZ);
% wi=zeros(N+2,NX,MZ);
% if step==1
% 
% gi(2:end-1,1:NX/3+1,[1:MZ/3+1 MZ-MZ/3+1:MZ])=ipermute(reshape(ICgkron1*ag1v+ICggkron1*g1v,[N 2*MZ/3+1 NX/3+1 ]),[1 3 2]);
% vi(2:end-1,1:NX/3+1,[1:MZ/3+1 MZ-MZ/3+1:MZ])=ipermute(reshape(ICvDvkron1*av1v+ICvkron1*v1v,[N 2*MZ/3+1 NX/3+1 ]),[1 3 2]);
% 
% elseif step==2
%     
% gi(2:end-1,1:NX/3+1,[1:MZ/3+1 MZ-MZ/3+1:MZ])=ipermute(reshape(ICgkron2*ag1v+ICggkron2*g1v,[N 2*MZ/3+1 NX/3+1 ]),[1 3 2]);
% vi(2:end-1,1:NX/3+1,[1:MZ/3+1 MZ-MZ/3+1:MZ])=ipermute(reshape(ICvDvkron2*av1v+ICvkron2*v1v,[N 2*MZ/3+1 NX/3+1 ]),[1 3 2]);
%     
% end
% 
% dydvihat=difY_F(vi,1);
% 
% % ghat=fft2_cube(g0);g1v=reshape(permute(ghat(2:end-1,1:NX/3+1,[1:MZ/3+1 MZ-MZ/3+1:MZ]),[1 3 2]),[N*(NX/3+1)*(2*MZ/3+1),1]);
% % vvhat=fft2_cube(v0);v1v=reshape(permute(vvhat(2:end-1,1:NX/3+1,[1:MZ/3+1 MZ-MZ/3+1:MZ]),[1 3 2]),[N*(NX/3+1)*(2*MZ/3+1),1]);
% % ag1hat=fft2_cube(advg1)*dt;ag1v=reshape(permute(ag1hat(2:end-1,1:NX/3+1,[1:MZ/3+1 MZ-MZ/3+1:MZ]),[1 3 2]),[N*(NX/3+1)*(2*MZ/3+1),1]);
% % av1hat=fft2_cube(advv1)*dt;av1v=reshape(permute(av1hat(2:end-1,1:NX/3+1,[1:MZ/3+1 MZ-MZ/3+1:MZ]),[1 3 2]),[N*(NX/3+1)*(2*MZ/3+1),1]);
% % 
% % vi=zeros(N+2,NX,MZ);
% % gi=zeros(N+2,NX,MZ);
% % ui=zeros(N+2,NX,MZ);
% % wi=zeros(N+2,NX,MZ);
% % 
% % gi(2:end-1,1:NX/3+1,[1:MZ/3+1 MZ-MZ/3+1:MZ])=ipermute(reshape(-ICgkron*ag1v+ICggkron*g1v,[N 2*MZ/3+1 NX/3+1 ]),[1 3 2]);
% % vi(2:end-1,1:NX/3+1,[1:MZ/3+1 MZ-MZ/3+1:MZ])=ipermute(reshape(-ICvDvkron*av1v+ICvkron*v1v,[N 2*MZ/3+1 NX/3+1 ]),[1 3 2]);
% % dydvihat=difY_F(vi,1);
% 
% for kk=1:NX/3+1
%     kp=k(kk);
%     
%     for jj=[1:MZ/3+1 MZ-MZ/3+1:MZ]
%         
%         lp=l(jj);
%         
%         if (kp==0) && (lp==0)
%         else
%             ui(:,kk,jj)=(-1i*lp*gi(:,kk,jj)+1i*kp*dydvihat(:,kk,jj))/(kp^2+lp^2);
%             wi(:,kk,jj)=(1i*kp*gi(:,kk,jj)+1i*lp*dydvihat(:,kk,jj))/(kp^2+lp^2);
%         end
%     end
% end
%
% for kk=2:NX/3+1
%     
%     gi(:,NX-kk+2,1)=conj(gi(:,kk,1));
%     vi(:,NX-kk+2,1)=conj(vi(:,kk,1));
%     ui(:,NX-kk+2,1)=conj(ui(:,kk,1));
%     wi(:,NX-kk+2,1)=conj(wi(:,kk,1));
%     
%     jj=[2:MZ/3+1 MZ-MZ/3+1:MZ];
%     gi(:,NX-kk+2,MZ-jj+2)=conj(gi(:,kk,jj));
%     vi(:,NX-kk+2,MZ-jj+2)=conj(vi(:,kk,jj));
%     ui(:,NX-kk+2,MZ-jj+2)=conj(ui(:,kk,jj));
%     wi(:,NX-kk+2,MZ-jj+2)=conj(wi(:,kk,jj));
%         
% end


