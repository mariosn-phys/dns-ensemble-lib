function [u0,v0,w0,g0]=noise_short_kron_en_kx(FFkron,F,sc,kx,kz,Fn,u0)

% noise from complete optimal base
% total energy scaling

global N NX MZ DYkronf llmf kkmf

vi=u0;
gi=vi;
ui=vi;
wi=vi;

Nk = length(kx);
Nm = length([kz MZ-flip(kz(2:end))+2]);

nn=FFkron*(randn(F*Nm*Nk,Fn)+1i*randn(F*Nm*Nk,Fn))/sqrt(2);

nn=reshape(nn,[F,Nm,Nk,Fn]); 

giv = reshape(nn(N+1:2*N,:,:,:),[N*Nk*Nm,Fn]);
viv = reshape(nn(1:N,:,:,:),[N*Nk*Nm,Fn]);

dviv= DYkronf*viv;
uiv = -llmf*giv+kkmf*dviv; 
wiv = kkmf*giv +llmf*dviv;

%gi=[zeros(1,NX,MZ);nn(N+1:2*N,:,:);zeros(1,NX,MZ)];
%vi=[zeros(1,NX,MZ);nn(1:N,:,:);zeros(1,NX,MZ)];

%%% 1:NX/3 NX-NX/3+2:NX   == 1 kx(2:end) NX-flip(kx(2:end))+2
%%% 1:MZ/3 MZ-MZ/3+2:MZ   == 1 kz(2:end) MZ-flip(kz(2:end))+2


%Rebuild
gi(2:end-1,kx,[kz MZ-flip(kz(2:end))+2:MZ],:)=ipermute(reshape(giv,[N Nm Nk Fn]),[1 3 2 4]);
vi(2:end-1,kx,[kz MZ-flip(kz(2:end))+2:MZ],:)=ipermute(reshape(viv,[N Nm Nk Fn]),[1 3 2 4]);
ui(2:end-1,kx,[kz MZ-flip(kz(2:end))+2:MZ],:)=ipermute(reshape(uiv,[N Nm Nk Fn]),[1 3 2 4]);
wi(2:end-1,kx,[kz MZ-flip(kz(2:end))+2:MZ],:)=ipermute(reshape(wiv,[N Nm Nk Fn]),[1 3 2 4]);

gi(2:end-1,NX-flip(kx(2:end))+2,[kz MZ-flip(kz(2:end))+2:MZ],:)=conj(flip(flip(gi(2:end-1,kx(2:end),[kz MZ-flip(kz(2:end))+2:MZ],:),3),2));
vi(2:end-1,NX-flip(kx(2:end))+2,[kz MZ-flip(kz(2:end))+2:MZ],:)=conj(flip(flip(vi(2:end-1,kx(2:end),[kz MZ-flip(kz(2:end))+2:MZ],:),3),2));
ui(2:end-1,NX-flip(kx(2:end))+2,[kz MZ-flip(kz(2:end))+2:MZ],:)=conj(flip(flip(ui(2:end-1,kx(2:end),[kz MZ-flip(kz(2:end))+2:MZ],:),3),2));
wi(2:end-1,NX-flip(kx(2:end))+2,[kz MZ-flip(kz(2:end))+2:MZ],:)=conj(flip(flip(wi(2:end-1,kx(2:end),[kz MZ-flip(kz(2:end))+2:MZ],:),3),2));

gi(2:end-1,NX-flip(kx(2:end))+2,1,:)=conj(flip(gi(2:end-1,kx(2:end),1,:),2));
vi(2:end-1,NX-flip(kx(2:end))+2,1,:)=conj(flip(vi(2:end-1,kx(2:end),1,:),2));
ui(2:end-1,NX-flip(kx(2:end))+2,1,:)=conj(flip(ui(2:end-1,kx(2:end),1,:),2));
wi(2:end-1,NX-flip(kx(2:end))+2,1,:)=conj(flip(wi(2:end-1,kx(2:end),1,:),2));

v0=ifft2_cuben(vi,Fn);
u0=ifft2_cuben(ui,Fn);
w0=ifft2_cuben(wi,Fn);
g0=ifft2_cuben(gi,Fn);

u0=[zeros(1,NX,MZ,Fn);u0(2:end-1,:,:,:);zeros(1,NX,MZ,Fn)];
w0=[zeros(1,NX,MZ,Fn);w0(2:end-1,:,:,:);zeros(1,NX,MZ,Fn)];

% Ei=Ener(u0,v0,w0);
% kap=sc/(Ei^(1/2));
% u0=kap*u0;
% v0=kap*v0;
% w0=kap*w0;
% g0=kap*g0;

end

