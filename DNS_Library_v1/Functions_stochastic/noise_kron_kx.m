function [u0,v0,w0,g0]=noise_kron_kx(FFkron,F,sc,kx,kz)

% noise from complete optimal base
% total energy scaling

global N NX MZ DYkron llm kkm

vi=zeros(N+2,NX,MZ);
gi=zeros(N+2,NX,MZ);
ui=zeros(N+2,NX,MZ);
wi=zeros(N+2,NX,MZ);

nn=FFkron*(randn(F*(2*MZ/3-1)*NX/3,1)+1i*randn(F*(2*MZ/3-1)*NX/3,1))/sqrt(2);

nn=reshape(nn,[F,(2*MZ/3-1),NX/3]); 

giv = reshape(nn(N+1:2*N,:,:),N*NX/3*(2*MZ/3-1));
viv = reshape(nn(1:N,:,:),N*NX/3*(2*MZ/3-1));

dviv= DYkron*viv;
uiv = -llm*giv+kkm*dviv; 
wiv = kkm*giv +llm*dviv;

%gi=[zeros(1,NX,MZ);nn(N+1:2*N,:,:);zeros(1,NX,MZ)];
%vi=[zeros(1,NX,MZ);nn(1:N,:,:);zeros(1,NX,MZ)];

%Rebuild
gi(2:end-1,1:NX/3,[1:MZ/3 MZ-MZ/3+2:MZ])=ipermute(reshape(giv,[N 2*MZ/3-1 NX/3]),[1 3 2]);
vi(2:end-1,1:NX/3,[1:MZ/3 MZ-MZ/3+2:MZ])=ipermute(reshape(viv,[N 2*MZ/3-1 NX/3]),[1 3 2]);
ui(2:end-1,1:NX/3,[1:MZ/3 MZ-MZ/3+2:MZ])=ipermute(reshape(uiv,[N 2*MZ/3-1 NX/3]),[1 3 2]);
wi(2:end-1,1:NX/3,[1:MZ/3 MZ-MZ/3+2:MZ])=ipermute(reshape(wiv,[N 2*MZ/3-1 NX/3]),[1 3 2]);

gi(2:end-1,NX-NX/3+2:NX,[2:MZ/3 MZ-MZ/3+2:MZ])=conj(flip(flip(gi(2:end-1,2:NX/3,[2:MZ/3  MZ-MZ/3+2:MZ]),3),2));
vi(2:end-1,NX-NX/3+2:NX,[2:MZ/3 MZ-MZ/3+2:MZ])=conj(flip(flip(vi(2:end-1,2:NX/3,[2:MZ/3  MZ-MZ/3+2:MZ]),3),2));
ui(2:end-1,NX-NX/3+2:NX,[2:MZ/3 MZ-MZ/3+2:MZ])=conj(flip(flip(ui(2:end-1,2:NX/3,[2:MZ/3  MZ-MZ/3+2:MZ]),3),2));
wi(2:end-1,NX-NX/3+2:NX,[2:MZ/3 MZ-MZ/3+2:MZ])=conj(flip(flip(wi(2:end-1,2:NX/3,[2:MZ/3  MZ-MZ/3+2:MZ]),3),2));

gi(2:end-1,NX-NX/3+2:NX,1)=conj(flip(gi(2:end-1,2:NX/3,1),2));
vi(2:end-1,NX-NX/3+2:NX,1)=conj(flip(vi(2:end-1,2:NX/3,1),2));
ui(2:end-1,NX-NX/3+2:NX,1)=conj(flip(ui(2:end-1,2:NX/3,1),2));
wi(2:end-1,NX-NX/3+2:NX,1)=conj(flip(wi(2:end-1,2:NX/3,1),2));

v0=ifft2_cube(vi);
u0=ifft2_cube(ui);
w0=ifft2_cube(wi);
g0=ifft2_cube(gi);

u0=[zeros(1,NX,MZ);u0(2:end-1,:,:);zeros(1,NX,MZ)];
w0=[zeros(1,NX,MZ);w0(2:end-1,:,:);zeros(1,NX,MZ)];

Ei=Ener(u0,v0,w0);
kap=sc/(Ei^(1/2));
u0=kap*u0;
v0=kap*v0;
w0=kap*w0;
g0=kap*g0;

end

