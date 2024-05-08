function NoiseStructure2_cheb(R,N,NX,MZ,a,b)

%Noise Structures / Diffusion Modes

L=2;

% increments of 12

A=2*pi/a;dx=A/NX;x=-dx+dx*(1:NX)';xE=[x;A]; 
B=2*pi/b;dz=B/MZ;z=-dz+dz*(1:MZ)';zE=[z;B];

%m=0;
I=eye(N);
I2=eye(2*N);
%K=1.25:1.25:10*1.25;%K=1;
%K=0.25;%K=0.46;
M=[0:(NX/2-1) 0 (1-NX/2):(-1)];
    K=2*pi*M(1:NX/3+1)/A; 
P=[0:(MZ/2-1) 0 (1-MZ/2):(-1)];
    M=2*pi*P([1:MZ/3+1 MZ-MZ/3+1:MZ])/B;  

%K=0.5;
%M=10/3;%M=1.9;M=1.66;
%K=linspace(0.1,0.4,4);
%M=linspace(1,3,20);

D1F=cheb(N+1);flip(flip(D1F,1),2);
[y,D4]=cheb4c(N+2);D4=flip(flip(D4,1),2);y=[-1;-y;1];

D2F=D1F*D1F;D2=D2F(2:end-1,2:end-1);
D1=D1F(2:end-1,2:end-1);

W2F=INTweights(N+2,2);
W2=W2F(2:end-1,2:end-1);
W1=sqrtm(W2);



% V=1-y.^2;
% V=y;
% Vy=D1F*V;
% Vyy=D2F*V;



TT=0.1;
for ik=1:length(K);
    for il=1:length(M);
%    ik
%    il
k=K(ik);
m=M(il);
alpha2=k^2+m^2;

if ik==1 && il==1
    FF(:,:,il,ik)=zeros(2*N);

else    
    
DEL=D2-alpha2*I;
INVDEL=DEL\I;
DEL4=D4-2*alpha2*D2+alpha2^2*I;

A=[INVDEL*DEL4 zeros(N);zeros(N) DEL];

%Symmetrize metric
%DELSYM=1/2*(DEL+DEL');
%INVDELSYM=DELSYM\I;

ij=sqrt(-1);

%Integral Weights
METRIC=-DEL;METRIC=sqrtm(METRIC)'*W2*sqrtm(METRIC);
MSQRT=sqrtm(METRIC);
%MIC=[METRIC zeros(N);zeros(N) W2F]/alpha2;
MIC=[METRIC zeros(N);zeros(N) W2]/alpha2;

MCSQR=sqrtm(MIC);
IMCSQR=MCSQR\I2;

%MIC=[METRIC zeros(N);zeros(N) I];
%MSQRT=sqrtm(MIC);
AM = MCSQR*A*IMCSQR;

MTC(:,:,ik,il)=AM;


    [U,S,V]=svd(AM);

    %FF(:,:,il,ik)=V*diag(1./sqrt(diag(S)))/(alpha2^(1/6));
    %FF(:,:,il,ik)=V*diag(1./sqrt(diag(S)))/(alpha2^(2/5));
       
    FF(:,:,il,ik)=IMCSQR*U;%'noise_structure_dfreedom_eqall_434848_DNS.mat','FF'
    %FF(:,:,il,ik)=IMCSQR*U*diag(1./sqrt(diag(S)));%'noise_structure_dfreedom_eqall_434848_DNS.mat','FF'
    %diag(FF(:,:,il,ik)'*MIC*FF(:,:,il,ik))
end

    end
end

%Test Energy 

save(['noise_structure_diffusion_eqall_',num2str(N),num2str(NX),num2str(MZ),'.mat'],'FF','MTC')

%save('noise_structure_dfreedom_lambdam1_km8_5_434848_DNS_correct.mat','FF')

