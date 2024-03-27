%Marios's Poiseuille linear optimals
clear all

R=450;
N=43;
m=0;
I=eye(N);
%K=1.25:1.25:5*1.25;%K=1;
K=2/1.75;
%K=linspace(1,3,19);

%dy=2/(N+1);
%y=linspace(-1+dy,1-dy,N)';

D1=cheb(N+1);flip(flip(D1,1),2);
[y,D4]=cheb4c(N+2);flip(flip(D4,1),2);
y=[-1;-y;1];

D2F=D1*D1;D2=D2F(2:end-1,2:end-1);

%V=1-y.^2;
V=y;
Vyy=D2F*V;

V=V(2:end-1);
Vyy=Vyy(2:end-1);

W2=INTweights(N+2,2);
W2=W2(2:end-1,2:end-1);
W1=sqrtm(W2);

% D2=zeros(N);
% for ia=2:N-1;
%     D2(ia,ia+1)=1;D2(ia,ia-1)=1;D2(ia,ia)=-2;
% end;
% D2(1,1)=-2;D2(1,2)=1;D2(N,N)=-2;D2(N,N-1)=1;
% D2=D2/(dy^2);


% D4=zeros(N);
% for ia=3:N-2;
%     D4(ia,ia+2)=1; D4(ia,ia+1)=-4;D4(ia,ia)=6;D4(ia,ia-2)=1; D4(ia,ia-1)=-4;
% end;
% 
% D4(1,1)=7;D4(1,2)=-4;D4(1,3)=1;D4(2,1)=-4;D4(2,2)=6;D4(2,3)=-4;D4(2,4)=1;
% D4(N-1,N-3)=1;D4(N-1,N-2)=-4;D4(N-1,N-1)=6;D4(N-1,N)=-4;D4(N,N)=7;D4(N,N-1)=-4;D4(N,N-2)=1;
% D4=D4/(dy^4);



TT=50;
%TT=40;

for ik=1:length(K);
    ik
k=K(ik);

alpha2=k^2+m^2;

DEL=D2-alpha2*I;
INVDEL=DEL\I;
DEL4=D4-2*alpha2*D2+alpha2^2*I;

%Symmetrize metric
DELSYM=1/2*(DEL+DEL');
INVDELSYM=DELSYM\I;


ij=sqrt(-1);

A=INVDEL*(-ij*k*diag(V)*DEL+ij*k*diag(Vyy)+(1/R)*DEL4);




%Integral Weights
METRIC=-DELSYM;METRIC=sqrtm(METRIC)*W2*sqrtm(METRIC);
MSQRT=sqrtm(METRIC);
%MSQRT=sqrtm(METRIC)*W1;
%METRIC=MSQRT*W2*MSQRT;
IMSQRT=MSQRT\I;
AM=MSQRT*A*IMSQRT;

end
%%

for it=1:length(TT);
    
topt=TT(it);
Phi=expm(AM*topt);
[u,s,v]=svd(Phi);
sd=diag(s);
v1=IMSQRT*v(:,1);
u1=IMSQRT*u(:,1);

figure(58)
subplot(121)
a=cartesianplotxy(v1,y,k,'contourf',20);
title([' k = ',num2str(k), ' t = ',num2str(topt)])

%figure(76)
subplot(122)
a=cartesianplotxy(u1,y,k,'contourf',20);
title([' E/E0 = ',num2str(sd(1)^2)])
hold off
drawnow

E(it,ik)=sd(1)^2;
end;

%%
%evolution of a specific optimal
tt=linspace(0,50,51);
dt=tt(2)-tt(1);
AE=expm(A*dt);
AET=expm(A'*dt);
vt=v1;
Ev(1)=real(vt'*METRIC*vt);
figure(90)
a=cartesianplotxy(vt,y,k,'contourf',20);
title([' t = ',num2str(tt(1)), ' E/E0 = ',num2str(Ev(1))])
drawnow
for ie=2:length(tt);
    vt=AE*vt;
%    v2=D2*(conj(vt).*vt+);
    Ev(ie)=real(vt'*METRIC*vt);
%    figure(74)
%    b=cartesianplotxy(v2,y,k,'contourf',20);
    figure(90)
    subplot(121)
    a=cartesianplotxy(vt,y,k,'contourf',20);
    
    title([' t = ',num2str(tt(1)), ' E/E0 = ',num2str(Ev(ie))])
    %drawnow
    %figure
    subplot(122)
    plot(tt(1:ie),Ev(1:ie),'.')
    drawnow
    
end;
%%

Evt(1)=real(vt'*METRIC*vt);
figure(91)
a=cartesianplotxy(vt,y,k,'contourf',20);
title([' t = ',num2str(tt(1)), ' E/E0 = ',num2str(Ev(1))])
drawnow
for ie=2:length(tt);
    vt=AET*vt;
%    v2=D2*(conj(vt).*vt+);
    Evt(ie)=real(vt'*METRIC*vt);
%    figure(74)
%    b=cartesianplotxy(v2,y,k,'contourf',20);
    figure(92)
    subplot(121)
    a=cartesianplotxy(vt,y,k,'contourf',20);
    
    title([' t = ',num2str(tt(1)), ' E/E0 = ',num2str(Ev(ie))])
    %drawnow
    %figure
    subplot(122)
    plot(-tt(1:ie),Evt(1:ie)/Evt(1),'.')
    drawnow
    
end;

%%
%end

% figure(10)
% %[V,LL]=eig(A);
% [V,L]=eig(A);
% [Ls,ind]=sort(diag(L));
% v1=V(:,ind(1));
% L=Ls;
% subplot(121)
% plot(i*L,'.','markersize',12)
% xlabel('c_r');ylabel('c_i');
% axis([0 1 -1 0])
% %subplot(122)
% %plot(V(:,end),y)
% %xlabel('v');
% subplot(122)
% a=cartesianplotxy(v1,y,k,'contourf',20);
% %ylabel('y')
% %title([' E/E0 = ',num2str(sd(1)^2)])
% hold off
% drawnow

