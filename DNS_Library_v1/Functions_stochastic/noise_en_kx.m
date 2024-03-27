function [u0,v0,w0,g0]=noise_en_kx(FF,F,sc,kx,kz,Fn)

% noise from complete optimal base
% total energy scaling

global N NX MZ k l DYF

vi=zeros(N+2,Fn,NX,MZ);
gi=zeros(N+2,Fn,NX,MZ);
ui=zeros(N+2,Fn,NX,MZ);
wi=zeros(N+2,Fn,NX,MZ);

%jjj=[1:MZ/3+1 MZ-MZ/3+1:MZ];

jjj=[kz MZ-flip(kz(2:end))+2];

for kk=kx%NX/3+1
    
    kp=k(kk);
    
    for ij=1:length(jjj)
        
        jj=jjj(ij);
        
        lp=l(jj);
        
        if kk==1 && ij==1
        xi=zeros(F,Fn);    
        else
        xi=(randn(F,Fn)+1i*randn(F,Fn))/sqrt(2);
        %xi=xi./abs(xi);
        end
        nn=FF(:,1:F,ij,kk)*xi;
        
        gihat=[zeros(1,Fn);nn(N+1:2*N,:);zeros(1,Fn)];
        vihat=[zeros(1,Fn);nn(1:N,:);zeros(1,Fn)];
        
        
        gi(:,:,kk,jj)=gihat;
        vi(:,:,kk,jj)=vihat;
                
        if (kp==0) && (lp==0)
        else
            dydvihat=DYF*vihat;
            ui(:,:,kk,jj)=(-1i*lp*gihat+1i*kp*dydvihat)/(kp^2+lp^2);
            wi(:,:,kk,jj)=(1i*kp*gihat+1i*lp*dydvihat)/(kp^2+lp^2);
        end
        
    end
    
end

gi=permute(gi,[1 3 4 2]);
vi=permute(vi,[1 3 4 2]);
ui=permute(ui,[1 3 4 2]);
wi=permute(wi,[1 3 4 2]);


for kk=2:kx(end)+1
    
    gi(:,NX-kk+2,1,:)=conj(gi(:,kk,1,:));
    vi(:,NX-kk+2,1,:)=conj(vi(:,kk,1,:));
    ui(:,NX-kk+2,1,:)=conj(ui(:,kk,1,:));
    wi(:,NX-kk+2,1,:)=conj(wi(:,kk,1,:));
    
    jj=[2:MZ/3+1 MZ-MZ/3+1:MZ];
    gi(:,NX-kk+2,MZ-jj+2,:)=conj(gi(:,kk,jj,:));
    vi(:,NX-kk+2,MZ-jj+2,:)=conj(vi(:,kk,jj,:));
    ui(:,NX-kk+2,MZ-jj+2,:)=conj(ui(:,kk,jj,:));
    wi(:,NX-kk+2,MZ-jj+2,:)=conj(wi(:,kk,jj,:));
        
end

v0=ifft2_cuben(vi,Fn);
u0=ifft2_cuben(ui,Fn);
w0=ifft2_cuben(wi,Fn);
g0=ifft2_cuben(gi,Fn);

u0=[zeros(1,NX,MZ,Fn);u0(2:end-1,:,:,:);zeros(1,NX,MZ,Fn)];
w0=[zeros(1,NX,MZ,Fn);w0(2:end-1,:,:,:);zeros(1,NX,MZ,Fn)];

% Ei=EnerFn(u0,v0,w0);
% kap=sc/(Ei^(1/2));
% u0=kap*u0;
% v0=kap*v0;
% w0=kap*w0;
% g0=kap*g0;

end
