function [FFkron,DYf,kkmf,llmf] = solvers_noise(FF,DY,Nk,Nm,kxband,kzband)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
global N MZ k l

DYf=kron(speye(Nk*Nm),DY);

FF=permute(FF,[1 2 4 3]);

%kxband = [1 2:3]; Nk = length(kxband);
%kzband = [1 2:MZ/3 MZ/3+3:2*MZ/3+1]; Nm = length(kzband);
%kzband = [1 2:MZ/3]; Nm = length([kzband MZ-flip(kzband(2:end))+2]);

FFkron=spalloc(2*N*(Nk)*(Nm),2*N*(Nk)*(Nm),4*N*N*(Nk)*(Nm));
kkmf=spalloc(N*(Nk)*(Nm),N*(Nk)*(Nm),N*(Nk)*(Nm));
llmf=kkmf;
lp([kzband 2*MZ/3-flip(kzband(2:end))+3]) = l([kzband end-flip(kzband(2:end))+2]);

for ii=kxband
    
    FFkronb=FF(:,:,ii,1);
    
    kpv=k(ii);
    lpv=lp(1);
    alp=kpv^2+lpv^2;
    
    lpb=zeros(N);
    if ii==1
        kpb=zeros(N);
    else
        kpb=1i*kpv*eye(N)/alp;
    end

    
%    for jk=[2:MZ/3  MZ-MZ/3+2:MZ]
    for jk=[kzband(2:end) 2*MZ/3-flip(kzband(2:end))+3]
        

        FFkronb=blkdiag(FFkronb,FF(:,:,ii,jk));
        
        lpv=lp(jk);
        alp=kpv^2+lpv^2;
        kpb=blkdiag(kpb,1i*kpv*eye(N)/alp);
        lpb=blkdiag(lpb,1i*lpv*eye(N)/alp);
        
    end
    
    ei=zeros(Nk,1);ei(ii)=1;
    
    IM=spdiags(ei,0,Nk,Nk);
    
    FFkron=FFkron+kron(IM,FFkronb);
    
    kkmf=kkmf+kron(IM,kpb);
    llmf=llmf+kron(IM,lpb);
    
end

