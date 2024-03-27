function [dFdz] = difZ_kn(F,k,n,Fn)
%%% Real fields
global D1z D2z N NX MZ

if n==1
    Dif=D1z;
else
    Dif=D2z;
end

%dFdz=zeros(N+2,NX,MZ);

% for ix=1:NX
% dFdz(:,ix,:)=(Dif*permute(F(:,ix,:),[1 3 2])')';
% end

%dFdz=F*Dif;

%dFdz=permute(reshape(Dif*reshape(permute(F,[3 1 2]),[MZ,(N+2)*NX]),[MZ,N+2,NX]),[2 3 1]);

dFdz=ipermute(reshape(Dif*reshape(permute(F,[2 1 3]),[MZ,(N+2)*Fn]),[MZ,N+2,Fn]),[2 1 3]);

end

