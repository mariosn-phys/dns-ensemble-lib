function [dFdz] = difZ_k(F,k,n)
%%% Real fields
global D1z D2z N NX MZ

if n==1
    Dif=D1z';
else
    Dif=D2z;
end

%dFdz=zeros(N+2,NX,MZ);

% for ix=1:NX
% dFdz(:,ix,:)=(Dif*permute(F(:,ix,:),[1 3 2])')';
% end

dFdz=F*Dif;

%dFdz=ipermute(Dif*permute(F,[2 1]),[2 1]);

end

