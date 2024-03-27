function [dFdy] = difY_k(F,k,n)

global DYF D2F N MZ

if n==1
    Dif=DYF;
else
    Dif=D2F;
end

%dFdy=zeros(N+2,MZ);
 
% for iz=1:MZ
% dFdy(:,:,iz)=Dif*F(:,:,iz);
% end

dFdy=Dif*F;


end

