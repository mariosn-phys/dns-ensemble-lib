function [dFdy] = difY_kn(F,k,n,Fn)

global DYF D2F N MZ

if n==1
    Dif=DYF;
else
    Dif=D2F;
end

%dFdy=zeros(N+2,MZ,Fn);
 
% for iz=1:MZ
% dFdy(:,:,iz)=Dif*F(:,:,iz);
% end

dFdy=reshape(Dif*reshape(F,[N+2, MZ*Fn]),[N+2, MZ, Fn]);


end

