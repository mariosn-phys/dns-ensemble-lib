function [dFdx] = difX_Fn(F,n,Fn)

% Ensemble state dx derivatives

global D1x D2x N NX MZ

if n==1
    Dif=D1x;
elseif n==2
    Dif=D2x;
end

%dFdx=zeros(N+2,NX,MZ);

% for iz=1:MZ
% dFdx(:,:,iz)=(Dif*F(:,:,iz)')';
% end
%DIFX Differentiation in X
%   Detailed explanation goes here

dFdx=permute(reshape(Dif*reshape(permute(F,[2 1 3 4]),[NX,(N+2)*MZ*Fn]),[NX,N+2,MZ,Fn]),[2 1 3 4]);

end

