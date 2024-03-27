function [dFdx] = difX_k(F,k,n)

global a

Dif=1i*a*k;
dFdx=(Dif^n)*F;

% for iz=1:MZ
% dFdx(:,:,iz)=(Dif*squeeze(F(:,:,iz))')';
% end

%DIFX Differentiation in X
%   Detailed explanation goes here


end

