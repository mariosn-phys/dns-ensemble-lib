function [ifft2_single ] = ifft2_cuben(spectral_field,Fn)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
if Fn==1
ifft2_single=real(permute(ifft2(permute(spectral_field,[2 3 1])),[3 1 2]));
else    
ifft2_single=real(permute(ifft2(permute(spectral_field,[2 3 1 4])),[3 1 2 4]));
end

end

