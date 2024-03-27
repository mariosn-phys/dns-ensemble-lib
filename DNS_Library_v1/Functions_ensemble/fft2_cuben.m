function [fft2_single ] = fft2_cuben(physical_field,Fn)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

if Fn==1
fft2_single=permute(fft2(permute(physical_field,[2 3 1])),[3 1 2]);
else
fft2_single=permute(fft2(permute(physical_field,[2 3 1 4])),[3 1 2 4]);
end

end

