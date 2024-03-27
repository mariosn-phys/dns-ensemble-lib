function output_ensemble(vi,gi,UP,WP,T,folder,stname,Fnn)
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here

for nn=1:length(Fnn)
    field_path=[folder,'N',num2str(Fnn(nn),'%03.f'),'/'];
    write_to_disk_compact(vi(:,:,:,nn),gi(:,:,:,nn),UP(:,nn),WP(:,nn),T,[field_path,stname]);
end

end