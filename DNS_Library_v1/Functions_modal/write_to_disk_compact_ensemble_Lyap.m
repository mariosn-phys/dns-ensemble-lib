function write_to_disk_compact_ensemble_Lyap(vi,gi,UP,WP,runtime,kLyap,NLyap,filename);
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

global NX MZ

vhat=vi(:,1:NX/3+1,[1:MZ/3+1 MZ-MZ/3+1:MZ],:,:,:);
ghat=gi(:,1:NX/3+1,[1:MZ/3+1 MZ-MZ/3+1:MZ],:,:,:);

save(filename,'vhat','ghat','UP','WP','NLyap','kLyap','runtime')

end

