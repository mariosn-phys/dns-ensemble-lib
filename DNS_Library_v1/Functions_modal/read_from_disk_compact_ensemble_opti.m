function [vi,gi,UP,WP,gr,NLyap,KL,runtime]=read_from_disk_compact_ensemble_opti(filename)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

global NX MZ N

load(filename,'vhat','ghat','UP','WP','gr','NLyap','kLyap','runtime')


Q=size(vhat);


KL=length(kLyap);

vi=zeros(N+2,NX,MZ,Q(4),KL);
gi=zeros(N+2,NX,MZ,Q(4),KL);

NXfill=min(Q(2),NX/3+1);
MZfill=min((Q(3)+1)/2,MZ/3+1);    


vi(:,1:NXfill,[1:MZfill MZ-MZfill+2:MZ],:,:)=vhat(:,1:NXfill,[1:MZfill Q(3)-MZfill+2:Q(3)],:,:);
gi(:,1:NXfill,[1:MZfill MZ-MZfill+2:MZ],:,:)=ghat(:,1:NXfill,[1:MZfill Q(3)-MZfill+2:Q(3)],:,:);
    
%    for jj=2:min(Q(2),NX/3+1)
        
%vi(:,1:NX/3+1,[1:MZ/3+1 MZ-MZ/3+1:MZ])=vhat;
%gi(:,1:NX/3+1,[1:MZ/3+1 MZ-MZ/3+1:MZ])=ghat;

%    end
%end

% Fill fields

end

