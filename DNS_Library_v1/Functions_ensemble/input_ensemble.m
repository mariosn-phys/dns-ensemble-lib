function [ue,ve,we,ge] = input_ensemble(folder,stname,Fnn)

% Compact ensemble state load function 

global b1 b2

for nn=Fnn
    start_file=[folder,'N',num2str(nn,'%03.f'),'/',stname];
    %    Load init and transform to physical space
    [vi,gi,UP1,WP1]=read_from_disk_compact(start_file);
    [u0,v0,w0,g0] = make_uw(gather(vi),gather(gi),gather(UP1),gather(WP1),b1,b2);
    if nn==Fnn(1)
    ue=u0;ve=v0;we=w0;ge=g0;
    else
    ue=cat(4,ue,u0);ve=cat(4,ve,v0);we=cat(4,we,w0);ge=cat(4,ge,g0); 
    end

end
