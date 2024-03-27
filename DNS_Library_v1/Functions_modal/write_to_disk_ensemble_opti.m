function write_to_disk_ensemble_opti(u0,v0,w0,g0,runtime,kLyap,NLyap,gr,filename)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

    UP=mean(mean(u0(:,:,:,1),3),2);
    WP=mean(mean(w0(:,:,:,1),3),2);

    for ik=1:length(kLyap)
    for ij=1:NLyap
    vi(:,:,:,ij,ik)=fft2_cube(v0(:,:,:,ij,ik));
    gi(:,:,:,ij,ik)=fft2_cube(g0(:,:,:,ij,ik));
    end
    end

    write_to_disk_compact_ensemble_opti(vi,gi,UP,WP,gr,runtime,kLyap,NLyap,filename);

end

