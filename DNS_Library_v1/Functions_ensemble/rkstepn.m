function [ui,vi,wi,gi,UP1,WP1] = rkstepn(u0,v0,w0,g0,advg1,advv1,advU1,advW1,step,h,Fn)

global S_mf S_mp Sol_m b1 b2 DYF dt gamma


% u0 the base flow
% u1 the advection flow
    b1n=repmat(b1,[1 Fn]);
    b2n=repmat(b2,[1 Fn]);
    
    z0n=zeros(1,Fn);
    
    if Fn==1
    UP=mean(mean(u0,3),2);
    WP=mean(mean(w0,3),2);  
    else
    UP=permute(mean(mean(u0,3),2),[1 4 2 3]);
    WP=permute(mean(mean(w0,3),2),[1 4 2 3]);
    end

    if Fn==1
    [gi,vi,ui,wi] = solv_vg_f_kron_zx_RK3(g0,v0,advg1,advv1,step,h);
    else
    [gi,vi,ui,wi] = solv_vgn_f_kron_zx_RK3(g0,v0,advg1,advv1,step,h,Fn);
    end

    dUP=S_mp(:,:,step)*UP-h*dt*advU1;
    dWP=S_mp(:,:,step)*WP-h*dt*advW1;
    
    UP1=[b1n;Sol_m(:,:,step)*(dUP(2:end-1,:)-S_mf(2:end-1,1,step)*b1n-S_mf(2:end-1,end,step)*b2n);b2n];  
    WP1=[z0n;Sol_m(:,:,step)*dWP(2:end-1,:);z0n];
        
end

