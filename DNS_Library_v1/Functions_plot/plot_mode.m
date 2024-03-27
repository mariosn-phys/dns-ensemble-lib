function plot_mode(kfig,uL,vL,wL,gL,Ln,Lk,um,vm,wm)

[un,vn,wn,gn]=pick_vec(uL,vL,wL,gL,Ln,Lk);

figure(kfig+1);clf
subplot(221)

    Uyz_t=mean(um,2);
%    Uyz_t=squeeze(mean(u0,2)-repmat(UP1,[1,1,MZ]));
    Vyz_t=mean(vm,2);
    Wyz_t=mean(wm,2);
       
   
    plot_contour(Uyz_t,Vyz_t,Wyz_t,'z','y','yz')
    
    subplot(222)
    
    Uyz_t=un(:,1,:);
%    Uyz_t=squeeze(mean(u0,2)-repmat(UP1,[1,1,MZ]));
    Vyz_t=vn(:,1,:);
    Wyz_t=wn(:,1,:);
    
    
    plot_contour(Uyz_t,Vyz_t,Wyz_t,'z','y','yz')
    
    subplot(223)
    
    Uyz_t=un(:,:,1);
%    Uyz_t=squeeze(mean(u0,2)-repmat(UP1,[1,1,MZ]));
    Vyz_t=vn(:,:,1);
    Wyz_t=wn(:,:,1);
    
    
    plot_contour(Uyz_t,Vyz_t,Wyz_t,'x','y','yx')
        subplot(224)
    
    Uyz_t=un(11,:,:);
%    Uyz_t=squeeze(mean(u0,2)-repmat(UP1,[1,1,MZ]));
    Vyz_t=vn(11,:,:);
    Wyz_t=wn(11,:,:);
    
    
    plot_contour(Uyz_t,Vyz_t,Wyz_t,'x','z','xz')