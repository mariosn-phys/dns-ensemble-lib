function [adv_v,adv_g] = advect_RNL_dtm(up,vp,wp,gp,u3d,v3d,w3d,g3d)

global N NX MZ

    %split 
    
%        uphat=fft(up,[],2);
%        vphat=fft(vp,[],2);
%        wphat=fft(wp,[],2);
%    
%        U_RNL=0*up;U_RNL(:,1,:)=uphat(:,1,:);U_RNL=ifft(U_RNL,[],2);
%        V_RNL=0*vp;V_RNL(:,1,:)=vphat(:,1,:);V_RNL=ifft(V_RNL,[],2);
%        W_RNL=0*wp;W_RNL(:,1,:)=wphat(:,1,:);W_RNL=ifft(W_RNL,[],2);
%    
%        u_RNL=uphat;u_RNL(:,1,:)=0*u_RNL(:,1,:);u_RNL=ifft(u_RNL,[],2);
%        v_RNL=vphat;v_RNL(:,1,:)=0*v_RNL(:,1,:);v_RNL=ifft(v_RNL,[],2);
%        w_RNL=wphat;w_RNL(:,1,:)=0*w_RNL(:,1,:);w_RNL=ifft(w_RNL,[],2);
        
%    replaced by

    U_RNL=repmat(mean(up,2),[1 NX 1]);
    V_RNL=repmat(mean(vp,2),[1 NX 1]);
    W_RNL=repmat(mean(wp,2),[1 NX 1]);

    u_RNL=up-U_RNL;
    v_RNL=vp-V_RNL;
    w_RNL=wp-W_RNL;


    dxdup=difX_F(u_RNL,1);
    dxdvp=difX_F(v_RNL,1);
    dxdwp=difX_F(w_RNL,1); 
    dxdUp=difX_F(U_RNL,1);
    dxdVp=difX_F(V_RNL,1);
    dxdWp=difX_F(W_RNL,1);
    
    
    dydup=difY_F(u_RNL,1);
    dydvp=difY_F(v_RNL,1);
    dydwp=difY_F(w_RNL,1); 
    dydUp=difY_F(U_RNL,1);
    dydVp=difY_F(V_RNL,1);
    dydWp=difY_F(W_RNL,1);
    
    
    dzdup=difZ_F(u_RNL,1);
    dzdvp=difZ_F(v_RNL,1);
    dzdwp=difZ_F(w_RNL,1);
    dzdUp=difZ_F(U_RNL,1);
    dzdVp=difZ_F(V_RNL,1);
    dzdWp=difZ_F(W_RNL,1);
    
    dxduu=repmat(mean(u_RNL.*dxdup,2),[1,NX,1]);
    dyduv=repmat(mean(v_RNL.*dydup,2),[1,NX,1]);
    dzduw=repmat(mean(w_RNL.*dzdup,2),[1,NX,1]);
    
    dxduv=repmat(mean(u_RNL.*dxdvp,2),[1,NX,1]);    
    dydvv=repmat(mean(v_RNL.*dydvp,2),[1,NX,1]);
    dzdvw=repmat(mean(w_RNL.*dzdvp,2),[1,NX,1]);
    
    dxduw=repmat(mean(u_RNL.*dxdwp,2),[1,NX,1]);    
    dydvw=repmat(mean(v_RNL.*dydwp,2),[1,NX,1]);
    dzdww=repmat(mean(w_RNL.*dzdwp,2),[1,NX,1]);
    
    %---------------------------------
    %v advection terms RNL add the means
    ad_u=(u_RNL+U_RNL).*dxdUp+U_RNL.*dxdup+dxduu+(v_RNL+V_RNL).*dydUp+V_RNL.*dydup+dyduv+(w_RNL+W_RNL).*dzdUp+W_RNL.*dzdup+dzduw+u3d;
    ad_v=(u_RNL+U_RNL).*dxdVp+U_RNL.*dxdvp+dxduv+(v_RNL+V_RNL).*dydVp+V_RNL.*dydvp+dydvv+(w_RNL+W_RNL).*dzdVp+W_RNL.*dzdvp+dzdvw+v3d;
    ad_w=(u_RNL+U_RNL).*dxdWp+U_RNL.*dxdwp+dxduw+(v_RNL+V_RNL).*dydWp+V_RNL.*dydwp+dydvw+(w_RNL+W_RNL).*dzdWp+W_RNL.*dzdwp+dzdww+w3d;
   
    
    %ad1=difX_F(ad_v,2)+difZ_F(ad_v,2);%+difY_F(ad_v,2)
    ad1=-difX_F(ad_w,1)+difZ_F(ad_u,1);
    ad2=-difY_F(difX_F(ad_u,1)+difZ_F(ad_w,1),1)+(difX_F(ad_v,2)+difZ_F(ad_v,2));

   
    
    adv_v=ad2;

    
    adv_g=ad1;


end

