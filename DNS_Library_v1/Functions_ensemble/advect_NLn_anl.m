function [adv_v,adv_g] = advect_NLn_anl(up,vp,wp,gp,af,Fn)

global N NX MZ

    %split 
%         uphat=fft(up,[],2);
%         vphat=fft(vp,[],2);
%         wphat=fft(wp,[],2);
%     
%         U_RNL=0*up;U_RNL(:,1,:,:)=uphat(:,1,:,:);U_RNL=ifft(U_RNL,[],2);
%         V_RNL=0*vp;V_RNL(:,1,:,:)=vphat(:,1,:,:);V_RNL=ifft(V_RNL,[],2);
%         W_RNL=0*wp;W_RNL(:,1,:,:)=wphat(:,1,:,:);W_RNL=ifft(W_RNL,[],2);
%     
%         u_RNL=uphat;u_RNL(:,1,:,:)=0*u_RNL(:,1,:,:);u_RNL=ifft(u_RNL,[],2);
%         v_RNL=vphat;v_RNL(:,1,:,:)=0*v_RNL(:,1,:,:);v_RNL=ifft(v_RNL,[],2);
%         w_RNL=wphat;w_RNL(:,1,:,:)=0*w_RNL(:,1,:,:);w_RNL=ifft(w_RNL,[],2);
        
    U_RNL=repmat(mean(up,2),[1 NX 1 1]);
    V_RNL=repmat(mean(vp,2),[1 NX 1 1]);
    W_RNL=repmat(mean(wp,2),[1 NX 1 1]);

    u_RNL=up-U_RNL;
    v_RNL=vp-V_RNL;
    w_RNL=wp-W_RNL;
        

    dxdup=difX_Fn(u_RNL,1,Fn);
    dxdvp=difX_Fn(v_RNL,1,Fn);
    dxdwp=difX_Fn(w_RNL,1,Fn); 
    dxdUp=difX_Fn(U_RNL,1,Fn);
    dxdVp=difX_Fn(V_RNL,1,Fn);
    dxdWp=difX_Fn(W_RNL,1,Fn);
    
    
    dydup=difY_Fn(u_RNL,1,Fn);
    dydvp=difY_Fn(v_RNL,1,Fn);
    dydwp=difY_Fn(w_RNL,1,Fn); 
    dydUp=difY_Fn(U_RNL,1,Fn);
    dydVp=difY_Fn(V_RNL,1,Fn);
    dydWp=difY_Fn(W_RNL,1,Fn);
    
    
    dzdup=difZ_Fn(u_RNL,1,Fn);
    dzdvp=difZ_Fn(v_RNL,1,Fn);
    dzdwp=difZ_Fn(w_RNL,1,Fn);
    dzdUp=difZ_Fn(U_RNL,1,Fn);
    dzdVp=difZ_Fn(V_RNL,1,Fn);
    dzdWp=difZ_Fn(W_RNL,1,Fn);
    
    dxduu=repmat(mean(u_RNL.*dxdup,2),[1,NX,1,1]);
    dyduv=repmat(mean(v_RNL.*dydup,2),[1,NX,1,1]);
    dzduw=repmat(mean(w_RNL.*dzdup,2),[1,NX,1,1]);
    
    dxduv=repmat(mean(u_RNL.*dxdvp,2),[1,NX,1,1]);    
    dydvv=repmat(mean(v_RNL.*dydvp,2),[1,NX,1,1]);
    dzdvw=repmat(mean(w_RNL.*dzdvp,2),[1,NX,1,1]);
    
    dxduw=repmat(mean(u_RNL.*dxdwp,2),[1,NX,1,1]);    
    dydvw=repmat(mean(v_RNL.*dydwp,2),[1,NX,1,1]);
    dzdww=repmat(mean(w_RNL.*dzdwp,2),[1,NX,1,1]);
    
    %---------------------------------
    %v advection terms RNL and DNS_a add the means
%     if af == 0
%     ad_u=(u_RNL+U_RNL).*dxdUp+U_RNL.*dxdup+dxduu+(v_RNL+V_RNL).*dydUp+V_RNL.*dydup+dyduv+(w_RNL+W_RNL).*dzdUp+W_RNL.*dzdup+dzduw;
%     ad_v=(u_RNL+U_RNL).*dxdVp+U_RNL.*dxdvp+dxduv+(v_RNL+V_RNL).*dydVp+V_RNL.*dydvp+dydvv+(w_RNL+W_RNL).*dzdVp+W_RNL.*dzdvp+dzdvw;
%     ad_w=(u_RNL+U_RNL).*dxdWp+U_RNL.*dxdwp+dxduw+(v_RNL+V_RNL).*dydWp+V_RNL.*dydwp+dydvw+(w_RNL+W_RNL).*dzdWp+W_RNL.*dzdwp+dzdww;
%       
%     else
    ad_u=(u_RNL+U_RNL).*dxdUp+U_RNL.*dxdup+af.*(u_RNL.*dxdup-dxduu)+dxduu+(v_RNL+V_RNL).*dydUp+V_RNL.*dydup+af.*(v_RNL.*dydup-dyduv)+dyduv+(w_RNL+W_RNL).*dzdUp+W_RNL.*dzdup+af.*(w_RNL.*dzdup-dzduw)+dzduw;
    ad_v=(u_RNL+U_RNL).*dxdVp+U_RNL.*dxdvp+af.*(u_RNL.*dxdvp-dxduv)+dxduv+(v_RNL+V_RNL).*dydVp+V_RNL.*dydvp+af.*(v_RNL.*dydvp-dydvv)+dydvv+(w_RNL+W_RNL).*dzdVp+W_RNL.*dzdvp+af.*(w_RNL.*dzdvp-dzdvw)+dzdvw;
    ad_w=(u_RNL+U_RNL).*dxdWp+U_RNL.*dxdwp+af.*(u_RNL.*dxdwp-dxduw)+dxduw+(v_RNL+V_RNL).*dydWp+V_RNL.*dydwp+af.*(v_RNL.*dydwp-dydvw)+dydvw+(w_RNL+W_RNL).*dzdWp+W_RNL.*dzdwp+af.*(w_RNL.*dzdwp-dzdww)+dzdww;
%     end
    
    %ad1=difX_F(ad_v,2)+difZ_F(ad_v,2);%+difY_F(ad_v,2)
    ad1=-difX_Fn(ad_w,1,Fn)+difZ_Fn(ad_u,1,Fn);
    ad2=-difY_Fn(difX_Fn(ad_u,1,Fn)+difZ_Fn(ad_w,1,Fn),1,Fn)+(difX_Fn(ad_v,2,Fn)+difZ_Fn(ad_v,2,Fn));

   
    
    adv_v=ad2;

    
    adv_g=ad1;


end

