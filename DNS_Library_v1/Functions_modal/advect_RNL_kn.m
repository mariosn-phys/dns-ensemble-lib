function [adv_v,adv_g] = advect_RNL_kn(up,vp,wp,gp,umean,vmean,wmean,gmean,Lk)

    global a DYF D1z N MZ

    aLk = permute(a * Lk, [1 3 2]);
    Fn = length(Lk);
    
    %split 
%         uphat=fft(up,[],2);
%         vphat=fft(vp,[],2);
%         wphat=fft(wp,[],2);
    
        U_RNL=umean;
        V_RNL=vmean;
        W_RNL=wmean;
    
        u_RNL=up;
        v_RNL=vp;
        w_RNL=wp;


%     dxdup=difX_k(u_RNL,Lk,1);
%     dxdvp=difX_k(v_RNL,Lk,1);
%     dxdwp=difX_k(w_RNL,Lk,1); 
%     dxdUp=difX_k(U_RNL,0,1);
%     dxdVp=difX_k(V_RNL,0,1);
%     dxdWp=difX_k(W_RNL,0,1);

%   Ensemble dx 

    U_RNL_dx=1i*aLk.*U_RNL;
    
%   U_RNL_dx.*up;U_RNL_dx.*vp;U_RNL_dx.*wp;

    dydup=difY_kn(u_RNL,Lk,1,Fn);
    dydvp=difY_kn(v_RNL,Lk,1,Fn);
    dydwp=difY_kn(w_RNL,Lk,1,Fn); 
    dydUp=difY_k(U_RNL,0,1);
    dydVp=difY_k(V_RNL,0,1);
    dydWp=difY_k(W_RNL,0,1);
    
    
    dzdup=difZ_kn(u_RNL,Lk,1,Fn);
    dzdvp=difZ_kn(v_RNL,Lk,1,Fn);
    dzdwp=difZ_kn(w_RNL,Lk,1,Fn);
    dzdUp=difZ_k(U_RNL,0,1);
    dzdVp=difZ_k(V_RNL,0,1);
    dzdWp=difZ_k(W_RNL,0,1);

    
%     dxduu=repmat(mean(u_RNL.*dxdup,2),[1,NX,1]);
%     dyduv=repmat(mean(v_RNL.*dydup,2),[1,NX,1]);
%     dzduw=repmat(mean(w_RNL.*dzdup,2),[1,NX,1]);
%     
%     dxduv=repmat(mean(u_RNL.*dxdvp,2),[1,NX,1]);    
%     dydvv=repmat(mean(v_RNL.*dydvp,2),[1,NX,1]);
%     dzdvw=repmat(mean(w_RNL.*dzdvp,2),[1,NX,1]);
%     
%     dxduw=repmat(mean(u_RNL.*dxdwp,2),[1,NX,1]);    
%     dydvw=repmat(mean(v_RNL.*dydwp,2),[1,NX,1]);
%     dzdww=repmat(mean(w_RNL.*dzdwp,2),[1,NX,1]);
    
    %---------------------------------
    %v advection terms RNL add the means
%     ad_u=(u_RNL+U_RNL).*dxdUp+U_RNL.*dxdup+dxduu+(v_RNL+V_RNL).*dydUp+V_RNL.*dydup+dyduv+(w_RNL+W_RNL).*dzdUp+W_RNL.*dzdup+dzduw;
%     ad_v=(u_RNL+U_RNL).*dxdVp+U_RNL.*dxdvp+dxduv+(v_RNL+V_RNL).*dydVp+V_RNL.*dydvp+dydvv+(w_RNL+W_RNL).*dzdVp+W_RNL.*dzdvp+dzdvw;
%     ad_w=(u_RNL+U_RNL).*dxdWp+U_RNL.*dxdwp+dxduw+(v_RNL+V_RNL).*dydWp+V_RNL.*dydwp+dydvw+(w_RNL+W_RNL).*dzdWp+W_RNL.*dzdwp+dzdww;
   
%     ad_u=U_RNL.*dxdup+(v_RNL+V_RNL).*dydUp+V_RNL.*dydup+(w_RNL+W_RNL).*dzdUp+W_RNL.*dzdup;
%     ad_v=U_RNL.*dxdvp+(v_RNL+V_RNL).*dydVp+V_RNL.*dydvp+(w_RNL+W_RNL).*dzdVp+W_RNL.*dzdvp;
%     ad_w=U_RNL.*dxdwp+(v_RNL+V_RNL).*dydWp+V_RNL.*dydwp+(w_RNL+W_RNL).*dzdWp+W_RNL.*dzdwp;
    
    ad_u=U_RNL_dx.*up+v_RNL.*dydUp+V_RNL.*dydup+w_RNL.*dzdUp+W_RNL.*dzdup;
    ad_v=U_RNL_dx.*vp+v_RNL.*dydVp+V_RNL.*dydvp+w_RNL.*dzdVp+W_RNL.*dzdvp;
    ad_w=U_RNL_dx.*wp+v_RNL.*dydWp+V_RNL.*dydwp+w_RNL.*dzdWp+W_RNL.*dzdwp;
    
%     %ad1=difX_F(ad_v,2)+difZ_F(ad_v,2);%+difY_F(ad_v,2)
%     ad1=-difX_k(ad_w,Lk,1)+difZ_k(ad_u,Lk,1);
%     ad2=-difY_k(difX_k(ad_u,Lk,1)+difZ_k(ad_w,Lk,1),Lk,1)+(difX_k(ad_v,Lk,2)+difZ_k(ad_v,Lk,2));
    
    %unwrap difX
    %ad1=difX_F(ad_v,2)+difZ_F(ad_v,2);%+difY_F(ad_v,2)
    ad1=-1i*aLk.*ad_w+difZ_kn(ad_u,Lk,1,Fn);
    ad2=-difY_kn(1i*aLk.*ad_u+difZ_kn(ad_w,Lk,1,Fn),Lk,1,Fn)+(-aLk.^2.*ad_v+difZ_kn(ad_v,Lk,2,Fn));
            
    adv_v=ad2;

    adv_g=ad1;

end

