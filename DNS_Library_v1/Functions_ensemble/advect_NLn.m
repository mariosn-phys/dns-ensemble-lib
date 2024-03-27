function [adv_v,adv_g] = advect_NLn(up,vp,wp,gp,Fn)


    dxdup=difX_Fn(up,1,Fn);
    dxdvp=difX_Fn(vp,1,Fn);
    dxdwp=difX_Fn(wp,1,Fn);
    %d2xdup=difX_Fn(up,2);
    %d2xdvp=difX_Fn(vp,2);
    %d2xdwp=difX_Fn(wp,2);
    
    
    dydup=difY_Fn(up,1,Fn);
    dydvp=difY_Fn(vp,1,Fn);
    dydwp=difY_Fn(wp,1,Fn);  
    %d2ydup=difY_Fn(up,2);
    %d2ydvp=difY_Fn(vp,2);
    %d2ydwp=difY_Fn(wp,2);
    
    
    dzdup=difZ_Fn(up,1,Fn);
    dzdvp=difZ_Fn(vp,1,Fn);
    dzdwp=difZ_Fn(wp,1,Fn);
    %d2zdup=difZ_Fn(up,2);
    %d2zdvp=difZ_Fn(vp,2);
    %d2zdwp=difZ_Fn(wp,2);
    
    %dxdgp=difX_F(gp,1);
    %dydgp=difY_F(gp,1);
    %dzdgp=difZ_F(gp,1);
    

    %D2gp=difX_F(gp,2)+difY_F(gp,2)+difZ_F(gp,2);

    %D2v=difX_F(vp,2)+difY_F(vp,2)+difZ_F(vp,2);

    %---------------------------------
    %v advection terms 
    ad_u=up.*dxdup+vp.*dydup+wp.*dzdup;
    ad_v=up.*dxdvp+vp.*dydvp+wp.*dzdvp;
    ad_w=up.*dxdwp+vp.*dydwp+wp.*dzdwp;
    %ad_g=up.*dxdgp+vp.*dydgp+wp.*dzdgp;
    
    %dzad_u=dzdup.*dxdup+dzdvp.*dydup+dzdwp.*dzdup;
    %dxad_w=dxdup.*dxdwp+dxdvp.*dydwp+dxdwp.*dzdwp;
    
    
    %ad1=difX_F(ad_v,2)+difZ_F(ad_v,2);%+difY_F(ad_v,2)
    ad1=-difX_Fn(ad_w,1,Fn)+difZ_Fn(ad_u,1,Fn);
    ad2=-difY_Fn(difX_Fn(ad_u,1,Fn)+difZ_Fn(ad_w,1,Fn),1,Fn)+(difX_Fn(ad_v,2,Fn)+difZ_Fn(ad_v,2,Fn));
   
    
    adv_v=ad2;

    
    adv_g=ad1;


end

