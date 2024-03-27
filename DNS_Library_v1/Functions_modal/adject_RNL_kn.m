function [adjv_v,adjv_g,adjU,adjW] = adject_RNL_kn(upa,vpa,wpa,gpa,umean,vmean,wmean,gmean,Lk)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

global a

    aLk = permute(a * Lk, [1 3 2]);
    Fn = length(Lk);

  %      upahat=fft(upa,[],2);
  %      vpahat=fft(vpa,[],2);
  %      wpahat=fft(wpa,[],2);
    
  %      Ua_RNL=zeros(N+2,NX,MZ);Ua_RNL(:,1,:)=upahat(:,1,:);Ua_RNL=ifft(Ua_RNL,[],2);
  %      Va_RNL=zeros(N+2,NX,MZ);Va_RNL(:,1,:)=vpahat(:,1,:);Va_RNL=ifft(Va_RNL,[],2);
  %      Wa_RNL=zeros(N+2,NX,MZ);Wa_RNL(:,1,:)=wpahat(:,1,:);Wa_RNL=ifft(Wa_RNL,[],2);

        ua_RNL=upa;%ua_RNL(:,1,:)=zeros(N+2,1,MZ);ua_RNL=ifft(ua_RNL,[],2);
        va_RNL=vpa;%va_RNL(:,1,:)=zeros(N+2,1,MZ);va_RNL=ifft(va_RNL,[],2);
        wa_RNL=wpa;%wa_RNL(:,1,:)=zeros(N+2,1,MZ);wa_RNL=ifft(wa_RNL,[],2);
        
 %       uphat=fft(up,[],2);
 %       vphat=fft(vp,[],2);
 %       wphat=fft(wp,[],2);
    
  %      U_RNL=zeros(N+2,NX,MZ);U_RNL(:,1,:)=uphat(:,1,:);U_RNL=ifft(U_RNL,[],2);
  %      V_RNL=zeros(N+2,NX,MZ);V_RNL(:,1,:)=vphat(:,1,:);V_RNL=ifft(V_RNL,[],2);
  %      W_RNL=zeros(N+2,NX,MZ);W_RNL(:,1,:)=wphat(:,1,:);W_RNL=ifft(W_RNL,[],2);
  
       U_RNL=umean;
       V_RNL=vmean;
       W_RNL=wmean;
    
  %      u_RNL=uphat;u_RNL(:,1,:)=zeros(N+2,1,MZ);u_RNL=ifft(u_RNL,[],2);
  %      v_RNL=vphat;v_RNL(:,1,:)=zeros(N+2,1,MZ);v_RNL=ifft(v_RNL,[],2);
  %      w_RNL=wphat;w_RNL(:,1,:)=zeros(N+2,1,MZ);w_RNL=ifft(w_RNL,[],2);
    
       U_RNL_dx = 1i * aLk .* U_RNL;
  
%    trapz(z,trapz(x,trapz(yE,(u0-Uback).*u0a+v0.*v0a+w0.*w0a)))
    
%     upavp=mean(mean(upa.*vp,3),2);
%     wpavp=mean(mean(wpa.*vp,3),2);
%     
%     UPa=mean(mean(upa,3),2);
%     WPa=mean(mean(wpa,3),2);
    
%     dxdupa=difX_F(upa,1);
%     dxdvpa=difX_F(vpa,1);
%     dxdwpa=difX_F(wpa,1);
    dxdU_RNL=difX_k(U_RNL,0,1);
    dxdV_RNL=difX_k(V_RNL,0,1);
    dxdW_RNL=difX_k(W_RNL,0,1);
%    dxdupa_RNL=difX_k(ua_RNL,Lk,1);
%    dxdvpa_RNL=difX_k(va_RNL,Lk,1);
%    dxdwpa_RNL=difX_k(wa_RNL,Lk,1);
    
%     dydupa=difY_F(upa,1);
%     dydvpa=difY_F(vpa,1);
%     dydwpa=difY_F(wpa,1);
    dydU_RNL=difY_k(U_RNL,0,1);
    dydV_RNL=difY_k(V_RNL,0,1);
    dydW_RNL=difY_k(W_RNL,0,1);
    dydupa_RNL=difY_kn(ua_RNL,Lk,1,Fn);
    dydvpa_RNL=difY_kn(va_RNL,Lk,1,Fn);
    dydwpa_RNL=difY_kn(wa_RNL,Lk,1,Fn);

%     dzdupa=difZ_F(upa,1);
%     dzdvpa=difZ_F(vpa,1);
%     dzdwpa=difZ_F(wpa,1);
    dzdU_RNL=difZ_k(U_RNL,0,1);
    dzdV_RNL=difZ_k(V_RNL,0,1);
    dzdW_RNL=difZ_k(W_RNL,0,1);
    dzdupa_RNL=difZ_kn(ua_RNL,Lk,1,Fn);
    dzdvpa_RNL=difZ_kn(va_RNL,Lk,1,Fn);
    dzdwpa_RNL=difZ_kn(wa_RNL,Lk,1,Fn);
    
    
    
    
    

    adj1_ua=U_RNL_dx.*ua_RNL+V_RNL.*dydupa_RNL+W_RNL.*dzdupa_RNL;
    adj1_va=U_RNL_dx.*va_RNL+V_RNL.*dydvpa_RNL+W_RNL.*dzdvpa_RNL;
    adj1_wa=U_RNL_dx.*wa_RNL+V_RNL.*dydwpa_RNL+W_RNL.*dzdwpa_RNL;
    
%     adj2_ua=U_RNL.*dxdupa_RNL+V_RNL.*dxdvpa_RNL+W_RNL.*dxdwpa_RNL;
%     adj2_va=U_RNL.*dydupa_RNL+V_RNL.*dydvpa_RNL+W_RNL.*dydwpa_RNL;
%     adj2_wa=U_RNL.*dzdupa_RNL+V_RNL.*dzdvpa_RNL+W_RNL.*dzdwpa_RNL;
    
    adj2_ua=-ua_RNL.*dxdU_RNL-va_RNL.*dxdV_RNL-wa_RNL.*dxdW_RNL;
    adj2_va=-ua_RNL.*dydU_RNL-va_RNL.*dydV_RNL-wa_RNL.*dydW_RNL;
    adj2_wa=-ua_RNL.*dzdU_RNL-va_RNL.*dzdV_RNL-wa_RNL.*dzdW_RNL;
    
    
%     adj3_ua=U_RNL.*dxdUpa+V_RNL.*dydUpa+W_RNL.*dzdUpa;
%     adj3_va=U_RNL.*dxdVpa+V_RNL.*dydVpa+W_RNL.*dzdVpa;
%     adj3_wa=U_RNL.*dxdWpa+V_RNL.*dydWpa+W_RNL.*dzdWpa;
%     
%     adj4_ua=U_RNL.*dxdUpa+V_RNL.*dxdVpa+W_RNL.*dxdWpa;
%     adj4_va=U_RNL.*dydUpa+V_RNL.*dydVpa+W_RNL.*dydWpa;
%     adj4_wa=U_RNL.*dzdUpa+V_RNL.*dzdVpa+W_RNL.*dzdWpa;
    
%     adj5_ua=repmat(mean(u_RNL.*dxdupa+v_RNL.*dydupa+w_RNL.*dzdupa,2),[1,NX,1]);
%     adj5_va=repmat(mean(u_RNL.*dxdvpa+v_RNL.*dydvpa+w_RNL.*dzdvpa,2),[1,NX,1]);
%     adj5_wa=repmat(mean(u_RNL.*dxdwpa+v_RNL.*dydwpa+w_RNL.*dzdwpa,2),[1,NX,1]);
%     
%     adj6_ua=repmat(mean(u_RNL.*dxdupa+v_RNL.*dxdvpa+w_RNL.*dxdwpa,2),[1,NX,1]);
%     adj6_va=repmat(mean(u_RNL.*dydupa+v_RNL.*dydvpa+w_RNL.*dydwpa,2),[1,NX,1]);
%     adj6_wa=repmat(mean(u_RNL.*dzdupa+v_RNL.*dzdvpa+w_RNL.*dzdwpa,2),[1,NX,1]);
    
%     adj3_ua=u_RNL.*dxdUpa+v_RNL.*dydUpa+w_RNL.*dzdUpa;
%     adj3_va=u_RNL.*dxdVpa+v_RNL.*dydVpa+w_RNL.*dzdVpa;
%     adj3_wa=u_RNL.*dxdWpa+v_RNL.*dydWpa+w_RNL.*dzdWpa;
%     
%     adj4_ua=u_RNL.*dxdUpa+v_RNL.*dxdVpa+w_RNL.*dxdWpa;
%     adj4_va=u_RNL.*dydUpa+v_RNL.*dydVpa+w_RNL.*dydWpa;
%     adj4_wa=u_RNL.*dzdUpa+v_RNL.*dzdVpa+w_RNL.*dzdWpa;
    
%     adj_u=adj1_ua+adj2_ua+adj3_ua+adj4_ua+adj5_ua+adj6_ua+adj7_ua+adj8_ua;
%     adj_v=adj1_va+adj2_va+adj3_va+adj4_va+adj5_va+adj6_va+adj7_va+adj8_va;
%     adj_w=adj1_wa+adj2_wa+adj3_wa+adj4_wa+adj5_wa+adj6_wa+adj7_wa+adj8_wa;


    adj_u=adj1_ua+adj2_ua;%+adj3_ua+adj4_ua+adj5_ua+adj6_ua+adj7_ua+adj8_ua;
    adj_v=adj1_va+adj2_va;%+adj3_va+adj4_va+adj5_va+adj6_va+adj7_va+adj8_va;
    adj_w=adj1_wa+adj2_wa;%+adj3_wa+adj4_wa+adj5_wa+adj6_wa+adj7_wa+adj8_wa;


%    adj1_ua=up.*dxdupa+vp.*dydupa+wp.*dzdupa;
% %   adj1_va=up.*dxdvpa+vp.*dydvpa+wp.*dzdvpa;
%    adj1_wa=up.*dxdwpa+vp.*dydwpa+wp.*dzdwpa;
%     
%     adj2_ua=up.*dxdupa+vp.*dxdvpa+wp.*dxdwpa;
% %    adj2_va=up.*dydupa+vp.*dydvpa+wp.*dydwpa;
%     adj2_wa=up.*dzdupa+vp.*dzdvpa+wp.*dzdwpa;
    
    adj1=-1i*aLk.*adj_w+difZ_kn(adj_u,Lk,1,Fn);
    adj2=-difY_kn(1i*aLk.*adj_u+difZ_kn(adj_w,Lk,1,Fn),Lk,1,Fn)+(-aLk.^2.*adj_v+difZ_kn(adj_v,Lk,2,Fn));
  
%    adjU=mean(mean(adj1_ua+adj2_ua,3),2);
%    adjW=mean(mean(adj1_wa+adj2_wa,3),2);
    
    adjv_g=adj1;
    adjv_v=adj2;

end

