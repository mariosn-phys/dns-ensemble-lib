function [uL,vL,wL,gL] = Adjoi_tt_module_kron_RK3_v2_k(uLn,vLn,wLn,gLn,umean,vmean,wmean,gmean,NT,h,Lkn)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

global kLyap NLyap N NX MZ dt


NN = sum(NLyap); NN2= NN/2;

uL = 0 * uLn;
vL = uL;
wL = uL;
gL = uL;

uLt=uLn(:,:,1:2:end);
vLt=vLn(:,:,1:2:end);
wLt=wLn(:,:,1:2:end);
gLt=gLn(:,:,1:2:end);

for tt=2:NT

[adv_v1,adv_g1] = adject_RNL_kn(uLt,vLt,wLt,gLt,umean,vmean,wmean,gmean,Lkn);


        [gi,vi,ui,wi] = solv_adj_vg_RK_kron_v2_k(gLt,vLt,adv_g1,adv_v1,1,1/2*h);
                        
%         %for qq=1:N+2
             vpw=ifft(vi,[],2);
             upw=ifft(ui,[],2);
             wpw=ifft(wi,[],2);
             gpw=ifft(gi,[],2);
%         %end
        
        upw=[zeros(1,MZ,NN2);upw(2:end-1,:,:);zeros(1,MZ,NN2)];
        wpw=[zeros(1,MZ,NN2);wpw(2:end-1,:,:);zeros(1,MZ,NN2)];
   
[adv_v2,adv_g2] = adject_RNL_kn(upw,vpw,wpw,gpw,umean,vmean,wmean,gmean,Lkn);

       
        [gi,vi,ui,wi] = solv_adj_vg_RK_kron_v2_k(gLt,vLt,-adv_g1+2*adv_g2,-adv_v1+2*adv_v2,2,1*h);  
        
%         %for qq=1:N+2
             vpw=ifft(vi,[],2);
             upw=ifft(ui,[],2);
             wpw=ifft(wi,[],2);
             gpw=ifft(gi,[],2);
%         %end
        
        upw=[zeros(1,MZ,NN2);upw(2:end-1,:,:);zeros(1,MZ,NN2)];
        wpw=[zeros(1,MZ,NN2);wpw(2:end-1,:,:);zeros(1,MZ,NN2)];      
        
[adv_v3,adv_g3] = adject_RNL_kn(upw,vpw,wpw,gpw,umean,vmean,wmean,gmean,Lkn);

        
        [gi,vi,ui,wi] = solv_adj_vg_RK_kron_v2_k(gLt,vLt,(adv_g1+4*adv_g2+adv_g3)/6,(adv_v1+4*adv_v2+adv_v3)/6,2,1*h);
       
                        
%         %for qq=1:N+2
             vp=ifft(vi,[],2);
             up=ifft(ui,[],2);
             wp=ifft(wi,[],2);
             gp=ifft(gi,[],2);
%         %end
                    
        up=[zeros(1,MZ,NN2);up(2:end-1,:,:);zeros(1,MZ,NN2)];
        wp=[zeros(1,MZ,NN2);wp(2:end-1,:,:);zeros(1,MZ,NN2)];
        
        uLt=up;
        vLt=vp;
        wLt=wp;
        gLt=gp;
end
        
        uL(:,:,1:2:NN)=up;
        vL(:,:,1:2:NN)=vp;
        wL(:,:,1:2:NN)=wp;
        gL(:,:,1:2:NN)=gp;   
        
        uL(:,:,2:2:NN)=1i*up;
        vL(:,:,2:2:NN)=1i*vp;
        wL(:,:,2:2:NN)=1i*wp;
        gL(:,:,2:2:NN)=1i*gp;
        

 
end

