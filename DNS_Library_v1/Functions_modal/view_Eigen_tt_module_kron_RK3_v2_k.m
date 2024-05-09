function [uL,vL,wL,gL] = view_Eigen_tt_module_kron_RK3_v2_k(uLn,vLn,wLn,gLn,umn,vmn,wmn,gmn,NT,h,Lkn,nn,kn)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

global kLyap NLyap N NX MZ dt

umean=squeeze(umn(:,1,:));
vmean=squeeze(vmn(:,1,:));
wmean=squeeze(wmn(:,1,:));
gmean=squeeze(gmn(:,1,:));

% [uLn,vLn,wLn,gLn]=compact_Lyap(uLn,vLn,wLn,gLn);
% 
NN = sum(NLyap); NN2= NN/2;

uL = 0 * uLn;
vL = uL;
wL = uL;
gL = uL;

% uL=zeros(N+2,MZ,NN);
% vL=zeros(N+2,MZ,NN);
% wL=zeros(N+2,MZ,NN);
% gL=zeros(N+2,MZ,NN);
% %adv_v1=zeros(N+2,MZ,sum(NLyap)/2);
% %adv_g1=zeros(N+2,MZ,sum(NLyap)/2);

uLt=uLn(:,:,1:2:end);
vLt=vLn(:,:,1:2:end);
wLt=wLn(:,:,1:2:end);
gLt=gLn(:,:,1:2:end);

for it=2:NT

[adv_v1,adv_g1] = advect_RNL_kn(uLt,vLt,wLt,gLt,umean,vmean,wmean,gmean,Lkn);


        [gi,vi,ui,wi] = solv_vg_lyap_RK_kron_v2_k(gLt,vLt,adv_g1,adv_v1,1,1/2*h);
                        
%         %for qq=1:N+2
             vpw=ifft(vi,[],2);
             upw=ifft(ui,[],2);
             wpw=ifft(wi,[],2);
             gpw=ifft(gi,[],2);
%         %end
        
        upw=[zeros(1,MZ,NN2);upw(2:end-1,:,:);zeros(1,MZ,NN2)];
        wpw=[zeros(1,MZ,NN2);wpw(2:end-1,:,:);zeros(1,MZ,NN2)];
        
       
[adv_v2,adv_g2] = advect_RNL_kn(upw,vpw,wpw,gpw,umean,vmean,wmean,gmean,Lkn);

        [gi,vi,ui,wi] = solv_vg_lyap_RK_kron_v2_k(gLt,vLt,-adv_g1+2*adv_g2,-adv_v1+2*adv_v2,2,1*h);
        
 %         %for qq=1:N+2
             vpw=ifft(vi,[],2);
             upw=ifft(ui,[],2);
             wpw=ifft(wi,[],2);
             gpw=ifft(gi,[],2);
%         %end
        
        upw=[zeros(1,MZ,NN2);upw(2:end-1,:,:);zeros(1,MZ,NN2)];
        wpw=[zeros(1,MZ,NN2);wpw(2:end-1,:,:);zeros(1,MZ,NN2)];    

 [adv_v3,adv_g3] = advect_RNL_kn(upw,vpw,wpw,gpw,umean,vmean,wmean,gmean,Lkn);
      

        [gi,vi,ui,wi] = solv_vg_lyap_RK_kron_v2_k(gLt,vLt,(adv_g1+4*adv_g2+adv_g3)/6,(adv_v1+4*adv_v2+adv_v3)/6,2,1*h);       
       
                        
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

        
        uL(:,:,1:2:NN)=up;
        vL(:,:,1:2:NN)=vp;
        wL(:,:,1:2:NN)=wp;
        gL(:,:,1:2:NN)=gp;   
        
        uL(:,:,2:2:NN)=1i*up;
        vL(:,:,2:2:NN)=1i*vp;
        wL(:,:,2:2:NN)=1i*wp;
        gL(:,:,2:2:NN)=1i*gp;
        
        if rem(it,10)==0
        
        Lku=zeros(1,2*length(Lkn));
        Lku(1:2:end)=Lkn;Lku(2:2:end)=Lkn;
        
        [u1,v1,w1,g1]=decompact_Lyap_Fn(uL,vL,wL,gL,Lku);

        plot_mode(37,u1,v1,w1,g1,nn,kn,umn,vmn,wmn)
        drawnow
        end
        
end

%     if rem(Ln,2)==0
% 	
%     [uLb,vLb,wLb,gLb]=pick_vec_k(uL,vL,wL,gL,Ln-1,Lkiter);
%     
% 	[uL(:,:,Ln,Lkiter),vL(:,:,Ln,Lkiter),wL(:,:,Ln,Lkiter),gL(:,:,Ln,Lkiter)]=Lyap_translate_k(uLb,vLb,wLb,gLb);	
% 	
%     else    


% %[uLn,vLn,wLn,gLn]=decompact_Lyap(uLn,vLn,wLn,gLn);
% %[uL,vL,wL,gL]=decompact_Lyap(uL,vL,wL,gL);
% 
%  for Lk=1:kLyap
%     
%     for Ln=1:NLyap(Lk)
%         
%         [uLi,vLi,wLi,gLi]=pick_vec_k(uLn,vLn,wLn,gLn,Ln,Lk);
%         [uLt,vLt,wLt,gLt]=pick_vec_k(uL,vL,wL,gL,Ln,Lk);
%         
%         
%         if Ln==1
%             
%             Ei=Ener_k(uLi,vLi,wLi);
%             Ef=Ener_k(uLt,vLt,wLt);
%         
%             gr(Ln,Lk)=log(Ef/Ei)/(2*dt);
%             
%             [uLn(:,:,Ln,Lk),vLn(:,:,Ln,Lk),wLn(:,:,Ln,Lk),gLn(:,:,Ln,Lk)]=Lyap_norm_k(uLt,vLt,wLt,gLt,1);
%             
%         else
%             [uL0,vL0,wL0,gL0]=pick_vec_k(uL,vL,wL,gL,Ln,Lk);
%             
%             for itr=1:Ln-1
%                 
%                 [uLb,vLb,wLb,gLb]=pick_vec_k(uLn,vLn,wLn,gLn,itr,Lk);
%                 %[uL0,vL0,wL0,gL0]=pick_vec(uL,vL,wL,gL,Ln,Lk)
%                 
%                 P1=project_Lyap_k(uLb,vLb,wLb,uL0,vL0,wL0);
%                 
%                 [uLt,vLt,wLt,gLt] = project_out(uLb,vLb,wLb,gLb,uLt,vLt,wLt,gLt,P1);
%                 
%             end
%                 Ei=Ener_k(uLi,vLi,wLi);
%                 Ef=Ener_k(uLt,vLt,wLt);
%         
%                 gr(Ln,Lk)=log(Ef/Ei)/(2*dt);
%             
%                 [uLn(:,:,Ln,Lk),vLn(:,:,Ln,Lk),wLn(:,:,Ln,Lk),gLn(:,:,Ln,Lk)]=Lyap_norm_k(uLt,vLt,wLt,gLt,1);
%             
%         end
%        
%     end
%  end 
% 
% [uLn,vLn,wLn,gLn]=decompact_Lyap(uLn,vLn,wLn,gLn);



end

