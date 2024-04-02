function [uLn,vLn,wLn,gLn,Lkn,grold,rtold]=initialize_dist(field_path,start_mode,noise_file,cont_old);

global b1 b2 kLyap NLyap NX MZ N

display('Generating random initial fields')
load(noise_file,'FF')
%load('noise_structure_dfreedom_eqall_5372120.mat','FF')
    
rtold=0;
grold=[];
F=2*N;
kxband=NX/3;kxband=1;
kzband=MZ/3;
FF=FF(:,1:F,:,:);   

Lkiter=0;
for Lk=kLyap
   Lkiter=Lkiter+1; 
    
   for Ln=1:NLyap(Lkiter) 
       
     [u_L,v_L,w_L,g_L]=noise_kx(FF,F,1,Lk,kzband);
 
     uL(:,:,:,Ln,Lkiter)=u_L;
     vL(:,:,:,Ln,Lkiter)=v_L;
     wL(:,:,:,Ln,Lkiter)=w_L;
     gL(:,:,:,Ln,Lkiter)=g_L;
     
   end
   
end

[uL,vL,wL,gL] = orth_Lyap(uL,vL,wL,gL); %initial normalization of disturbances

if cont_old==1

display('Replacing with previous disturbances')    
    
[vi,gi,UP,WP,grold,NL,KL,rtold]=read_from_disk_compact_ensemble_opti([field_path,start_mode]);
[u1,v1,w1,g1]=make_uw_ens_Lyap(vi,gi,UP,WP,0*b1,0*b2,NL,KL);

%% %init if 1
  Lkiter=0;
  for Lk=1:min(length(KL),length(kLyap))
    Lkiter=Lkiter+1;
    %
     for Ln=1:min(NL(Lkiter),NLyap(Lkiter))

       Lk1=find(KL==kLyap(Lkiter));

       if Lk1, 
       uL(:,:,:,Ln,Lk)=u1(:,:,:,Ln,Lk1);
       vL(:,:,:,Ln,Lk)=v1(:,:,:,Ln,Lk1);
       wL(:,:,:,Ln,Lk)=w1(:,:,:,Ln,Lk1);
       gL(:,:,:,Ln,Lk)=g1(:,:,:,Ln,Lk1);
       else
       end
 %
     end
 %
  end
 
end

[uLn,vLn,wLn,gLn,Lkn]=compact_Lyap_Fn(uL,vL,wL,gL);
