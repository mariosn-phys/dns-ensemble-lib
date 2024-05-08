
%    Ensemble-Integration DNS Script of plane parallel flows (Couette/ Poiseuille) for MATLAB   
%    Copyright (C) 2024 Marios-Andreas Nikolaidis
%    Developed during the author's thesis at the University of Athens (NKUA)

%    This program is free software: you can redistribute it and/or modify
%    it under the terms of the GNU Affero General Public License as published
%    by the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.

%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU Affero General Public License for more details.

%    You should have received a copy of the GNU Affero General Public License
%    along with this program.  If not, see <https://www.gnu.org/licenses/>.

% DNS derived from the code used in Nikolaidis & Ioannou (2022)
% A pseudo-spectal Navier-Stokes solver for plane parallel Couette flow
% in wall-normal velocity/ vorticity formulation (Kim,Moin,Moser 1987 - 
% an example of a spectral mpi DNS has been developed by the Fluid Dynamics 
% group at UPM/ J.Jimenez et al.), that is capable of simulating turbulent flow
% dynamics.
% Differentiations are performed with the pseudo-spectral Chebyshev
% (Weideman & Reddy 2000) and Fourier (e.g Trefethen 2000) matrices. 
% Tested with NVIDIA gpus, enabled with igpu flag. Exceeding available GPU memory
% will result in crashes.
% Wherever possible during the time-stepping for-loops have been restructured to 
% matrix multiplication operations, which utilizes the build-in vectorization of 
% MATLAB (see also the DNS codes of Vuorinen & Keskinen 2015)

	
	
	
	
clear all;

global M N Re NX D1x D2x D1z D2z DYF D2F D2 D4 y gamma
global MZ k l dt xE yE zE A B L dx dz kkm llm DYkron DYkronf kkmf llmf ksolv msolv Nksolv Nmsolv 
 
global ICvkron1 ICvDvkron1 ICgkron1 ICggkron1
global ICvkron2 ICvDvkron2 ICgkron2 ICggkron2

% global ICvkron1c ICvDvkron1c ICgkron1c ICggkron1c
% global ICvkron2c ICvDvkron2c ICgkron2c ICggkron2c
% 
% global ICvCvRK1 ICvIDELvRK1 ICgRK1 ICgCgRK1
% global ICvCvRK2 ICvIDELvRK2 ICgRK2 ICgCgRK2

global S_mf S_mp Sol_m b1 b2
% global kLyap NLyap kkmm llmm DYkronm
global a b

%Function paths or run init_path in main folder
%addpath('../Functions/')
%addpath('../Functions_ensemble/')
%addpath('Functions_modal/')

solv=1; % 1 to calculate and save solver matrices, 0 to load stored 
%precalculated ones for the same parameters
psolv=1; % 1 enables parallel pool for solver matrices, 0 serial mode
npc=0; % number of parallel workers in pool for solvers, increases 
       %memory requirements during calculation

af=permute(linspace(6.2,8,10),[1 3 4 2]); af=1; 
% Nonlinear term modulation, 0 RNL 1 DNS // If used only as DNS 
% advect_NL_anl is changed to advect_NL for faster execution


init_ens=1; 
% initializes members of the ensemble with stochastic excitations on 
init_file=0;
% single state

stoch=1; 
% Stochastic forcing, 0 off 1 on.
% stochastic parameters for excitation   
% load('noise_structure_diffusion_eqall_313030.mat','FF')
% %load('noise_structure_dfreedom_realeqall_313030.mat','FF')
if stoch
    en2=35*0.5*10^-6;  % Noise energy level ~ 50 critical  
    %en2=0.3/16*70*10^-6;  %%% 1: 32, 2: 8
end

% F=2*N;
% kxband=1;
% kzband=MZ/3;
% FF=FF(:,1:F,:,:);   
% en2=64*70*10^-6;  %%% 1: 32, 2: 8

%%% GPU ---------------------------------
igpu=1; % 1 gpu on, 0 gpu off
if igpu
gpuDevice(1) % Assign to gpudevice
end
%%% -------------------------------------

%%% !!!!!!!!!!!!
field_path = 'Data/Re1000_n61_transition_ens/'  %% save and restart file path
diag_file  = 'diagnostics1'; %% save elementary diagnostics
%%% !!!!!!!!!!!!


Re=1000; % Reynolds number 

% Poiseuille (test) 'p' or Couette 'c'
modf='c';
it0=2; % First advance step
dt=0.02; % Plot and save intervals assume that 1/dt = integer. Some compatible 
         % choices for dt=[0.025,0.02,0.0125,0.01,0.008,0.00625,0.005] 
sdt=sqrt(dt); % Gudunov noise normalized step         
Ti=0; % Initial Time (should match time on restart file name)
Tf=2000; % Final Time
T=Ti:dt:Tf; 
NT=length(T);
g=dt/(2*Re); % CN coefficient

%%% File name %% Save file defined at the end of the time-stepping loop
fmt='%04.2f'; % Load and save time label format
%start_file=[field_path,'state_Re',num2str(Re),'_',num2str(T(1),fmt),'.mat']
start_filen=['state_Re',num2str(Re),'_',num2str(T(1),fmt),'.mat']

Fnn = 1:10; % List of ensemble members
tsav=10; % interval of saves
tplot=10; % plot basic diagnostics
tp = round(1/dt);
Fn = length(Fnn);

if init_ens == 1 %% create folders for ensemble member data
for ff = Fnn, mkdir([field_path,'N',num2str(ff,'%03.f')]), end
end

% Exclude ensemble members
% load([field_path,'lam2.mat'],'lam1')
% load([field_path,'lam2n.mat'],'lam2')
% load([field_path,'lam3n.mat'],'lam3')
% load([field_path,'lam4n.mat'],'lam4')
% load([field_path,'lam5n.mat'],'lam5')
% %load([field_path,'lam6n.mat'],'lam6')
% Fnn = find(not(sum(Fnn==lam1)+sum(Fnn==lam2)+sum(Fnn==lam3)+sum(Fnn==lam4)+sum(Fnn==lam5)));
% %Fnn = find(not(sum((Fnn)==lam1)));
% Fn = length(Fnn);

% Fundamental wavenumbers
%Re3000 and Re3250
a=2/2;
b=2/1.2;

%Re600 and Re2250
%a=2/1.75;
%b=2/(1.2);


%%%% NX,MZ 12,24,48 preferably increment at steps of 12 point, also works
%%%% with increments of 6. N is odd for centerline point.
%%%% Boundary points not counted. 
	
NX=72;
MZ=72;
N=61;

A=2*pi/a;dx=A/NX;x=-dx+dx*(1:NX)';xE=[x;A]; 
B=2*pi/b;dz=B/MZ;z=-dz+dz*(1:MZ)';zE=[z;B];
L=2;

%ksolv=1:3; 
ksolv=1:NX/3+1;

Nksolv = length(ksolv);
msolv=1:MZ/3+1; 

Nmsolv = length([msolv MZ-flip(msolv(2:end))+2]); 

%save([field_path,'parameters_Re',num2str(Re),'_n',num2str(N),'.mat'],'a','b','Re','N','NX','MZ','dt')
%Parameters for 600, 2250 and 3000 can be found in their respective folders '../Data/Re*/

% Form Differentiation matrices
% x matrix (see eg Trefethen 2000)
    column=[0 .5*(-1).^(1:NX-1).*cot((1:NX-1)*a*dx/2)];
    D1x=a*toeplitz(column,column([1 NX:-1:2]));
    D2x=D1x^2;
% z matrix    
    column=[0 .5*(-1).^(1:MZ-1).*cot((1:MZ-1)*b*dz/2)];
    D1z=b*toeplitz(column,column([1 MZ:-1:2]));
    D2z=D1z^2;
  
% y matrix (see Weideman and Reddy 2000)
    DYF=cheb(N+1);DYF=flip(flip(DYF,1),2);
    DY=DYF(2:end-1,2:end-1);
    DYkron=kron(speye((Nksolv)*(Nmsolv)),DY);

    D2F=DYF^2;
    D2=D2F(2:end-1,2:end-1);

    [y,D4]=cheb4c(N+2);D4=flip(flip(D4,1),2);
    y=-y;
    yE=[-1;y;1];
    
%-------%
%Required only to plot the instantaneous flow fields
%[xx,yy]=meshgrid(x,y);
%[xxE,yyE]=meshgrid(x,yE);
%-------%
    
%

M=[0:(NX/2-1) 0 (1-NX/2):(-1)];
    k=2*pi*M/A;
        
P=[0:(MZ/2-1) 0 (1-MZ/2):(-1)];
    l=2*pi*P/B;     
    
    if modf=='p'
    U1=4/3*(1-yE.^2);gm=-8/3/Re; %Poiseuille (constant pressure) 
    elseif modf=='c'
    U1=yE;gm=0; %Couette
    elseif modf=='z'
    U1=yE*0;gm=-2/Re*0;    
    end
    b1=U1(1);b2=U1(end);
  
    Uback=repmat(U1,[1,NX,MZ]);
    gamma=gm*ones(N+2,1);


if stoch
FF=NoiseStructure2_cheb(Re,N,NX,MZ,a,b,0)
%load('noise_structure_diffusion_eqall_617272.mat','FF')
%load('noise_structure_dfreedom_realeqall_313030.mat','FF')
%load('noise_short_kron_structure_diffusion_eqall_313030.mat','FFkron')
%load('noise_structure_diffusion_eqall_717272.mat','FF')

F=2*N;
kxband=1:2;Nk=length(kxband);
kzband=1:MZ/3+1;Nm=length([kzband MZ-flip(kzband(2:end))+2]);
FF=en2*FF(:,1:F,:,:);  FF(:,:,:,1)=0*FF(:,:,:,1); 



%save('noise_short_kron_structure_diffusion_eqall_313030.mat','FFkron')
[FFkron,DYkronf,kkmf,llmf]=solvers_noise(FF,DY,Nk,Nm,kxband,kzband);


end
    


%----------Build Solvers  
if solv==1
    if psolv==1    
    parpool(npc) 

    [ICvkron1,ICvDvkron1,ICgkron1,ICggkron1,~,~,~,~,S_mf(:,:,1),S_mp(:,:,1),Sol_m(:,:,1),kkm,llm]=parsolvers(g,1/2);
      
    [ICvkron2,ICvDvkron2,ICgkron2,ICggkron2,~,~,~,~,S_mf(:,:,2),S_mp(:,:,2),Sol_m(:,:,2),~,~]=parsolvers(g,1);

    delete(gcp)
    else
  
    [ICvkron1,ICvDvkron1,ICgkron1,ICggkron1,~,~,~,~,S_mf(:,:,1),S_mp(:,:,1),Sol_m(:,:,1),kkm,llm]=solvers(g,1/2);
      
    [ICvkron2,ICvDvkron2,ICgkron2,ICggkron2,~,~,~,~,S_mf(:,:,2),S_mp(:,:,2),Sol_m(:,:,2),~,~]=solvers(g,1);   
    end

save([field_path,'Solvers.mat'],'ICvkron1','ICvkron2','ICvDvkron1','ICvDvkron2','ICgkron1','ICgkron2','ICggkron1','ICggkron2','S_mf','S_mp','Sol_m','kkm','llm','-v7.3')
else 
%----------Load Solvers (only if calculated with same parameteres!!!)
load([field_path,'Solvers.mat'])
display('Solvers loaded')
end

display('Solver initialization complete')

% Transfer Variables to GPU 
if igpu==1
      ICvkron1=gpuArray(ICvkron1);
      ICvDvkron1=gpuArray(ICvDvkron1);
      ICgkron1=gpuArray(ICgkron1);
      ICggkron1=gpuArray(ICggkron1);
      kkm=gpuArray(kkm);
      llm=gpuArray(llm);
      DYkron=gpuArray(DYkron);

      ICvkron2=gpuArray(ICvkron2);
      ICvDvkron2=gpuArray(ICvDvkron2);
      ICgkron2=gpuArray(ICgkron2);
      ICggkron2=gpuArray(ICggkron2);
      display('GPU upload complete')
      
      if stoch
          FFkron=gpuArray(FFkron);
          DYkronf=gpuArray(DYkronf);
      end 
end
%----------Build Solvers End   

% Restart from ensemble folder

if init_ens == 0

[u0,v0,w0,g0]=input_ensemble(field_path,start_filen,Fnn);

else
    
if init_file == 1
    %    Load init and transform to physical space
[vi,gi,UP1,WP1]=read_from_disk_compact([field_path,start_filen]);
[u0,v0,w0,g0] = make_uw(gather(vi),gather(gi),gather(UP1),gather(WP1),b1,b2);

else    
 u0 = Uback;
 v0 = 0 * u0;
 w0 = 0 * u0;
 g0 = 0 * u0;
end

    u0 = repmat(u0,[1 1 1 Fn]);
    v0 = repmat(v0,[1 1 1 Fn]);
    w0 = repmat(w0,[1 1 1 Fn]);
    g0 = repmat(g0,[1 1 1 Fn]);

end

    if igpu==1
         v0=gpuArray(v0);
         u0=gpuArray(u0);
         w0=gpuArray(w0);
         g0=gpuArray(g0);
    end

 if stoch
 %for it=1:10
    
[un,vn,wn,~]=noise_short_kron_en_kx(FFkron,F,en2,kxband,kzband,Fn,0*u0);
Ei=EnerFn(un,vn,wn);

% Rescale noise
ffsc=(en2/mean(Ei))^(1/2);
FFkron=ffsc*FFkron;
 
%   [u0,v0,w0,g0]=noise_kron_en_kx(FFkron,F,en2,Fn);
 [un,vn,wn,gn]=noise_short_kron_en_kx(FFkron,F,en2,kxband,kzband,Fn,0*u0);
    
 %end

 u0 = u0 + sdt*un;
 v0 = v0 + sdt*vn;
 w0 = w0 + sdt*wn;
 g0 = g0 + sdt*gn;
 
 
 [u0,v0,w0,g0] = replace_mean_Fn(u0,v0,w0,g0,Fn);

end
 

% Mean flow        
%    p1 = 5 % Nth member
%    vmean1=repmat(mean(v0(:,:,:,p1),2),[1 NX 1 1]); % vmean=0*vmean;
%    gmean1=repmat(mean(g0(:,:,:,p1),2),[1 NX 1 1]); % gmean=0*gmean;
%    umean1=repmat(mean(u0(:,:,:,p1),2),[1 NX 1 1]);
%    wmean1=repmat(mean(w0(:,:,:,p1),2),[1 NX 1 1]);
    
%     figure(91);
%     pcolor(squeeze(umean1(:,1,:)))
    

    
    Efm=zeros(Fn,NT);CFL=Efm;O_bot=Efm;O_top=Efm;
% 
%    if stoch
%  
%     [un,vn,wn,gn]=noise_en_kx(FF,F,en2,1,kzband,Fn);
%  
%  	v0=v0+sdt*vn;
%  	u0=u0+sdt*un;
%  	w0=w0+sdt*wn;
%  	g0=g0+sdt*gn;
%      
%    end


    Efm(:,1)=gather(EnerFn(u0-Uback,v0,w0)); % Energy of deviations from the laminar flow
    %CFL(:,1)=gather((max(u0(:))/dx+max(abs(max(v0(2:end,:),[],2)./diff(yE)))+max(w0(:))/dz)*dt);
	
    CFL(:,1)=gather((max(reshape(u0,[(N+2)*NX*MZ,Fn]))/dx...
        +max(abs(max(reshape(permute(v0(2:end,:,:,:),[1 4 2 3]),[(N+1),Fn,NX*MZ]),[],3)./repmat(diff(yE),[1 Fn])))...
        +max(reshape(w0,[(N+2)*NX*MZ,Fn]))/dz)*dt);
   

    % Shear at top and bottom walls	
    O_bot(:,1)=gather(DYF(1,:)*permute(mean(mean(u0,3),2),[1 4 2 3])); 
    O_top(:,1)=gather(DYF(end,:)*permute(mean(mean(u0,3),2),[1 4 2 3]));

% 
% sum = 0;
% for it=1:100
%     tic;
%     [advv1,advg1] = advect_NLn(ue,ve,we,ge,20);   
%     toc
%     %sum = sum + toc
% end

%%% speed up of advect by ~ 3 times
    
    %%
    tic;
    tdt=[];
for it=it0:NT
    
   %if rem(it,125)==1
   tic;
   %end
   
    % u init
    
    up=u0;
    vp=v0;
    wp=w0;
    gp=g0;
    
%     if af==1
     [advv1,advg1] = advect_NLn(up,vp,wp,gp,Fn);   
%     else
%     [advv1,advg1] = advect_NLn_anl(up,vp,wp,gp,af,Fn);
%     end
    
    [advU1,advW1] = advect_mean(up,vp,wp,Fn);

    [ui,vi,wi,gi,UP1,WP1] = rkstepn(u0,v0,w0,g0,advg1,advv1,advU1,advW1,1,1/2,Fn);
    
    v1=ifft2_cuben(vi,Fn);
    u1=ifft2_cuben(ui,Fn);
    w1=ifft2_cuben(wi,Fn);
    g1=ifft2_cuben(gi,Fn);

    [u1,w1]=boundary_uw(u1,w1,UP1,WP1,Fn);
   
    
    [u1,v1,w1,g1] = replace_mean_Fn(u1,v1,w1,g1,Fn);
    
    % u step 1
    
    up=u1;
    vp=v1;
    wp=w1;
    gp=g1;
    
    
%     if af==1
     [advv2,advg2] = advect_NLn(up,vp,wp,gp,Fn);   
%     else
%    [advv2,advg2] = advect_NLn_anl(up,vp,wp,gp,af,Fn);
%     end 
        
       
    [advU2,advW2] = advect_mean(up,vp,wp,Fn);

    [ui,vi,wi,gi,UP2,WP2] = rkstepn(u0,v0,w0,g0,2*advg2-advg1,2*advv2-advv1,2*advU2-advU1,2*advW2-advW1,2,1,Fn);

    v2=ifft2_cuben(vi,Fn);
    u2=ifft2_cuben(ui,Fn);
    w2=ifft2_cuben(wi,Fn);
    g2=ifft2_cuben(gi,Fn);

    [u2,w2]=boundary_uw(u2,w2,UP2,WP2,Fn);
    
    [u2,v2,w2,g2] = replace_mean_Fn(u2,v2,w2,g2,Fn);

    
    % u step 2
    
    up=u2;
    vp=v2;
    wp=w2;
    gp=g2;
    
%     if af==1
     [advv3,advg3] = advect_NLn(up,vp,wp,gp,Fn);   
%     else
%    [advv3,advg3] = advect_NLn_anl(up,vp,wp,gp,af,Fn);
%     end

    [advU3,advW3] = advect_mean(up,vp,wp,Fn);
        
    [ui,vi,wi,gi,UP3,WP3] = rkstepn(u0,v0,w0,g0,(advg1+4*advg2+advg3)/6,(advv1+4*advv2+advv3)/6,(advU1+4*advU2+advU3)/6,(advW1+4*advW2+advW3)/6,2,1,Fn);

    v0=ifft2_cuben(vi,Fn);
    u0=ifft2_cuben(ui,Fn);
    w0=ifft2_cuben(wi,Fn);
    g0=ifft2_cuben(gi,Fn);

    % u final
    
    [u0,w0]=boundary_uw(u0,w0,UP3,WP3,Fn);
    
    [u0,v0,w0,g0] = replace_mean_Fn(u0,v0,w0,g0,Fn);

            
%     [~,u0]=kx_filter(u0,srg);
%     [~,v0]=kx_filter(v0,srg);
%     [~,w0]=kx_filter(w0,srg);
%     [~,g0]=kx_filter(g0,srg);

%     if stoch
%     [un,vn,wn,gn]=noise_en_kx(FF,F,en2,1,kzband,Fn);
% 
%     un=gpuArray(un);
%     vn=gpuArray(vn);
%     wn=gpuArray(wn);
%     gn=gpuArray(gn);
% 
%     Gud step
% 
%     v0=v0+sdt*vn;
% 	  u0=u0+sdt*un;
% 	  w0=w0+sdt*wn;
% 	  g0=g0+sdt*gn;
%     end

    if stoch
 
    %[un,vn,wn,gn]=noise_en_kx(FF,F,en2,1,kzband,Fn);
    [un,vn,wn,gn]=noise_short_kron_en_kx(FFkron,F,en2,kxband,kzband,Fn,0*u0);
 
 	v0=v0+sdt*vn;
 	u0=u0+sdt*un;
 	w0=w0+sdt*wn;
 	g0=g0+sdt*gn;
     
    end


    Efm(:,it)=gather(EnerFn(u0-Uback,v0,w0)); 
    
    CFL(:,it)=gather((max(reshape(u0,[(N+2)*NX*MZ,Fn]))/dx...
        +max(abs(max(reshape(permute(v0(2:end,:,:,:),[1 4 2 3]),[(N+1),Fn,NX*MZ]),[],3)./repmat(diff(yE),[1 Fn])))+max(reshape(w0,[(N+2)*NX*MZ,Fn]))/dz)*dt);
   
  %  O_bot(:,it)=gather(DYF(1,:)*UP3);
  %  O_top(:,it)=gather(DYF(end,:)*UP3);


    O_bot(:,it)=gather(DYF(1,:)*permute(mean(mean(u0,3),2),[1 4 2 3])); 
    O_top(:,it)=gather(DYF(end,:)*permute(mean(mean(u0,3),2),[1 4 2 3]));
    
    

    if rem(T(it),tplot)==0
    figure(99);
    plot(UP3,yE,'+',U1,yE,'r')
    xlabel('U'),ylabel('y')
%     figure(105);
%     contourf(x,yE,v0(:,:,1),20)
%     xlabel('x');ylabel('y');title('v')
    figure(111);clf
    plot(T(1:tp:it),Efm(:,1:tp:it))

    xlabel('T');ylabel('E')
    figure(112);clf
    plot(T(1:tp:it),CFL(:,1:tp:it))
    xlabel('T');ylabel('CFL')
    figure(113);clf
    	if modf=='c'
    	plot(T(1:tp:it),sqrt(Re*(O_bot(:,1:tp:it)+O_top(:,1:tp:it))/2))
    	elseif modf=='p'
    	plot(T(1:tp:it),sqrt(Re*(O_bot(:,1:tp:it)-O_top(:,1:tp:it))/2))
	    end
    xlabel('T');ylabel('Re_{\tau}')
    figure(1111);
    bar(1:Fn,Efm(:,it))
    ylabel('E');xlabel('n')
    title(['n_{lam}=',num2str(sum(Efm(:,it)<=10^-3))])
    
%     figure(122);
%     contourf(z,yE,squeeze(mean(u0,2)))
%     xlabel('z');ylabel('y');title('U(y,z,t)')
%     drawnow
%     figure(123);
%     contourf(z,yE,squeeze(mean(u0-repmat(mean(mean(u0,3),2),[1 NX MZ]),2)))
%     xlabel('z');ylabel('y');title('streak')
    drawnow
    end
    

   
    %if rem(it,125)==0
    tdt=[tdt toc]; % time step bench
    %end
     
    if rem(T(it),tsav)==0
    
    vi=fft2_cuben(v0,Fn);
    %u0=ifft2_cube(ui);
    %w0=ifft2_cube(wi);
    gi=fft2_cuben(g0,Fn);   
     
    UP3=mean(mean(u0,3),2);
    WP3=mean(mean(w0,3),2);
        
    % Write Fourier coefficients to disk on separate files and folders    
    output_ensemble(gather(vi),gather(gi),gather(UP3),gather(WP3),T(it),field_path,['state_Re',num2str(Re),'_',num2str(T(it),fmt),'.mat'],Fnn)


    end

    if rem(T(it),3*tsav)==0

    save([field_path,diag_file,'.mat'],'Efm','CFL','O_bot','O_top','tdt','af','T')

    end

    
    % Code Breakdown check
    %if max(abs(u0(:)))>=1.2*abs(max(U1))
    %   display('error')
    %	break
    %end
    
end
    
if init_ens == 1 %% save as initial field
output_ensemble(gather(vi),gather(gi),gather(UP3),gather(WP3),0,field_path,['state_Re',num2str(Re),'_',num2str(0,fmt),'.mat'],Fnn)
end
