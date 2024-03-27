
%    Ensemble-Integration DNS Script of plane parallel flows (Couette/ Poiseuille) for MATLAB   
%    Copyright (C) 2023 Marios-Andreas Nikolaidis
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

% Base DNS code used in Nikolaidis & Ioannou (2022)
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
global MZ k l dt xE yE zE A B L dx dz kkm llm DYkron DYkronf kkmf llmf
 
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

%Function paths or add main folder ('../') and subfolders to path
addpath('../Functions/')
%addpath('Functions_modal/')

solv=1; % 1 to calculate and save solver matrices, 0 to load stored 
%precalculated ones for the same parameters
psolv=0; % 1 enables parallel pool for solver matrices, 0 serial mode
npc=2; % number of parallel workers in pool for solvers, increases 
       %memory requirements during calculation

af=1;  
% Nonlinear term modulation, 0 RNL 1 DNS // If used only as DNS 
% advect_NL_anl is changed to advect_NL for faster execution


init_ens=0; 
% initializes members of the ensemble with stochatstic excitations on 
% single state

stoch=1; 
% Stochastic forcing, 0 off 1 on.
% stochastic parameters for excitation   


%%% GPU ---------------------------------
igpu=1; % 1 gpu on, 0 gpu off
if igpu
gpuDevice(2) % Assign to gpudevice
end
%%% -------------------------------------

%%% !!!!!!!!!!!!
field_path = '../Data/Re400_n31_nsblue/'  %% save and restart file path
diag_file  = 'diagnostics2'; %% save elementary diagnostics
%%% !!!!!!!!!!!!

Re=400; % Reynolds number 

% Poiseuille (test) 'p' or Couette 'c'
modf='c';
it0=2; % First advance step
dt=0.025; % Plot and save intervals assume that 1/dt = integer. Some compatible 
         % choices for dt=[0.025,0.02,0.0125,0.01,0.008,0.00625,0.005] 
sdt=sqrt(dt); % Gudunov noise normalized step         
Ti=70000; % Initial Time (should match time on restart file name)
Tf=80000; % Final Time
T=Ti:dt:Tf; 
NT=length(T);
g=dt/(2*Re); % CN coefficient

%%% File name %% Save file defined at the end of the time-stepping loop
fmt='%04.2f'; % Load and save time label format
%start_file=[field_path,'state_Re',num2str(Re),'_',num2str(T(1),fmt),'.mat']
start_filen=['state_Re',num2str(Re),'_',num2str(T(1),fmt),'.mat']

Fnn = 1:100; % List of ensemble members
tsav=1; % interval of saves
tplot=50; % plot basic diagnostics
tp = round(1/dt);
Fn = length(Fnn);

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
%a=2/2;
%b=2/1;

%Re600 and Re2250
a=2/1.75;
b=2/(1.2);


%%%% NX,MZ 12,24,48 preferably increment at steps of 12 point, also works
%%%% with increments of 6. N is odd for centerline point.
%%%% Boundary points not counted. 
	
NX=30;
MZ=30;
N=31;

A=2*pi/a;dx=A/NX;x=-dx+dx*(1:NX)';xE=[x;A]; 
B=2*pi/b;dz=B/MZ;z=-dz+dz*(1:MZ)';zE=[z;B];
L=2;

if stoch
load('noise_structure_diffusion_eqall_313030.mat','FF')
%load('noise_structure_dfreedom_realeqall_313030.mat','FF')
load('noise_short_kron_structure_diffusion_eqall_313030.mat','FFkron')

F=2*N;
kxband=1:3;Nk=length(kxband);
kzband=1:MZ/3;Nm=length([kzband MZ-flip(kzband(2:end))+2]);
FF=FF(:,1:F,:,:);  FF(:,:,:,1)=0*FF(:,:,:,1); 
en2=64*70*10^-6;  %%% 1: 32, 2: 8



%save('noise_short_kron_structure_diffusion_eqall_313030.mat','FFkron')


end

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
    DYkron=kron(speye((NX/3)*(2*MZ/3-1)),DY);

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

    [FFkron,DYkronf,kkmf,llmf]=solvers_noise(FF,DY,Nk,Nm,kxband,kzband);


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

    %save([field_path,'Solvers.mat'],'ICvkron1','ICvkron2','ICvDvkron1','ICvDvkron2','ICgkron1','ICgkron2','ICggkron1','ICggkron2','S_mf','S_mp','Sol_m','kkm','llm','-v7.3')
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
          FFkrong=gpuArray(FFkron);
          DYkronf=gpuArray(DYkronf);
      end
          
end
%----------Build Solvers End   

% Restart from ensemble folder

%[u0,v0,w0,g0]=input_ensemble(field_path,start_filen,Fnn);

% Mean flow        
%    p1 = 5 % Nth member
%    vmean1=repmat(mean(v0(:,:,:,p1),2),[1 NX 1 1]); % vmean=0*vmean;
%    gmean1=repmat(mean(g0(:,:,:,p1),2),[1 NX 1 1]); % gmean=0*gmean;
%    umean1=repmat(mean(u0(:,:,:,p1),2),[1 NX 1 1]);
%    wmean1=repmat(mean(w0(:,:,:,p1),2),[1 NX 1 1]);
    
%     figure(91);
%     pcolor(squeeze(umean1(:,1,:)))
    
%    if igpu==1
%         v0=gpuArray(v0);
%         u0=gpuArray(u0);
%         w0=gpuArray(w0);
%         g0=gpuArray(g0);
%    end
    
%    Efm=zeros(Fn,NT);CFL=Efm;O_bot=Efm;O_top=Efm;
%%
    
   if stoch
 for it=1:10
    
%   [u0,v0,w0,g0]=noise_kron_en_kx(FFkron,F,en2,Fn);
   [u0,v0,w0,g0]=noise_short_kron_en_kx(FFkrong,F,en2,kxband,kzband,Fn);
    
 end
 
 for it=1:10
    
    [un,vn,wn,gn]=noise_en_kx(FF,F,en2,kxband,kzband,Fn);
    
 end
 
 
 
 	v0=v0+sdt*vn;
 	u0=u0+sdt*un;
 	w0=w0+sdt*wn;
 	g0=g0+sdt*gn;
     
  end


