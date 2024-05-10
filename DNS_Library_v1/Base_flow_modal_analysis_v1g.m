
%    Power method eigenanalysis of plane parallel flows (Couette/ Poiseuille) 
%    for MATLAB   
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

clear all;

global M N Re NX D1x D2x D1z D2z DYF D2F D2 D4 y gamma
global MZ k l dt xE yE zE A B L dx dz kkm llm kkmm llmm DYkron DYkronm msolv Nmsolv
 
global ICvkron1 ICvDvkron1 ICgkron1 ICggkron1
global ICvkron2 ICvDvkron2 ICgkron2 ICggkron2

global ICvkron1c ICvDvkron1c ICgkron1c ICggkron1c
global ICvkron2c ICvDvkron2c ICgkron2c ICggkron2c

global ICvCvRK1 ICvIDELvRK1 ICgRK1 ICgCgRK1
global ICvCvRK2 ICvIDELvRK2 ICgRK2 ICgCgRK2

global S_mf S_mp Sol_m b1 b2 igpu
global kLyap NLyap
global a b z x
global xx_yx yy_yx yy_yz zz_yz xx_xz zz_xz

%load kaw_parameters.mat
%addpath('Functions/')
%addpath('Functions_modal/')

field_path='Data/Re1500_n61_alpha_ens/';%field_path='';

%%% GPU ---------------------------------
igpu=0; % 1 gpu on, 0 gpu off
if igpu
gpuDevice(1) % Assign to gpudevice
end
%%% -------------------------------------

%%%% Optimization parameters 
kLyap=1:3; % disturbance k_x range
NLyap=[10;10;10]; % N disturbances per k_x (these include streamwise translations) 

%Iteration settigns
cstep=10; % number of full iterations
nstep=101; % timesteps before normalization
h=10; %dt speed factor


Re=1500;%Re=dum(2);
% Poiseuille or Couette 'p' or 'c'
mod='c';
it0=2;
dt=0.0125;
Tfopt=15;
T=0:h*dt:Tfopt;
NT=length(T);
g=dt/(2*Re);



%%% Load initial files
init_file=1; Tinit=335;
ty = "file";
start_filen=['state_Re',num2str(Re),'_',num2str(Tinit(1),'%04.2f'),'.mat']
%ty = "tseries"
%file_series=['mean_states_1500_N1.mat']
%tt = 47
fmt = '%04.2f';
cont_old=0; % initializes disturbances from previous fields, 
            % if 0 starts with random fields 
%noise_file='noise_structure_diffusion_eqall_617272.mat';
            %Noi            
%if cont_old == 0 % continues from old disturbance file
start_mode=['modes_T',num2str(Tinit(1),'%04.2f'),'_k123.mat']
start_adjoint=['adjoi_T',num2str(Tinit(1),'%04.2f'),'_k123.mat']
%end


% Fundamental wavenumbers
%Re3000 and Re3250
a=2/2;
b=2/1;

%Re600 and Re2250
%a=2/1.75;
%b=2/(1.2);


%%%% NX,MZ 12,24,48 preferably increment at steps of 12 point, also works
%%%% with increments of 6. N is odd for centerline point.
%%%% Boundary points not counted. 

A=2*pi/a;NX=72;dx=A/NX;x=-dx+dx*(1:NX)';xE=[x;A]; 
B=2*pi/b;MZ=72;dz=B/MZ;z=-dz+dz*(1:MZ)';zE=[z;B];
L=2;

N=61;

msolv=1:MZ/3+1; 

Nmsolv = length([msolv MZ-flip(msolv(2:end))+2]); 

% x matrix
    column=[0 .5*(-1).^(1:NX-1).*cot((1:NX-1)*a*dx/2)];
    D1x=a*toeplitz(column,column([1 NX:-1:2]));
    D2x=D1x^2;
% z matrix    
    column=[0 .5*(-1).^(1:MZ-1).*cot((1:MZ-1)*b*dz/2)];
    D1z=b*toeplitz(column,column([1 MZ:-1:2]));
    D2z=D1z^2;
  
% y matrix 
    DYF=cheb(N+1);DYF=flip(flip(DYF,1),2);
    DY=DYF(2:end-1,2:end-1);
    DYkron=kron(speye((NX/3)*Nmsolv),DY);
    DYkronm=kron(speye(sum(NLyap/2)*Nmsolv),DY);

    D2F=DYF^2;
    D2=D2F(2:end-1,2:end-1);

    [y,D4]=cheb4c(N+2);D4=flip(flip(D4,1),2);
    y=-y;
    yE=[-1;y;1];
    
%-------%
%Required only to plot the instantaneous flow fields
[xx,yy]=meshgrid(x,y);
[xxE,yyE]=meshgrid(x,yE);
[yy_yx,xx_yx]=meshgrid(yE,x);
[yy_yz,zz_yz]=meshgrid(yE,z);
[zz_xz,xx_xz]=meshgrid(z,x);
%-------%
    

M=[0:(NX/2-1) 0 (1-NX/2):(-1)];
    k=2*pi*M/A; 
P=[0:(MZ/2-1) 0 (1-MZ/2):(-1)];
    l=2*pi*P/B;  
  
    if mod=='p'
    U1=(1-yE.^2);gm=-2/Re; %Poiseuille
    elseif mod=='c'
    U1=yE;gm=0; %Couette
    elseif mod=='z'
    U1=yE*0;gm=-2/Re*0;    
    end
    b1=U1(1);b2=U1(end);

    %Initial Conditions vortex+optimal n_x=1,n_z=1
  
    Uback=repmat(U1,[1,NX,MZ]);
    gamma=gm*ones(N+2,1);

   %%
    
%----------Build Solvers   
[ICvCvRK1,ICvIDELvRK1,ICgRK1,ICgCgRK1,~,~,~,kkmm,llmm]=solvers_module(g,1/2*h);
% 


[ICvCvRK2,ICvIDELvRK2,ICgRK2,ICgCgRK2,~,~,~,~,~]=solvers_module(g,1*h);
%    
if igpu == 1
ICvCvRK1 = gpuArray(ICvCvRK1 );
ICvIDELvRK1 = gpuArray(ICvIDELvRK1 );
ICgRK1 = gpuArray(ICgRK1 );
ICgCgRK1 = gpuArray(ICgCgRK1 );

ICvCvRK2 = gpuArray(ICvCvRK2 );
ICvIDELvRK2 = gpuArray(ICvIDELvRK2 );
ICgRK2 = gpuArray(ICgRK2 );
ICgCgRK2 = gpuArray(ICgCgRK2 );

kkmm = gpuArray(kkmm);
llmm = gpuArray(llmm);
DYkronm = gpuArray(DYkronm);
end
%----------Build Solvers End   

%% Load disturbances

FF=NoiseStructure2_cheb(Re,N,NX,MZ,a,b,0);    
[uLn,vLn,wLn,gLn,Lkn,gmold,rtold]=initialize_dist(field_path,start_mode,FF,cont_old);
[uLj,vLj,wLj,gLj,Lkj,gaold,rtold]=initialize_dist(field_path,start_adjoint,FF,cont_old);

    if igpu == 1
    uLn = gpuArray(uLn);
    vLn = gpuArray(vLn);
    wLn = gpuArray(wLn);
    gLn = gpuArray(gLn);

    uLj = gpuArray(uLj);
    vLj = gpuArray(vLj);
    wLj = gpuArray(wLj);
    gLj = gpuArray(gLj);
    end

%% Load base flow 
if ty == "file"
if init_file == 1
    %    Load init and transform to physical space
[vi,gi,UP1,WP1]=read_from_disk_compact([field_path,start_filen]);
[u0,v0,w0,g0] = make_uw(gather(vi),gather(gi),gather(UP1),gather(WP1),b1,b2);

else
    % Laminar flow
    u0=Uback;
    v0=0*u0;
    w0=0*u0;
    g0=0*u0;
end    

    vmean1=repmat(mean(v0(:,:,:,1),2),[1 NX 1 1]); % vmean=0*vmean;
    gmean1=repmat(mean(g0(:,:,:,1),2),[1 NX 1 1]); % gmean=0*gmean;
    umean1=repmat(mean(u0(:,:,:,1),2),[1 NX 1 1]);
    wmean1=repmat(mean(w0(:,:,:,1),2),[1 NX 1 1]);
 
elseif ty == "tseries" 

    % load from timeseries
    load(file_series)

   %for 
   %tt=47:47;
    
    vmean1=repmat(permute(vmean(:,:,tt),[1 3 2]),[1 NX 1 1]); % vmean=0*vmean;
    umean1=repmat(permute(umean(:,:,tt),[1 3 2]),[1 NX 1 1]);
    wmean1=repmat(permute(wmean(:,:,tt),[1 3 2]),[1 NX 1 1]);
    gmean1=difZ_F(umean1,1)-difX_F(wmean1,1); % gmean=0*gmean;
    
   %end


    % load from velocity snapshot
end

    figure(91);
    pcolor(squeeze(umean1(:,1,:)))
    drawnow
    pause(.1)
  

  
%% Iteration loop


    itc=0;   

    stln={'.-','-','--'} % plot line style for each k_x, increase as required
    
  tic;  
    for cvstep=1:cstep   
        
    itc=itc+1;    
    
    [uLn,vLn,wLn,gLn,uLj,vLj,wLj,gLj,cgr(itc,:),cgrj(itc,:)] = mode_tt_v2_module(uLn,vLn,wLn,gLn,uLj,vLj,wLj,gLj,umean1,vmean1,wmean1,gmean1,nstep,h,Lkn);
    
    if rem(cvstep,5)==0

   % P11 = test_base_convergence(uLjb(:,:,:,:,tb),vLjb(:,:,:,:,tb),wLjb(:,:,:,:,tb),gLjb(:,:,:,:,tb),uLnb(:,:,:,:,tb-1),vLnb(:,:,:,:,tb-1),wLnb(:,:,:,:,tb-1),gLnb(:,:,:,:,tb-1));

   % stop_flags     

        figure(21);clf
        for Lkk=1:length(kLyap)
        Lin = find( Lkn == kLyap(Lkk));
        plot((101-1)*h*dt*[1:itc-1],cgr(1:itc-1,Lin),stln{Lkk})
        hold on
       
        ylim([-0.3 0.3])
        end

        drawnow
    end
    
    end
    toc
    
    %% evolution of selected distrubance

    view_mode_time_v2_k(uLn,vLn,wLn,gLn,umean1,vmean1,wmean1,gmean1,NT,h,Lkn(1:2:end),1,1);
    %%

        for Lkk=1:length(kLyap)
        Lin = find( Lkn == kLyap(Lkk));
        gr(1,:,Lkk)=cgr(end,Lin);
        agr(1,:,Lkk)=cgrj(end,Lin);
        end

    %gr(1,:,:)=cgr(end,[1:2:NLyap],:);  
    
    [uL,vL,wL,gL]=decompact_Lyap_Fn(uLn,vLn,wLn,gLn,Lkn);
    [uJ,vJ,wJ,gJ]=decompact_Lyap_Fn(uLj,vLj,wLj,gLj,Lkn);
 
%     % test orth
% 
% for kcount1=1:kLyap
% for kcount2=1:kLyap
%     for lcount1=1:NLyap
%         for lcount2=1:NLyap
%             
%             [uL0,vL0,wL0,gL0]=pick_vec(uL,vL,wL,gL,lcount1,kcount1);
%             [uL1,vL1,wL1,gL1]=pick_vec(uJ,vJ,wJ,gJ,lcount2,kcount2);
%             P11(lcount1,lcount2,kcount1,kcount2)=project_Lyap(uL0,vL0,wL0,uL1,vL1,wL1,kLyap(kcount1)+1);
%         end
%     end
% end
% end
    
plot_mode(37,uL,vL,wL,gL,1,1,umean1,vmean1,wmean1)

rtime=rtold+h*cvstep*dt*(nstep-1);

write_to_disk_ensemble_opti(gather(uL),gather(vL),gather(wL),gather(gL),rtime,kLyap,NLyap,gather(gr),[field_path,'modes_T',num2str(Tinit(1),fmt),'_k123.mat'])
write_to_disk_ensemble_opti(gather(uJ),gather(vJ),gather(wJ),gather(gJ),rtime,kLyap,NLyap,gather(agr),[field_path,'adjoi_T',num2str(Tinit(1),fmt),'_k123.mat'])

%start_mode=[field_path,'modesa150p0k12.',num2str(T(it),'%04.2f'),'.mat']
% start_adjoint=[field_path,'adjoi.vw.k1.h5.',num2str(tt,'%04.2f'),'.mat']
%     
%          [vi,gi,UP,WP,NL,KL]=read_from_disk_compact_ensemble_Lyap(start_adjoint
%          [u1,v1,w1,g1]=make_uw_ens_Lyap(vi,gi,UP,WP,0*b1,0*b2,NL,KL);
%          
%          plot_mode(35,u0,v0,w0,g0,1,1,umean1,vmean1,wmean1)
% 
%         [vi,gi,UP,WP,NL,KL]=read_from_disk_compact_ensemble_Lyap(start_mode);
%         [u1,v1,w1,g1]=make_uw_ens_Lyap(vi,gi,UP,WP,0*b1,0*b2,NL,KL);
%         
%         plot_mode(36,u1,v1,w1,g1,4,1,umean1,vmean1,wmean1)
