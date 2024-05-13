# dns-ensemble-lib

A collection of scripts that perform Direct Numerical Simulations of single or ensemble Couette and Poiseuille flow turbulence 
in plane parallel,streamwise and spanwise periodic 3D domains. The DNS ensemble time-stepping routines have been vectorized 
to run simultaneously, reducing the temporal cost of each ensemble member. Executable scripts are located in the 'DNS_Library_v1/' folder.
Additional diagnostic scripts that perform the optimal and modal stability analysis on streamwise-mean base flows are also included.

An initialization script 'NL3D_Couette_gpu_init_ensemble.m' is provided, which advances multiple copies of a single initial 
state for a short time interval with a stochastic forcing term and stores the resulting states. This may be utilized to 
perform transition experiments from a laminar base flow, using the script 'NL3D_Couette_gpu_init_transition_ensemble.m'. 

The main script 'NL3D_Couette_gpu_ensemble.m' continues simulation of these states for a specified number of ensemble members. 
Stochastic forcing can also be enabled in this script. 

To construct a Noise Structure for the appropriate problem resolution we utilize the function 'Functions_stochastic/ NoiseStructure2_cheb.m'.
The last input parameter of this function exports the structure to a save file ( 1 == save , 0 == off) 

The DNS settings are identical to the single flow case (refer to # dns-matlab-gpu). 

Additional options for the ensemble simulations are listed below

Ensemble options
------------
start_filen | file name of single initial state,
init_ens | initialize ensemble member save folders on field_path, 1 == on, 0 == off,
init_file | start ensemble from single file on field_path, 1 == on, 0 == off,
stoch | enables stochastic forcing, 1 == on, 0 == off,
af | List of nonlinear term modulation for perturbation-perturbation interactions in the perturbation equation, 1 for DNS,

Ensemble parameters
------------
Fnn | List of simulated states from the field_path ensemble folder
Fn | Ensemble members / Depends on Fnn
sdt | Gudunov step for stochastic terms / square root of dt

(When solver matrices are loaded these parameters have to be the same with the ones used to precalculate them in order for the DNS to work correctly!)
ksolv | Number of streamwise harmonics solved after dealiasing  
msolv | Number of spanwise harmonics solved after dealiasing  

Stability scripts
------------
Base_flow_modal_analysis_v1g
Base_flow_optimal_analysis_v1g

Modal analysis of a (U,V,W) base flow with underlying symmetries requires considerable amounts of time to distinguish modes with comparable growth.

Stability options
------------
field_path | set path of DNS save state folder,
Tfopt | Target time for optimization,
ty | file type 'file' extracts mean flow from state or 'tseries' selects from time-series of mean flows,
init_file | if ty = 'file' the file name of the base flow containing state,
file_series | if ty = 'tseries' the file name of a series of mean flows,
cont_old | restart from previously obtained stability disturbances,
start_opt | file name of the previous stability disturbances, 
cstep | iterations of stability loop,
nstep | completed timesteps before normalization

Stability parameters
------------
kLyap | range of streamwise harmonic dependence of perturbations,
NLyap | Number of traced modes for each wavenumber,
h | Time-step acceleration factor
 



