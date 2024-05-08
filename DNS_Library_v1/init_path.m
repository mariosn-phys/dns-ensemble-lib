init_structure = 1; %% Generates the noise basis for a 61x72x72 grid with provided parameters, set to 0 after 1st run
                    
% % add current folder and subfolders 
%addpath(genpath('../'))
path_list=["Functions/";"Functions_cheb/";"Functions_ensemble/";"Functions_metrics/";"Functions_modal/";"Functions_stochastic/";"Functions_plot/"];
for ij=1:length(path_list)
addpath(path_list(ij))
end

if init_structure
R=1000;a=2/2;b=2/1;
N=61;NX=72;MZ=72;
NoiseStructure2_cheb(R,N,NX,MZ,a,b,1);
end
