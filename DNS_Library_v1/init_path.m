init_structure = 1; %% Generates the noise basis for a 61x72x72 grid, set to 0 after 1st run
addpath(genpath('../'))

if init_structure
NoiseStructure2_cheb
end
