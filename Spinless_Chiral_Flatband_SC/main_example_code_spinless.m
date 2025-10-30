% System: Nx by Ny k-grid with a inverse symmetric diamond shape. Nx=Ny and
% both are even. 
Nx = 4;
Ny = 4;

M = 1; % For M > 1, the finite size effect becomes significant, thus a 
% larger system size is needed for a larger M.
s = 1; % This is the k0 in the manuscrierpt.

plotResults = true; % set this to false if no plot needed

% Calculating the stiffness:
[nu,stiffness] = stiffness_chiralSC(Nx,Ny,M,s,plotResults);