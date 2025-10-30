% This is the main example code for showing the comparison between
% bootstrap and ED on the Model I discussed in the manuscript.

% System size: Nx by Ny, where both Nx and Ny are odd numbers
Nx = 3;
Ny = 5; % for 3x5 lattice, it takes 20s for one filling. for 5x5, it takes roughly 400s.

% Particle number: N, an even number (from 2 to 2*Nx*Ny) to ensure the 
% superconducting groud state.
N = 6;

% Flat gauge insertion: A=[Ax,Ay]. For Nx = Ny, there is no isotropy. But if Nx is
% not equal to Ny, we shall see slightly different results along two
% directions. The difference vanishes when Nx,Ny -> infinity.
A = [0.01,0]; % this should be sufficiently small.

% Let us first compare the bootstrap result for one filling with ED
% Bootstrap
[Stiffness_boot,runtime] = bootstrapSC(Nx, Ny, N, A);
fprintf('The stiffness for Nx=%d, Ny=%d, and N=%d from bootstrap is %.12f, with runtime %.4fs\n',Nx,Ny,N,Stiffness_boot,runtime)
% ED
[stiffness_ED,runtime] = ED_gamma(Nx,Ny,N,A);
fprintf('The stiffness for Nx=%d, Ny=%d, and N=%d from ED is %.12f, with runtime %.4fs\n',Nx,Ny,N,stiffness_ED,runtime)


% Filling dependence
N_list = 2:2:Nx*Ny; % due to the particle-hole symmetry, we only need to calculate the stiffness for nu<1/2
nu_list = N_list/(2*Nx*Ny);
runtime_boot_list = zeros(size(nu_list));
stiffness_boot_list = zeros(size(N_list));

for i = 1:length(N_list)
    fprintf('Bootstrapping stiffness for Nx=%d, Ny=%d, and N=%d \t',Nx,Ny,N_list(i))
    [stiffness_boot,run_time] = bootstrapSC(Nx, Ny, N_list(i), A);
    stiffness_boot_list(i) = stiffness_boot;
    runtime_boot_list(i) = run_time;
    fprintf('Runtime: %.4fs\n',run_time)
end

num_ED = 4;
runtime_ED_list = zeros(num_ED,1);
stiffness_ED_list = zeros(num_ED,1);
for i = 1:num_ED
    fprintf('ED stiffness for Nx=%d, Ny=%d, and N=%d \t',Nx,Ny,N_list(i))
    [stiffness_ED,run_time] = ED_gamma(Nx, Ny, N_list(i), A);
    stiffness_ED_list(i) = stiffness_ED;
    runtime_ED_list(i) = run_time;
    fprintf('Runtime: %.4fs\n',run_time)
end

% Get the 2-particle pair mass from ED
nu_2particle = 2/(2*Nx*Ny);
inverse_mass_2particle_ED = ED_gamma(Nx,Ny,2,A)/(nu_2particle*(1-nu_2particle));

x = 0:0.01:1;
% plot the results
figure;
s1 = scatter(nu_list, stiffness_boot_list, 50, 's', 'filled', 'MarkerFaceColor','b', 'DisplayName','Bootstrap');
hold on
s2 = scatter(1-nu_list, stiffness_boot_list, 50, 's', 'filled', 'MarkerFaceColor','b');  % no DisplayName
s2.HandleVisibility = 'off';   % <-- exclude from legend
s3 = scatter(nu_list(1:num_ED), stiffness_ED_list, 80, 's', 'MarkerFaceColor','none','MarkerEdgeColor','r', 'DisplayName','ED');
hold on
p1 = plot(x, inverse_mass_2particle_ED.*x.*(1-x), '--', 'Color','k', 'LineWidth',1.5, 'DisplayName','$m_{pair}^{-1}\nu(1-\nu)$');
box on
legend('Location','south',Interpreter='latex')
xlabel('$\nu$',Interpreter='latex')
ylabel('Stiffness $D/|U|$',Interpreter='latex')