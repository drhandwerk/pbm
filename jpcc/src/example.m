%% Example solve
mech = '3step_alt';
kargs = [11.3; 3.6e-2; 7.27e4; 64512; 16294; 5566; 274];
ic = 0.0012;
tend = 5.05;
maxsize = 2500;
sol = solve_mepbm(mech,kargs,ic,tend,maxsize);
TEMtimes = [0.918 1.710 2.336 4.838];
plot_histogram(sol,S2,maxsize,TEMtimes(1));
plot_histogram(sol,S3,maxsize,TEMtimes(2));
plot_histogram(sol,S4,maxsize,TEMtimes(3));
plot_histogram(sol,S5,maxsize,TEMtimes(4));
plot_precursor(sol,CH);

%% Example fit
% Initial values:
fit_type = 'histogram';                 % 'precursor' or 'histogram'
mech = '3step_alt';                     % Mechanism name (ODE RHS)
kconst = [10.1; 0.007];                 % Fixed parameters
kvar = [1000; 49024; 13496; 393; 149];  % Variable parameters
lb = [9090.9; 4800; 10; 10; 10];        % Lower kvar bounds
ub = [14286; 8e+07; 85000; 25000; 800]; % Upper kvar bounds
ic = 0.0012;                            % Initial precursor
tend = 7.9;                             % Final time
maxsize = 2000;                         % Maximum Particle size

% Run the fit with "dust_filt" data.
[k,fval,exitflag,output] = fit_mepbm(mech,kconst,kvar,ic,tend,maxsize,'histogram',dust_filt,lb,ub);