%%% Made available under the folling license:
%%% GNU GENERAL PUBLIC LICENSE
%%%                       Version 3, 29 June 2007
%%%
%%% Author: Derek Handwerk

%% Declare globals and load data %%
global maxsize; global r; global x;
global CH; global S2; global S3; global S4; global S5; global atoms; global diam;
[CH,S2,S3,S4,S5,atoms,diam] = loadTEMdata();
fprintf('TEM data loaded.\n');

%% Parameters %%
maxsize = 2500; % Maximum particle size
x = 1:maxsize;
r = 2.677.*(x.^(.72))./x; % Schmidt and Smirnov function

%% Solve and plot %%
sol = solve('3step_alt'); % Change input to 2step or 3step_alt

function sol = solve(mech)
switch mech
    case '2step'
        Kargs = [0.0067; 2524];
    case '3step_alt'
        % Kargs = [65024; 16046; 5787; 197]; % Published
        Kargs = [65536; 16526; 5635; 274]; % Corrected
    otherwise
        error('Incorrect mechanism name.')
end

sol = solvePBM(mech,Kargs);
plotPrecursor(sol);
plotHistograms(sol);

end

%% Solver function %%
function sol = solvePBM(mech,Kargs)
global maxsize;

tspan = [0, 5.05]; % Start and end time for ODE solver (hours)
options = odeset('Stats','off','RelTol',1e-13,'AbsTol',1e-13);
IC = 0.0012; % Precursor initial condition (Moles)
switch mech
    case '2step'
        nuc_order = 1; % MAKE SURE THIS IS WHAT YOU NEED IT TO BE
        ic = zeros(1,maxsize+1);
        ic(1) = IC;
        ode_time = tic;
        sol = ode15s(@(t,n) fw2step(t,n,Kargs(1),Kargs(2),nuc_order),tspan,ic,options);
    case '3step_alt'
        ic = zeros(1,maxsize+2);
        ic(1) = IC;
        S = 11.3; kf = 3.6e-2; kb = 7.27e4;
        ode_time = tic;
        sol = ode15s(@(t,n) fw3step_alttermol(t,n,S,kf,kb,Kargs(1),Kargs(2),Kargs(3),Kargs(4)),tspan,ic,options);
    otherwise
        error('Incorrect mechanism name.')
end
elapsedTime = toc(ode_time);
fprintf(strcat(mech,' elapsed time is %f seconds.\n'), elapsedTime);
end

%% Helper functions %%
% Load Data
function [CH,S2,S3,S4,S5,atoms,diam] = loadTEMdata()
CH = load('CH.mat');
S2 = load('S2.mat');
S3 = load('S3.mat');
S4 = load('S4.mat');
S5 = load('S5.mat');
CH = CH.CH; S2 = S2.S2; S3 = S3.S3; S4 = S4.S4; S5 = S5.S5;
atoms =  load('atoms.mat');
diam = load('diam.mat');
end

% Convert particle of x atoms into corresponding nanometers size
function diam = atomstodiam(x)
diam = 0.3000805*x.^(1/3);
end

%% Plotting functions %%
function plotPrecursor(sol)
global CH;
t = sol.x;
n = sol.y;
figure()
plot(CH(:,1),CH(:,2)./1375,'o'); hold on
plot(t,n(1,:)+n(2,:)); xlabel('Time (hours)'); ylabel('[A] (M)'); % Need to add 1 and 2 together for when there is a_solv. Need total precursor!
xlim([0 5]);
end

function plotHistograms(sol)
% TEM histogram comparison
global maxsize; global x;
global S2; global S3; global S4; global S5;
TEMtimes = [0.918 1.710 2.336 4.838];
TEMsol = deval(sol,TEMtimes);
TEMsol1 = TEMsol(3:maxsize,1);
TEMsol2 = TEMsol(3:maxsize,2);
TEMsol3 = TEMsol(3:maxsize,3);
TEMsol4 = TEMsol(3:maxsize,4);
figure()
% subplot(2,2,1)
title('0.918 hours');
yyaxis left; histogram(S2,25); xlim([0 4.5]); ylim([0 Inf]); xlabel('Size (nm)'); ylabel('TEM data counts');
yyaxis right; plot(atomstodiam(x(3:maxsize)),TEMsol1,'.'); xlim([0 4.5]); ylim([0 Inf]); ylabel('Simualted Ir particles (M)');
hold on; line([mean(TEMsol1) mean(TEMsol1)], [0 max(TEMsol1)], 'LineWidth', 2, 'Color', 'r'); hold off;
figure()
% subplot(2,2,2)
title('1.170 hours');
yyaxis left; histogram(S3,25); xlim([0 4.5]); ylim([0 Inf]); xlabel('Size (nm)'); ylabel('TEM data counts');
yyaxis right; plot(atomstodiam(x(3:maxsize)),TEMsol2,'.'); xlim([0 4.5]); ylim([0 Inf]); ylabel('Simualted Ir particles (M)');
hold on; line([mean(TEMsol2) mean(TEMsol2)], [0 max(TEMsol2)], 'LineWidth', 2, 'Color', 'r'); hold off;
figure()
% subplot(2,2,3)
title('2.336 hours');
yyaxis left; histogram(S4,25); xlim([0 4.5]); ylim([0 Inf]); xlabel('Size (nm)'); ylabel('TEM data counts');
yyaxis right; plot(atomstodiam(x(3:maxsize)),TEMsol3,'.'); xlim([0 4.5]); ylim([0 Inf]); ylabel('Simualted Ir particles (M)');
hold on; line([mean(TEMsol3) mean(TEMsol3)], [0 max(TEMsol3)], 'LineWidth', 2, 'Color', 'r'); hold off;
figure()
% subplot(2,2,4)
title('4.838 hours');
yyaxis left; histogram(S5,25); xlim([0 4.5]); ylim([0 Inf]); xlabel('Size (nm)'); ylabel('TEM data counts');
yyaxis right; plot(atomstodiam(x(3:maxsize)),TEMsol4,'.'); xlim([0 4.5]); ylim([0 Inf]); ylabel('Simualted Ir particles (M)');
hold on; line([mean(TEMsol4) mean(TEMsol4)], [0 max(TEMsol4)], 'LineWidth', 2, 'Color', 'r'); hold off;
end


%% RHS equations
% FW Two-step Mechanism
% k1 nucleation rate
% k2 growth rate
% k3 nucleation order
function dn = fw2step(t,n,k1,k2,k3)
global maxsize; global r;
dn = zeros(maxsize+1,1);
dn(1) = -k3.*k1.*n(1).^k3 - k2.*n(1).*sum(r(3:maxsize).*n(3:maxsize)'.*(3:maxsize)); % Precursor A
dn(2) = 0;
dn(3) = k1.*n(1).^k3 - k3.*k2.*n(1).*r(3).*n(3); % Nucleus
for j = 4:maxsize % All other particle sizes
    dn(j) = k2.*n(1).*(r(j-1).*n(j-1).*(j-1) - r(j).*n(j).*j);
end
dn(maxsize+1) = k2.*n(1).*(r(maxsize).*n(maxsize).*(maxsize)); % keep track of losses 'off-screen'.
end

% New Three-step Mechanism with alternative nucleation
% S  Solvent
% kf kDiss forward rate
% kb KDiss backward rate
% k1 nucleation rate
% k2 growth rate of small particles
% k3 growth rate of large particles
% N  Cutoff between small and large particles
function dn = fw3step_alttermol(t,n,S,kf,kb,k1,k2,k3,N)
N = floor(N);
global maxsize; global r;
dn = zeros(maxsize+2,1);
dn(1) = -k1.*n(1).*n(2).^2 - k2.*n(1).*sum(r(3:N).*n(3:N)'.*(3:N)) - k3.*n(1).*sum(r(N+1:maxsize).*n(N+1:maxsize)'.*(N+1:maxsize)) - kf.*n(1).*S^2 + kb.*n(2).*n(maxsize+2); % precursor A
dn(2) = kf.*n(1).*S^2 - kb.*n(2).*n(maxsize+2) - 2*k1.*n(1).*n(2).^2; % A_solv
dn(3) = k1.*n(1).*n(2).^2 - 3.*k2.*n(1).*r(3).*n(3); % nucleus
for j = 4:N % Small particles
    dn(j) = k2.*n(1).*(r(j-1).*n(j-1).*(j-1) - r(j).*n(j).*j);
end
for j = N+1 % Borderline particle size (formed via k2, lost via k3)
    dn(j) = k2.*n(1).*r(j-1).*n(j-1).*(j-1) - k3.*n(1).*r(j).*n(j).*j;
end
for j = N+2:maxsize % Large particles
    dn(j) = k3.*n(1).*(r(j-1).*n(j-1).*(j-1) - r(j).*n(j).*j);
end
dn(maxsize+1) = k2.*n(1).*(r(maxsize).*n(maxsize).*(maxsize)); % keep track of losses 'off-screen'.
dn(maxsize+2) = kf.*n(1).*S^2 - kb.*n(2).*n(maxsize+2) + k1.*n(1).*n(2).^2 + k2.*n(1).*sum(r(3:N).*n(3:N)'.*(3:N)) + k3.*n(1).*sum(r(N+1:maxsize).*n(N+1:maxsize)'.*(N+1:maxsize)); % POM
end
