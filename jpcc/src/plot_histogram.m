function plot_histogram(sol,histdata,maxsize,time)
x = 1:maxsize;
tem_sol = deval(sol,time);
figure;
title([num2str(time) ' hours']);
yyaxis left; histogram(histdata,25); xlim([0 4.5]); ylim([0 Inf]); xlabel('Size (nm)'); ylabel('TEM data counts');
yyaxis right; plot(atoms_to_diam(x(3:maxsize)),tem_sol(3:maxsize),'.'); xlim([0 4.5]); ylim([0 Inf]); ylabel('Simualted Ir particles (M)');
end