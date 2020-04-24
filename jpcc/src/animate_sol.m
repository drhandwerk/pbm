function animate_sol(sol,maxsize)
x = 1:maxsize+1;
t = sol.x;
n = sol.y;
scrn_sz = get(groot,'ScreenSize');
figure('Position',[100 100 scrn_sz(3)/3 scrn_sz(4)/2]);
for i = 1:length(t)
    subplot(1,2,1)
    plot(x(3:maxsize),n(3:maxsize,i),'b.'); ylim([0 1.5*max(n(3:maxsize,end))]); xlim([0 maxsize]); xlabel('number of monomers'); ylabel('Moles of particles');
    subplot(1,2,2)
    plot(atoms_to_diam(x(3:maxsize)),n(3:maxsize,i),'b.'); ylim([0 2*max(n(3:maxsize,end))]); xlim([0 4.5]); xlabel('nm'); ylabel('Moles of particles');
    pause(.01);
end
end