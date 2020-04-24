function plot_precursor(sol,pcdata) 
t = sol.x;
n = sol.y;
figure;
plot(pcdata(:,1),pcdata(:,2)./1375,'o'); hold on
plot(t,n(1,:)+n(2,:)); xlabel('Time (hours)'); ylabel('[A] (M)'); % Need to add 1 and 2 together for when there is a_solv. Need total precursor!
xlim([0 5]);
end