function test_conservation_3step_alt()
addpath('../src')
mech = '3step_alt';
kargs = [11.3; 3.6e-2; 7.27e4; 64512; 16294; 5566; 274];
ic = 0.0012;
tend = 1.0;
maxsize = 2500;
x = 1:maxsize;
sol = solve_mepbm(mech,kargs,ic,tend,maxsize);
total = sum(sol.y(1:end-1,:).*[1 1 x(3:end)]',1);
assert(max(abs(total-ic)) < 1e-8,'3step_alt conservation failed.');

disp('3step_alt conservation test passed.')

% check if assert fails
% figure(); plot(sol.x,total,'b.'); % conservation of monomer check
end