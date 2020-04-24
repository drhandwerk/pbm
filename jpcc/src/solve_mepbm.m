function sol = solve_mepbm(mech,kargs,ic_in,tend,maxsize)
x = 1:maxsize+1;
tspan = [0 tend];
ic = zeros(1,maxsize+1);
ic(1) = ic_in;
options = odeset('Stats','off','RelTol',1e-13,'AbsTol',1e-13);
r = 2.677.*(x.^(.72))./x;
sol = ode15s(@(t,n) mepbm_rhs(t,n,mech,kargs,maxsize,r),tspan,ic,options);
log_mepbm(mech,kargs,ic_in,tend,maxsize,sol)
end