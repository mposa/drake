function [gmult_v,gmult_vdot,gmult_rad,gmult_ab] = initQuadraticControlAlternationWithReset(x_mss,u_mss,f,V0,T,rho0,a,b,d)
% trace_val = trace(S0);
t = msspoly('t',1);
V = V0;
rho = rho0;
GO = V/dmsubs(rho,t,T);


prog = spotsosprog;
prog = prog.withIndeterminate(x_mss);
prog = prog.withIndeterminate(t);

[prog,u] = prog.newFreePoly(monomials([t;x_mss],0:3),length(u_mss));
Vdot = diff(V,x_mss)*subs(f,u_mss,u) + diff(V,t);
rhodot = diff(rho,t);

[prog, rad_sos,gmult_rad] = spotless_add_sprocedure(prog, b^2-4*a*d, subs(1-GO,t,T),x_mss,2);
[prog, rad_sos] = spotless_add_sprocedure(prog, rad_sos, b,x_mss,2);
prog = prog.withSOS(rad_sos);

[prog, ab_sos,gmult_ab] = spotless_add_sprocedure(prog, 2*a-b, subs(rho-V,t,T),x_mss,2);
[prog, ab_sos] = spotless_add_sprocedure(prog, ab_sos, b,x_mss,2);
prog = prog.withSOS(ab_sos);



[prog, Vdot_sos,gmult_vdot,coeff] = spotless_add_sprocedure(prog, rhodot-Vdot, 1-GO,[t;x_mss],4);
[prog, Vdot_sos] = spotless_add_sprocedure(prog, Vdot_sos, t*(T-t),[t;x_mss],4);
prog = prog.withSOS(Vdot_sos);


[prog, V_sos,gmult_v,coeff] = spotless_add_eq_sprocedure(prog, V-rho, 1-GO,[t;x_mss],4);
[prog, V_sos] = spotless_add_sprocedure(prog, V_sos, t*(T-t),[t;x_mss],4);
prog = prog.withSOS(V_sos);

for i=1:length(u),
  [prog, u_sos,upmult{i}] = spotless_add_sprocedure(prog, 1-u(i), rho-V,[t;x_mss],2);
  [prog, u_sos] = spotless_add_sprocedure(prog, u_sos, t*(T-t),[t;x_mss],2);
  prog = prog.withSOS(u_sos);
  [prog, u_sos,ummult{i}] = spotless_add_sprocedure(prog, u(i)+1, rho-V,[t;x_mss],2);
  [prog, u_sos] = spotless_add_sprocedure(prog, u_sos, t*(T-t),[t;x_mss],2);
  prog = prog.withSOS(u_sos);
end

spot_options = spotprog.defaultOptions;
spot_options.verbose = true;
spot_options.do_fr = true;
solver = @spot_mosek;
sol = prog.minimize(0,solver,spot_options);

gmult_v = sol.eval(gmult_v);
gmult_vdot = sol.eval(gmult_vdot);
gmult_rad = sol.eval(gmult_rad);
gmult_ab = sol.eval(gmult_ab);
end


