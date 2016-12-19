clear all
x = msspoly('x',1);
t = msspoly('t',1);
u_mss = msspoly('u',1);
V = 5*x^2;
f = x + u_mss;
rho = 1 + 5*t;

T = 1;
xT = 5;

GO = V/dmsubs(rho,t,T);


%%
prog = spotsosprog();
prog = prog.withIndeterminate(x);
prog = prog.withIndeterminate(t);
[prog,u] = prog.newFreePoly(monomials([t;x],0:3));

Vdot_sos = diff(rho,t) - diff(V,x)*subs(f,u_mss,u) - diff(V,t);
[prog,Vdot_sos,gmult1] = spotless_add_sprocedure(prog, Vdot_sos, 1-GO,[t;x],2);
[prog,Vdot_sos] = spotless_add_sprocedure(prog, Vdot_sos, t*(T-t),[t;x],4);
prog = prog.withSOS(Vdot_sos);

[prog,V_sos,gmult2] = spotless_add_eq_sprocedure(prog, V-rho, 1-GO,[t;x],2);
[prog,V_sos] = spotless_add_sprocedure(prog, V_sos, t*(T-t),[t;x],4);
prog = prog.withSOS(V_sos);

for i=1:length(u),
  [prog, u_sos,upmult{i}] = spotless_add_sprocedure(prog, 1-u(i), rho-V,[t;x],0);
  [prog, u_sos] = spotless_add_sprocedure(prog, u_sos, t*(T-t),[t;x],2);
  prog = prog.withSOS(u_sos);
  [prog, u_sos,ummult{i}] = spotless_add_sprocedure(prog, u(i)+1, rho-V,[t;x],0);
  [prog, u_sos] = spotless_add_sprocedure(prog, u_sos, t*(T-t),[t;x],2);
  prog = prog.withSOS(u_sos);
end

spot_options = spotprog.defaultOptions;
spot_options.verbose = true;
spot_options.do_fr = true;
solver = @spot_mosek;
sol = prog.minimize(0,solver,spot_options);

gmult1 = sol.eval(gmult1);
gmult2 = sol.eval(gmult2);

%%
for i=1:15,

prog = spotsosprog();
prog = prog.withIndeterminate(x);
prog = prog.withIndeterminate(t);
[prog,u] = prog.newFreePoly(monomials([t;x],0:3));
[prog,gamma] = prog.newPos(1);

[prog,GO] = prog.newFreePoly(reshape(monomials(x,(1:2))*monomials(t,(0:2))',[],1));

Vdot_sos = diff(rho,t) - diff(V,x)*subs(f,u_mss,u) - diff(V,t) - gmult1*(1-GO);
[prog,Vdot_sos,mult2] = spotless_add_sprocedure(prog, Vdot_sos, t*(T-t),[t;x],4);
prog = prog.withSOS(Vdot_sos);

V_sos = V-rho-gamma-gmult2*(1-GO);
[prog,V_sos,mult2] = spotless_add_sprocedure(prog, V_sos, t*(T-t),[t;x],4);
prog = prog.withSOS(V_sos);

for i=1:length(u),
  [prog, u_sos,upmult{i}] = spotless_add_sprocedure(prog, 1-u(i), rho-V,[t;x],0);
  [prog, u_sos] = spotless_add_sprocedure(prog, u_sos, t*(T-t),[t;x],2);
  prog = prog.withSOS(u_sos);
  [prog, u_sos,ummult{i}] = spotless_add_sprocedure(prog, u(i)+1, rho-V,[t;x],0);
  [prog, u_sos] = spotless_add_sprocedure(prog, u_sos, t*(T-t),[t;x],2);
  prog = prog.withSOS(u_sos);
end

spot_options = spotprog.defaultOptions;
spot_options.verbose = true;
spot_options.do_fr = true;
solver = @spot_mosek;
sol = prog.minimize(-gamma,solver,spot_options);

% mult = sol.eval(mult);
for i=1:length(u),
  upmult{i} = sol.eval(upmult{i});
  ummult{i} = sol.eval(ummult{i});
end

u = sol.eval(u);
GO = sol.eval(GO);

prog = spotsosprog();
prog = prog.withIndeterminate(x);
prog = prog.withIndeterminate(t);

% [prog,V] = prog.newFreePoly(monomials(x,1:2));
[prog,V] = prog.newFreePoly(reshape(monomials(x,(1:2))*monomials(t,(0:2))',[],1));
% [prog,V] = prog.newFreePoly(monomials([t;x],0:4));
S0 = diff(diff(subs(V,t,0),x)',x)/2;
S1 = diff(diff(subs(V,t,T),x)',x)/2;
[prog,rho] = prog.newFreePoly(monomials(t,1:1));
rho = rho + 1;

[prog,z] = prog.withPos(S1*xT^2-subs(rho,t,T));

  


Vdot_sos = diff(rho,t) - diff(V,x)*subs(f,u_mss,u) - diff(V,t);
[prog,Vdot_sos,gmult1] = spotless_add_sprocedure(prog, Vdot_sos, 1-GO,[t;x],2);
[prog,Vdot_sos] = spotless_add_sprocedure(prog, Vdot_sos, t*(T-t),[t;x],4);
prog = prog.withSOS(Vdot_sos);

[prog,V_sos,gmult2] = spotless_add_eq_sprocedure(prog, V-rho, 1-GO,[t;x],2);
[prog,V_sos] = spotless_add_sprocedure(prog, V_sos, t*(T-t),[t;x],4);
prog = prog.withSOS(V_sos);


for i=1:length(u),
  u_sos = upmult{i}*(V-rho) + 1 - u(i);
  [prog, u_sos] = spotless_add_sprocedure(prog, u_sos, t*(T-t),[t;x],2);
  prog = prog.withSOS(u_sos);
  u_sos = ummult{i}*(V-rho) + 1 + u(i);
  [prog, u_sos] = spotless_add_sprocedure(prog, u_sos, t*(T-t),[t;x],2);
  prog = prog.withSOS(u_sos);
end

spot_options = spotprog.defaultOptions;
spot_options.verbose = true;
spot_options.do_fr = true;
solver = @spot_mosek;
sol = prog.minimize(S0,solver,spot_options);

V = sol.eval(V);
rho = sol.eval(rho);
gmult1 = sol.eval(gmult1);
gmult2 = sol.eval(gmult2);

end
t = msspoly('t',1);