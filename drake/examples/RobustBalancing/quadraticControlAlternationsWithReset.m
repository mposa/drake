function [ V,u,rho,gmult_v,gmult_vdot,gmult_rad,gmult_ab] = quadraticControlAlternationsWithReset(x_mss,u_mss,f,V0,T,rho0,a,b,d,gmult_v,gmult_vdot,gmult_rad,gmult_ab)
% trace_val = trace(S0);

%TODO: try using GO for umin and umax too

N = 1;
nX = length(x_mss);
V = V0;
rho = rho0;
  t = msspoly('t',1);

  prog = spotsosprog;
  prog = prog.withIndeterminate(x_mss);
  [prog,t] = prog.newIndeterminate('t',1);
  [prog,gamma] = prog.newSOSPoly(monomials(t,0:3));
  [prog,u] = prog.newFreePoly(monomials([t;x_mss],0:3),length(u_mss));
  Vdot = diff(V,x_mss)*subs(f,u_mss,u) + diff(V,t);
  rhodot = diff(rho,t);
  
  [prog,GO] = prog.newFreePoly(reshape(monomials(x_mss,(1:2))*monomials(t,(0:2))',[],1));
  
  rad_sos = b^2-4*a*d - gmult_rad*(1-GO);
%   [prog, rad_sos,rad_mult] = spotless_add_sprocedure(prog, b^2-4*a*d, subs(rho-V,t,T),x_mss,2);
  [prog, rad_sos] = spotless_add_sprocedure(prog, rad_sos, b,x_mss,2);
  prog = prog.withSOS(rad_sos);
  
  ab_sos = 2*a-b - gmult_ab*(1-GO);
%   [prog, ab_sos,ab_mult] = spotless_add_sprocedure(prog, 2*a-b, subs(rho-V,t,T),x_mss,2);
  [prog, ab_sos] = spotless_add_sprocedure(prog, ab_sos, b,x_mss,2);
  prog = prog.withSOS(ab_sos);
  
  
  Vdot_sos = rhodot-Vdot-gmult_vdot*(1-GO) + ([t;x_mss]'*[t;x_mss])^3*0e-6;
%   [prog, Vdot_sos,mult,coeff] = spotless_add_eq_sprocedure(prog, rhodot-Vdot, rho-V,[t;x_mss],4);
  [prog, Vdot_sos] = spotless_add_sprocedure(prog, Vdot_sos, t*(T-t),[t;x_mss],4);
  prog = prog.withSOS(Vdot_sos);
  
  V_sos = V-rho-gamma-gmult_v*(1-GO) + ([t;x_mss]'*[t;x_mss])^3*0e-6;
%   [prog, V_sos,gmult_v,coeff] = spotless_add_eq_sprocedure(prog, V-rho, 1-GO,[t;x_mss],4);
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
  
  
  cost = spotlessIntegral(prog,-gamma,[],[],t,[0 T]);
  
  spot_options = spotprog.defaultOptions;
  spot_options.verbose = true;
  spot_options.do_fr = true;
  solver = @spot_mosek;
  sol = prog.minimize(cost,solver,spot_options);
  u = sol.eval(u)
  GO = sol.eval(GO);
  
  for i=1:length(u),
    upmult{i} = sol.eval(upmult{i});
    ummult{i} = sol.eval(ummult{i});
  end
  %   keyboard
  
  
  prog = spotsosprog;
  [prog,t] = prog.newIndeterminate('t',1);
  [prog,rho] = prog.newFreePoly(monomials(t,1:2));
  rho = rho  + 1; % set rho(-) = 1
  
  prog = prog.withIndeterminate(x_mss);
%   [prog,gamma] = prog.newPos(1);
  %   [prog,S] = prog.newPSD(nX);
  
  [prog,V] = prog.newFreePoly(reshape(monomials(x_mss,(1:2))*monomials(t,(0:2))',[],1));
  S = diff(diff(subs(V,t,T),x_mss)',x_mss)/2;
  S0 = diff(diff(subs(V,t,0),x_mss)',x_mss)/2;
  %   V = x_mss'*S*x_mss;
  %   prog = prog.withEqs(-trace(S) + trace_val);
  
  
  [prog, rad_sos,gmult_rad] = spotless_add_sprocedure(prog, b^2-4*a*d, 1-GO,x_mss,2);
  [prog, rad_sos] = spotless_add_sprocedure(prog, rad_sos, b,x_mss,2);
  prog = prog.withSOS(rad_sos);
  
  [prog, ab_sos,gmult_ab] = spotless_add_sprocedure(prog, 2*a-b, 1-GO,x_mss,2);
  [prog, ab_sos] = spotless_add_sprocedure(prog, ab_sos, b,x_mss,2);
  prog = prog.withSOS(ab_sos);
  
%     prog = prog.withPSD(S - Qp*subs(rho,t,T));
  
  Vdot = diff(V,x_mss)*subs(f,u_mss,u) + diff(V,t);
  rhodot = diff(rho,t);
  
  Vdot_sos = rhodot-Vdot + ([t;x_mss]'*[t;x_mss])^3*0e-6;
  [prog, Vdot_sos,gmult_vdot] = spotless_add_sprocedure(prog, Vdot_sos, 1-GO,[t;x_mss],4);
  [prog, Vdot_sos] = spotless_add_sprocedure(prog, Vdot_sos, t*(T-t),[t;x_mss],4);
%   [prog, Vdot_sos] = spotless_add_sprocedure(prog, Vdot_sos, x_mss'*Qp*x_mss - r,[t;x_mss],4);
  prog = prog.withSOS(Vdot_sos);
  
  [prog, V_sos,gmult_v] = spotless_add_eq_sprocedure(prog, V-rho, 1-GO,[t;x_mss],4);
  [prog, V_sos] = spotless_add_sprocedure(prog, V_sos, t*(T-t),[t;x_mss],4);
  prog = prog.withSOS(V_sos);
  
  
  for i=1:length(u),
    u_sosp = upmult{i}*(V-rho) + 1 - u(i);
    [prog, u_sosp] = spotless_add_sprocedure(prog, u_sosp, t*(T-t),[t;x_mss],4);
    prog = prog.withSOS(u_sosp);
    u_sosm = ummult{i}*(V-rho) + 1 + u(i);
    [prog, u_sosm] = spotless_add_sprocedure(prog, u_sosm, t*(T-t),[t;x_mss],4);
    prog = prog.withSOS(u_sosm);
  end
  
  spot_options = spotprog.defaultOptions;
  spot_options.verbose = true;
  spot_options.do_fr = true;
  solver = @spot_mosek;
  sol = prog.minimize(trace(S0),solver,spot_options);
  
%   keyboard
  
%   solver = @spot_mosek;
%   prog = prog.withPos(1.03*sol.eval(trace(Q)) - trace(Q));
%   sol = prog.minimize(gamma,solver,spot_options);
  V = sol.eval(V);
  rho = sol.eval(rho)
  gmult_v = sol.eval(gmult_v);
  gmult_vdot = sol.eval(gmult_vdot);
  gmult_ab = sol.eval(gmult_ab);
  gmult_rad = sol.eval(gmult_rad);
%   S = sol.eval(S);
  %   keyboard
  
end

