function [ V,u,rho] = quadraticControlAlternationsWithResetNoG(x_mss,u_mss,f,V0,T,rho0,a,b,d)
N = 1;
nX = length(x_mss);
V = V0;
rho = rho0;
  t = msspoly('t',1);
for i=1:1,
  prog = spotsosprog;
  prog = prog.withIndeterminate(x_mss);
  prog = prog.withIndeterminate(t);
  [prog,gamma] = prog.newPos(1);
  prog = prog.withPos(1 - gamma);
  [prog,u] = prog.newFreePoly(monomials([t;x_mss],0:3),length(u_mss));
  Vdot = diff(V,x_mss)*subs(f,u_mss,u) + diff(V,t);
  rhodot = diff(rho,t);
  
  [prog, Vdot_sos,mult] = spotless_add_sprocedure(prog, rhodot-Vdot, rho-V,[t;x_mss],4);
  [prog, Vdot_sos] = spotless_add_sprocedure(prog, Vdot_sos, t*(T-t),[t;x_mss],4);
  prog = prog.withSOS(Vdot_sos - gamma*(x_mss'*x_mss)^3);
  for i=1:length(u),
    [prog, u_sos,upmult{i}] = spotless_add_sprocedure(prog, 1-u(i), rho-V,[t;x_mss],2);
    [prog, u_sos] = spotless_add_sprocedure(prog, u_sos, t*(T-t),[t;x_mss],2);
    prog = prog.withSOS(u_sos - gamma*(x_mss'*x_mss)^2);
    [prog, u_sos,ummult{i}] = spotless_add_sprocedure(prog, u(i)+1, rho-V,[t;x_mss],2);
    [prog, u_sos] = spotless_add_sprocedure(prog, u_sos, t*(T-t),[t;x_mss],2);
    prog = prog.withSOS(u_sos - gamma*(x_mss'*x_mss)^2);
  end
  
  [prog, rad_sos,rad_mult] = spotless_add_sprocedure(prog, (b^2-4*a*d)*(1+x_mss'*x_mss), subs(rho-V,t,T),x_mss,2);
  [prog, rad_sos] = spotless_add_sprocedure(prog, rad_sos, b,x_mss,2);
  prog = prog.withSOS(rad_sos);  
  
  [prog, ab_sos,ab_mult] = spotless_add_sprocedure(prog, (2*a-b)*(1+x_mss'*x_mss), subs(rho-V,t,T),x_mss,2);
  [prog, ab_sos] = spotless_add_sprocedure(prog, ab_sos, b,x_mss,2);
  prog = prog.withSOS(ab_sos);  
  
  spot_options = spotprog.defaultOptions;
  spot_options.verbose = true;
  spot_options.do_fr = true;
  solver = @spot_mosek;
  sol = prog.minimize(-gamma,solver,spot_options);
  u = sol.eval(u)
  
  for i=1:length(u),
    upmult{i} = sol.eval(upmult{i});
    ummult{i} = sol.eval(ummult{i});
  end
  ab_mult = sol.eval(ab_mult);
  rad_mult = sol.eval(rad_mult);
%   keyboard
  
  mult = sol.eval(mult);
  
  prog = spotsosprog;
  prog = prog.withIndeterminate(x_mss);
  prog = prog.withIndeterminate(t);
%   [prog,gamma] = prog.newPos(1);
  
  [prog,rho] = prog.newFreePoly(monomials(t,1:2));
  rho = rho  + 1; % set rho(0) = 1

%   [prog,V] = prog.newFreePoly(reshape(monomials(x_mss,(1:2))*monomials(t,(0:2))',[],1));
%   S0 = diff(diff(subs(V,t,0),x_mss)',x_mss)/2;

  [prog,Q] = prog.newPSD(nX);
  V = x_mss'*Q*x_mss;
  S0 = Q;
  
  Vdot = diff(V,x_mss)*subs(f,u_mss,u) + diff(V,t);
  rhodot = diff(rho,t);
    
  Vdot_sos = rhodot-Vdot + mult*(V-rho) + (x_mss'*x_mss)^3*1e-6;  
  [prog, Vdot_sos] = spotless_add_sprocedure(prog, Vdot_sos, t*(T-t),[t;x_mss],4);
  prog = prog.withSOS(Vdot_sos);
  
  for i=1:length(u),
    u_sos = upmult{i}*(V-rho) + 1 - u(i);
    [prog, Vdot_sos] = spotless_add_sprocedure(prog, u_sos, t*(T-t),[t;x_mss],2);
    prog = prog.withSOS(u_sos);
    
    u_sos = ummult{i}*(V-rho) + 1 + u(i);
    [prog, Vdot_sos] = spotless_add_sprocedure(prog, u_sos, t*(T-t),[t;x_mss],2);
    prog = prog.withSOS(u_sos);
  end
  
  
  [prog, rad_sos] = spotless_add_sprocedure(prog, (b^2-4*a*d)*(1+x_mss'*x_mss) - rad_mult*subs(rho-V,t,T), b,x_mss,2);
  prog = prog.withSOS(rad_sos);    
  
  [prog, ab_sos] = spotless_add_sprocedure(prog, (2*a-b)*(1+x_mss'*x_mss) - ab_mult*subs(rho-V,t,T), b,x_mss,2);
  prog = prog.withSOS(ab_sos);
  
  spot_options = spotprog.defaultOptions;
  spot_options.verbose = true;
  spot_options.do_fr = true;
  solver = @spot_mosek;
  sol = prog.minimize(trace(S0)+5*S0(2,2),solver,spot_options);
  
%   keyboard
  
%   solver = @spot_mosek;
%   prog = prog.withPos(1.03*sol.eval(trace(Q)) - trace(Q));
%   sol = prog.minimize(gamma,solver,spot_options);
  
  V = sol.eval(V);
  rho = sol.eval(rho);
  %   keyboard
  
end
end

