function [ V,u,rho] = quadraticControlAlternationsWithResetNoGThreeSteps(x_mss,u_mss,f,V0,T,rho0,a,b,d,constraint)
if nargin < 10
  constraint = [];
end

N = 1;
nX = length(x_mss);
V = V0;

rho = rho0;
  t = msspoly('t',1);
Q_init = double(diff(diff(subs(V,t,0),x_mss)',x_mss))/2;
Q_init_T = double(diff(diff(subs(V,t,T),x_mss)',x_mss))/2;
for i=1:1,
  %% step 1
  prog = spotsosprog;
  prog = prog.withIndeterminate(x_mss);
  prog = prog.withIndeterminate(t);
  [prog,gamma] = prog.newFree(1);
  prog = prog.withPos(1 - gamma);
  [prog,u] = prog.newFreePoly(monomials([t;x_mss],0:3),length(u_mss));
  Vdot = diff(V,x_mss)*subs(f,u_mss,u) + diff(V,t);
  rhodot = diff(rho,t);
  
  [prog, Vdot_sos,mult] = spotless_add_eq_sprocedure(prog, rhodot-Vdot, rho-V,[t;x_mss],4);
  [prog, Vdot_sos] = spotless_add_sprocedure(prog, Vdot_sos, t*(T-t),[t;x_mss],4);
  
  if ~isempty(constraint)
    [prog, Vdot_sos] = spotless_add_sprocedure(prog, Vdot_sos, constraint,[t;x_mss],[]);
  end
  
  Vdot_degree = even_degree(Vdot_sos,[t;x_mss]);
  
  prog = prog.withSOS(Vdot_sos - gamma*(x_mss'*x_mss)^(Vdot_degree/2));
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
  spot_options.do_fr = false;
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
  
  
  

  %% step 2
  prog = spotsosprog;
  prog = prog.withIndeterminate(x_mss);
  prog = prog.withIndeterminate(t);
  [prog,u] = prog.newFreePoly(monomials([t;x_mss],0:3),length(u_mss));
  Vdot = diff(V,x_mss)*subs(f,u_mss,u) + diff(V,t);
  
  [prog,rho] = prog.newFreePoly(monomials(t,0:2));
  rhodot = diff(rho,t);
  
  Vdot_sos = rhodot - Vdot - (rho-V)*mult;
  [prog, Vdot_sos] = spotless_add_sprocedure(prog, Vdot_sos, t*(T-t),[t;x_mss],4);
    if ~isempty(constraint)
    [prog, Vdot_sos] = spotless_add_sprocedure(prog, Vdot_sos, constraint,[t;x_mss],[]);
    end  
        
  Vdot_degree = even_degree(Vdot_sos,[t;x_mss]);
  prog = prog.withSOS(Vdot_sos + 1e-6*(x_mss'*x_mss)^(Vdot_degree/2));
  
  for i=1:length(u),
    u_sos = 1 - u(i) - upmult{i}*(rho-V);
    [prog, u_sos] = spotless_add_sprocedure(prog, u_sos, t*(T-t),[t;x_mss],2);
    prog = prog.withSOS(u_sos);
    
    u_sos = 1 + u(i) - ummult{i}*(rho-V);
    [prog, u_sos] = spotless_add_sprocedure(prog, u_sos, t*(T-t),[t;x_mss],2);
    prog = prog.withSOS(u_sos);
  end
  
  rad_sos = (b^2-4*a*d)*(1+x_mss'*x_mss) -  subs(rho-V,t,T)*rad_mult;
  [prog, rad_sos] = spotless_add_sprocedure(prog, rad_sos, b,x_mss,2);
  prog = prog.withSOS(rad_sos);  
  
  ab_sos = (2*a-b)*(1+x_mss'*x_mss) - subs(rho-V,t,T)*ab_mult;
  [prog, ab_sos] = spotless_add_sprocedure(prog, ab_sos, b,x_mss,2);
  prog = prog.withSOS(ab_sos);  
  
  spot_options = spotprog.defaultOptions;
  spot_options.verbose = true;
  spot_options.do_fr = false;
  solver = @spot_mosek;
  sol = prog.minimize(-subs(rho,t,0),solver,spot_options);
  u = sol.eval(u)
  rho = sol.eval(rho);
  
  %% step 3
  prog = spotsosprog;
  prog = prog.withIndeterminate(x_mss);
  prog = prog.withIndeterminate(t);
%   [prog,gamma] = prog.newPos(1);

  rho_T_nom = double(subs(rho,t,T));

  [prog,rho] = prog.newFreePoly(monomials(t,1:2));
  rho = rho  + 1; % set rho(0) = 1

  [prog,V] = prog.newFreePoly(reshape(monomials(x_mss,2:2)*monomials(t,(0:2))',[],1));
  S0 = diff(diff(subs(V,t,0),x_mss)',x_mss)/2;
  ST = diff(diff(subs(V,t,T),x_mss)',x_mss)/2;
%   [prog,Q] = prog.newPSD(nX);
%   V = x_mss'*Q*x_mss;
%   S0 = Q;
%   ST = Q;
  
  Vdot = diff(V,x_mss)*subs(f,u_mss,u) + diff(V,t);
  rhodot = diff(rho,t);
    
  Vdot_sos = rhodot - Vdot - (rho-V)*mult;
  [prog, Vdot_sos] = spotless_add_sprocedure(prog, Vdot_sos, t*(T-t),[t;x_mss],4);
  if ~isempty(constraint)
    [prog, Vdot_sos] = spotless_add_sprocedure(prog, Vdot_sos, constraint,[t;x_mss],[]);
  end  
  
  Vdot_degree = even_degree(Vdot_sos,[t;x_mss]);
  prog = prog.withSOS(Vdot_sos + 1e-6*(x_mss'*x_mss)^(Vdot_degree/2));
  
  for i=1:length(u),
    u_sos = 1 - u(i) - upmult{i}*(rho-V);
    [prog, u_sos] = spotless_add_sprocedure(prog, u_sos, t*(T-t),[t;x_mss],2);
    prog = prog.withSOS(u_sos);
    
    u_sos = 1 + u(i) - ummult{i}*(rho-V);
    [prog, u_sos] = spotless_add_sprocedure(prog, u_sos, t*(T-t),[t;x_mss],2);
    prog = prog.withSOS(u_sos);
  end
  
  
  rad_sos = (b^2-4*a*d)*(1+x_mss'*x_mss) -  subs(rho-V,t,T)*rad_mult;
  [prog, rad_sos] = spotless_add_sprocedure(prog, rad_sos, b,x_mss,2);
  prog = prog.withSOS(rad_sos);  
  
  ab_sos = (2*a-b)*(1+x_mss'*x_mss) - subs(rho-V,t,T)*ab_mult;
  [prog, ab_sos] = spotless_add_sprocedure(prog, ab_sos, b,x_mss,2);
  prog = prog.withSOS(ab_sos);  
  
  spot_options = spotprog.defaultOptions;
  spot_options.verbose = true;
  spot_options.do_fr = false;
  solver = @spot_mosek;
  
%   cost = 0.01*trace(S0)+5*S0(2,2);

  scale_mat = eye(length(x_mss));

  det_init = det(scale_mat*Q_init*scale_mat');
  det_init_T = det(scale_mat*Q_init_T*scale_mat');

  % linearization of determinant
  cost_coeffs = det_init*inv(scale_mat*Q_init*scale_mat');
  cost_coeffs_T = det_init*inv(scale_mat*Q_init_T*scale_mat');
  
  cost_rho = -length(x_mss)/rho_T_nom*det_init_T*subs(rho,t,T); 
   
  cost = 1*sum(sum(scale_mat*(S0-Q_init)*scale_mat'.*cost_coeffs));
  cost = cost + .01*(sum(sum(scale_mat*(ST-Q_init_T)*scale_mat'.*cost_coeffs_T)) + cost_rho);
  
  cost = cost/norm(cost_coeffs(:),inf);
  
  sol = prog.minimize(cost,solver,spot_options);
  
 
  det_new = det(scale_mat*sol.eval(S0)*scale_mat');
  display(sprintf('Determinant from %f to %f, percent change %f',det_init,det_new,100-100*det_new/det_init));
  
%   keyboard

  
  V = sol.eval(V);
  rho = sol.eval(rho);
  %   keyboard
  
end

