function [ V,u,rho] = quadraticControlAlternationsWithResetNoGThreeStepsAltBounds(x_mss,u_mss,f,V0,T,rho0,a,b,d,constraint,sample_time)
do_backoff = false;
backoff_ratio = .08;

if nargin < 10
  constraint = [];
end

if nargin < 11
  sample_time = false;
end
t_sample = linspace(0,T,5);

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
  if ~sample_time
    prog = prog.withIndeterminate(t);
    sproc_vars = [t;x_mss];
  else
    sproc_vars = x_mss;
  end
  [prog,gamma] = prog.newPos(1);
  prog = prog.withPos(1 - gamma);
  [prog,u] = prog.newFreePoly(monomials(sproc_vars,0:3),length(u_mss));
  Vdot = diff(V,x_mss)*subs(f,u_mss,u) + diff(V,t);
  rhodot = diff(rho,t);
  
  [prog, Vdot_sos,mult] = spotless_add_eq_sprocedure(prog, rhodot-Vdot, rho-V,sproc_vars,4);
  if ~sample_time
    [prog, Vdot_sos] = spotless_add_sprocedure(prog, Vdot_sos, t*(T-t),sproc_vars,4);
  end
  
  if ~isempty(constraint)
    [prog, Vdot_sos] = spotless_add_sprocedure(prog, Vdot_sos, constraint,sproc_vars,[]);
  end
  
  Vdot_degree = even_degree(Vdot_sos,sproc_vars)-2;
  
  if ~sample_time
    prog = prog.withSOS(Vdot_sos - gamma*(x_mss'*x_mss)^(Vdot_degree/2));
  else
    for i=1:length(t_sample),
      prog = prog.withSOS(subs(Vdot_sos - gamma*(x_mss'*x_mss)^(Vdot_degree/2),t,t_sample(i)));
    end
  end
  
  for i=1:length(u),
    [prog, up_sos,upmult{i}] = spotless_add_sprocedure(prog, 1-u(i), rho-V,sproc_vars,2);
    [prog, um_sos,ummult{i}] = spotless_add_sprocedure(prog, u(i)+1, rho-V,sproc_vars,2);
    
    if ~sample_time
      [prog, up_sos] = spotless_add_sprocedure(prog, up_sos, t*(T-t),sproc_vars,2);
      prog = prog.withSOS(up_sos - 0*gamma*(x_mss'*x_mss)^2);
      
      [prog, um_sos] = spotless_add_sprocedure(prog, um_sos, t*(T-t),sproc_vars,2);
      prog = prog.withSOS(um_sos - 0*gamma*(x_mss'*x_mss)^2);
    else
      for j=1:length(t_sample),
        prog = prog.withSOS(subs(up_sos - 0*gamma*(x_mss'*x_mss)^2,t,t_sample(j)));
        prog = prog.withSOS(subs(um_sos - 0*gamma*(x_mss'*x_mss)^2,t,t_sample(j)));
      end
    end
  end
  
  rad_sos = subs(V-rho,t,T)*(1+x_mss'*x_mss);
  [prog, rad_sos] = spotless_add_sprocedure(prog, rad_sos, 4*a*d-b^2,x_mss,2);
%   [prog, rad_sos] = spotless_add_sprocedure(prog, rad_sos, b,x_mss,2);
  prog = prog.withSOS(rad_sos);
  
  ab_sos = subs(V-rho,t,T)*(1+x_mss'*x_mss);
  [prog, ab_sos] = spotless_add_sprocedure(prog, ab_sos, -2*a+b,x_mss,2);
  [prog, ab_sos] = spotless_add_sprocedure(prog, ab_sos, a-b+d,x_mss,2);
%   [prog, ab_sos] = spotless_add_sprocedure(prog, ab_sos, -b,x_mss,2);
  prog = prog.withSOS(ab_sos);  

%   [prog, rad_sos,rad_mult] = spotless_add_sprocedure(prog, (b^2-4*a*d)*(1+x_mss'*x_mss), subs(rho-V,t,T),x_mss,2);
%   [prog, rad_sos] = spotless_add_sprocedure(prog, rad_sos, b,x_mss,2);
%   prog = prog.withSOS(rad_sos);  
  
%   [prog, ab_sos,ab_mult] = spotless_add_sprocedure(prog, (2*a-b)*(1+x_mss'*x_mss), subs(rho-V,t,T),x_mss,2);
%   [prog, ab_sos] = spotless_add_sprocedure(prog, ab_sos, b,x_mss,2);
%   prog = prog.withSOS(ab_sos);  
  
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
%   ab_mult = sol.eval(ab_mult);
%   rad_mult = sol.eval(rad_mult);
%   keyboard
  
  mult = sol.eval(mult);
  
  
  %%
% keyboard
  %% step 2
  prog = spotsosprog;
  prog = prog.withIndeterminate(x_mss);
  if ~sample_time
    prog = prog.withIndeterminate(t);
  end
  
  [prog,u] = prog.newFreePoly(monomials(sproc_vars,0:3),length(u_mss));
  Vdot = diff(V,x_mss)*subs(f,u_mss,u) + diff(V,t);
  
  [prog,rho] = prog.newFreePoly(monomials(t,0:2));
  rhodot = diff(rho,t);
  
  Vdot_sos = rhodot - Vdot - (rho-V)*mult;
  if ~sample_time
    [prog, Vdot_sos] = spotless_add_sprocedure(prog, Vdot_sos, t*(T-t),sproc_vars,4);
  end
  if ~isempty(constraint)
    [prog, Vdot_sos] = spotless_add_sprocedure(prog, Vdot_sos, constraint,sproc_vars,[]);
  end
  
  Vdot_degree = even_degree(Vdot_sos,sproc_vars);
  if ~sample_time
    prog = prog.withSOS(Vdot_sos + 1e-6*(x_mss'*x_mss)^(Vdot_degree/2));
  else
    for i=1:length(t_sample),
        prog = prog.withSOS(subs(Vdot_sos + 1e-6*(x_mss'*x_mss)^(Vdot_degree/2),t,t_sample(i)));
    end
  end
  
  for i=1:length(u),
    up_sos = 1 - u(i) - upmult{i}*(rho-V);
    um_sos = 1 + u(i) - ummult{i}*(rho-V);
    
    if ~sample_time
      [prog, up_sos] = spotless_add_sprocedure(prog, up_sos, t*(T-t),sproc_vars,2);
      [prog, um_sos] = spotless_add_sprocedure(prog, um_sos, t*(T-t),sproc_vars,2);
      prog = prog.withSOS(up_sos);
      prog = prog.withSOS(um_sos);
    else
      for j=1:length(t_sample)
        prog = prog.withSOS(subs(up_sos,t,t_sample(j)));
        prog = prog.withSOS(subs(um_sos,t,t_sample(j)));
      end
    end
  end
  
  rad_sos = subs(V-rho,t,T)*(1+x_mss'*x_mss);
  [prog, rad_sos] = spotless_add_sprocedure(prog, rad_sos, 4*a*d-b^2,x_mss,2);
  prog = prog.withSOS(rad_sos);
  
  ab_sos = subs(V-rho,t,T)*(1+x_mss'*x_mss);
  [prog, ab_sos] = spotless_add_sprocedure(prog, ab_sos, -2*a+b,x_mss,2);
  [prog, ab_sos] = spotless_add_sprocedure(prog, ab_sos, a-b+d,x_mss,2);
  prog = prog.withSOS(ab_sos);
  
%   rad_sos = (b^2-4*a*d)*(1+x_mss'*x_mss) -  subs(rho-V,t,T)*rad_mult;
%   [prog, rad_sos] = spotless_add_sprocedure(prog, rad_sos, b,x_mss,2);
%   prog = prog.withSOS(rad_sos);  
%   
%   ab_sos = (2*a-b)*(1+x_mss'*x_mss) - subs(rho-V,t,T)*ab_mult;
%   [prog, ab_sos] = spotless_add_sprocedure(prog, ab_sos, b,x_mss,2);
%   prog = prog.withSOS(ab_sos);  
  
  spot_options = spotprog.defaultOptions;
  spot_options.verbose = true;
  spot_options.do_fr = false;
  spot_options.sos_slack = 0e-6;
  solver = @spot_mosek;
  cost = -subs(rho,t,0) - subs(rho,t,T);
  sol = prog.minimize(cost,solver,spot_options);
%   keyboard
  if do_backoff
    prog_bkp = prog;
    [prog,slack] = prog.newPos(1);
    display('Backing off in iter 2 and resolving');
    % cost < sol.eval(cost) + ratio*abs(sol.eval(cost))
%     prog = prog.withPos(-cost + double(sol.eval(cost)) + .2*abs(double(sol.eval(cost))));
    sol = prog.minimize(slack,solver,spot_options);
%     keyboard
  end
  
  u = sol.eval(u)
  rho = sol.eval(rho);
  %%
%   keyboard
  
  %% step 3
  prog = spotsosprog;
  prog = prog.withIndeterminate(x_mss);
  if ~sample_time
    prog = prog.withIndeterminate(t);
  end
  
%   [prog,gamma] = prog.newPos(1);

  rho_T_nom = double(subs(rho,t,T));

  [prog,rho] = prog.newFreePoly(monomials(t,1:2));
  rho = rho  + 1; % set rho(0) = 1

  [prog,V] = prog.newFreePoly(reshape(monomials(x_mss,2:2)*monomials(t,(0:2))',[],1));
% V = V0
  S0 = diff(diff(subs(V,t,0),x_mss)',x_mss)/2;
  ST = diff(diff(subs(V,t,T),x_mss)',x_mss)/2;
%   [prog,Q] = prog.newPSD(nX);
%   V = x_mss'*Q*x_mss;
%   S0 = Q;
%   ST = Q;
  
  Vdot = diff(V,x_mss)*subs(f,u_mss,u) + diff(V,t);
  rhodot = diff(rho,t);
    
  Vdot_sos = rhodot - Vdot - (rho-V)*mult;
  if ~sample_time
    [prog, Vdot_sos] = spotless_add_sprocedure(prog, Vdot_sos, t*(T-t),sproc_vars,4);
  end
  if ~isempty(constraint)
    [prog, Vdot_sos] = spotless_add_sprocedure(prog, Vdot_sos, constraint,sproc_vars,[]);
  end  
  
  Vdot_degree = even_degree(Vdot_sos,sproc_vars);
  if ~sample_time
    prog = prog.withSOS(Vdot_sos + 1e-6*(x_mss'*x_mss)^(Vdot_degree/2));
  else
    for i=1:length(t_sample)
      prog = prog.withSOS(subs(Vdot_sos + 1e-6*(x_mss'*x_mss)^(Vdot_degree/2),t,t_sample(i)));
    end
  end
  
  for i=1:length(u),
    up_sos = 1 - u(i) - upmult{i}*(rho-V);
    um_sos = 1 + u(i) - ummult{i}*(rho-V);
    
    if ~sample_time
      [prog, up_sos] = spotless_add_sprocedure(prog, up_sos, t*(T-t),sproc_vars,2);
      [prog, um_sos] = spotless_add_sprocedure(prog, um_sos, t*(T-t),sproc_vars,2);
      prog = prog.withSOS(up_sos);
      prog = prog.withSOS(um_sos);
    else
      for j=1:length(t_sample)
        prog = prog.withSOS(subs(up_sos,t,t_sample(j)));
        prog = prog.withSOS(subs(um_sos,t,t_sample(j)));
      end
    end
  end
  
  rad_sos = subs(V-rho,t,T)*(1+x_mss'*x_mss);
  [prog, rad_sos] = spotless_add_sprocedure(prog, rad_sos, 4*a*d-b^2,x_mss,2);
  prog = prog.withSOS(rad_sos);
  
  ab_sos = subs(V-rho,t,T)*(1+x_mss'*x_mss);
  [prog, ab_sos] = spotless_add_sprocedure(prog, ab_sos, -2*a+b,x_mss,2);
  [prog, ab_sos] = spotless_add_sprocedure(prog, ab_sos, a-b+d,x_mss,2);
  prog = prog.withSOS(ab_sos);
  
%   rad_sos = (b^2-4*a*d)*(1+x_mss'*x_mss) -  subs(rho-V,t,T)*rad_mult;
%   [prog, rad_sos] = spotless_add_sprocedure(prog, rad_sos, b,x_mss,2);
%   prog = prog.withSOS(rad_sos);  
%   
%   ab_sos = (2*a-b)*(1+x_mss'*x_mss) - subs(rho-V,t,T)*ab_mult;
%   [prog, ab_sos] = spotless_add_sprocedure(prog, ab_sos, b,x_mss,2);
%   prog = prog.withSOS(ab_sos);  
  
  spot_options = spotprog.defaultOptions;
  spot_options.verbose = true;
  spot_options.do_fr = false;
  spot_options.sos_slack = 0e-6;
  solver = @spot_mosek;
  
%   cost = 0.01*trace(S0)+5*S0(2,2);

  scale_mat = eye(length(x_mss));
  scale_mat(1) = 2;
  
  det_init = det(scale_mat*Q_init*scale_mat');
  det_init_T = det(scale_mat*Q_init_T*scale_mat');

  % linearization of determinant
  cost_coeffs = det_init*inv(scale_mat*Q_init*scale_mat');
  cost_coeffs_T = det_init_T*inv(scale_mat*Q_init_T*scale_mat');
  
  cost_rho = -length(x_mss)/rho_T_nom*det_init_T*subs(rho,t,T); 
   
  cost = 1*sum(sum(scale_mat*(S0-Q_init)*scale_mat'.*cost_coeffs));
  cost = cost + 10*rho_T_nom^(-length(x_mss))*(sum(sum(scale_mat*(ST-Q_init_T)*scale_mat'.*cost_coeffs_T)) + cost_rho);
  
  cost = cost/norm(cost_coeffs(:),inf);
  
  xstar = [sqrt(10);0;0;-1;zeros(2,1)];
  cost = cost + .1*xstar'*S0*xstar;
  
%   cost = cost + S0(1,1);
  
  sol = prog.minimize(cost,solver,spot_options);
  
  if do_backoff
    display('Backing off in iter 2 and resolving');
    prog = prog.withPos(cost + backoff_ratio*abs(double(sol.eval(cost))));
    sol = prog.minimize(0,solver,spot_options);
  end
  
 
  det_new = det(scale_mat*double(sol.eval(S0))*scale_mat');
  display(sprintf('Determinant from %f to %f, percent change %f',det_init,det_new,100-100*det_new/det_init));

  
  V = sol.eval(V);
  rho = sol.eval(rho);
  %%
%     keyboard
  
end

