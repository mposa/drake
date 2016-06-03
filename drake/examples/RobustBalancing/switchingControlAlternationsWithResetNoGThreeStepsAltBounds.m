function [ V,u,rho] = switchingControlAlternationsWithResetNoGTwoStepsAltBounds(x_mss,f,g,V0,rho0,B0,ta,b,d)
do_backoff = false;
backoff_ratio = .08;

N = 1;
nX = length(x_mss);
nU = size(g,2);
V = V0;
B = B0;
ndgrid_arg = mat2cell(repmat([-1;1],1,nU),2,ones(1,nU)');
[ugrid{1:nU}] = ndgrid(ndgrid_arg{:});
umat = zeros(2^nU,nU);
for i=1:nU,
  umat(:,i) = ugrid{i}(:);
end

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
  sproc_vars = [t;x_mss];

  [prog,gamma] = prog.newPos(1);
  prog = prog.withPos(1 - gamma);
  
  rhodot = diff(rho,t);  
  
  for j=1:2^nU
    Vdot = diff(V,x_mss)*(f + g*umat(j,:)') + diff(V,t);
    Vdot_degree = even_degree(Vdot_sos,sproc_vars)-2;
    [prog, Vdot_sos,mult{j},coeff] = spotless_add_eq_sprocedure(prog, rhodot-Vdot, 1-V,x_mss,4);
    [prog, Vdot_sos] = spotless_add_sprocedure(prog, Vdot_sos, t*(T-t),sproc_vars,4);
    for k=1:nU,
      [prog, Vdot_sos,bmult{j}{k},coeff] = spotless_add_sprocedure(prog, Vdot_sos, umat(j,k)*B(k),x_mss,4);
    end
    prog = prog.withSOS(Vdot_sos - gamma*(x_mss'*x_mss)^(Vdot_degree/2));
  end  
  
  rad_sos = subs(V-rho,t,T)*(1+x_mss'*x_mss);
  [prog, rad_sos] = spotless_add_sprocedure(prog, rad_sos, 4*a*d-b^2,x_mss,2);
  [prog, rad_sos] = spotless_add_sprocedure(prog, rad_sos, b,x_mss,2);
  prog = prog.withSOS(rad_sos);
  
  ab_sos = subs(V-rho,t,T)*(1+x_mss'*x_mss);
  [prog, ab_sos] = spotless_add_sprocedure(prog, ab_sos, 2*a-b,x_mss,2);
  [prog, ab_sos] = spotless_add_sprocedure(prog, ab_sos, -a+b-d,x_mss,2);
  [prog, ab_sos] = spotless_add_sprocedure(prog, ab_sos, -b,x_mss,2);
  prog = prog.withSOS(ab_sos);  

  spot_options = spotprog.defaultOptions;
  spot_options.verbose = true;
  spot_options.do_fr = false;
  solver = @spot_mosek;
  sol = prog.minimize(-gamma,solver,spot_options);
  
  for j=1:2^nU
    mult{j} = sol.eval(mult{j});
    for k=1:nU,
      bmult{j}{k} = sol.eval(bmult{j}{k});
    end
  end
  
  

  
  %% step 2 of 2
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
  [prog, ab_sos] = spotless_add_sprocedure(prog, ab_sos, 2*a-b,x_mss,2);
  [prog, ab_sos] = spotless_add_sprocedure(prog, ab_sos, -a+b-d,x_mss,2);
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

