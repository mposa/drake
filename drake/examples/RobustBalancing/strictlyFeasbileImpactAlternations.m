function [ V,B,rho] = strictlyFeasbileImpactAlternations(x,s,f,g,V0,rho0,B0,T,r,reset_constraint,Vprev)
nX = length(x);
nU = size(g,2);
ndgrid_arg = mat2cell(repmat([-1;1],1,nU),2,ones(1,nU)');
[ugrid{1:nU}] = ndgrid(ndgrid_arg{:});
umat = zeros(2^nU,nU);
for i=1:nU,
  umat(:,i) = ugrid{i}(:);
end

V = V0;
t = msspoly('t',1);


Q_init = double(diff(diff(subs(V,t,0),x)',x))/2;
Q_init_T = double(diff(diff(subs(V,t,T),x)',x))/2;


% scale problem data based on t=0
% y = T*x
% x = T^-1 y
y = msspoly('y',nX);
scale_transform = Q_init^(.5);
Q_init_y = eye(nX);

% scale_transform = eye(nX);
% Q_init_y = Q_init;


x_y = inv(scale_transform)*y;
f_y = scale_transform*subs(f,x,x_y);
g_y = scale_transform*subs(g,x,x_y);
r_y = scale_transform*subs(r,x,x_y);
reset_constraint_y = subs(reset_constraint,x,x_y);
V0_y = subs(V0,x,x_y);
B0_y = subs(B0,x,x_y);

Q_init_T_y = inv(scale_transform)'*Q_init_T*inv(scale_transform);
Vprev_y = subs(Vprev,x,x_y);


[mult,bmult,rho_y] = binarySearchRho(t,y,f_y,g_y,T,umat,V0_y,B0_y,rho0);
[V0_y,B0_y,rho] = binarySearchVandB(t,y,s,f_y,g_y,T,r_y,reset_constraint_y,Vprev_y,umat,V0_y,Q_init_y,Q_init_T_y,scale_transform,mult,bmult,rho_y);
V = subs(V0_y,y,scale_transform*x);
B = subs(B0_y,y,scale_transform*x);

% [mult,bmult] = iter1(t,x,f,g,T,umat,V0,B0,rho0);
% [V,B,rho] = iter2(t,x,s,f,g,T,r,reset_constraint,Vprev,umat,Q_init,Q_init_T,mult,bmult,rho0);


end


function [mult_opt,bmult_opt,rho_opt] = binarySearchRho(t,x,f,g,T,umat,V,B,rho)
% maxing rho_mult, a scalar multiplier on rho(t)
% s.t. V=rho ==> Vdot < rhodot 
%
% what about jump equation? I feel like I need that in here too
outer_radius = 2;

max_iter = 2;
rho_mult = 1;
rho_mult_min = rho_mult;
rho_mult_max = inf;

nX = length(x);
nU = size(g,2);

for i=1:max_iter
  
  prog = spotsosprog;
  prog = prog.withIndeterminate(x);
  prog = prog.withIndeterminate(t);
  [prog,gamma] = prog.newFree(1);
  prog = prog.withPos(gamma+10);
  sproc_vars = [t;x];
  
  rhodot = diff(rho,t);
  
  for j=1:2^nU
    Vdot = diff(V,x)*(f + g*umat(j,:)') + diff(V,t);
    
    [prog, Vdot_sos,mult{j},coeff] = spotless_add_eq_sprocedure(prog, rho_mult*rhodot-Vdot, rho_mult*rho-V,sproc_vars,2);
    [prog, Vdot_sos] = spotless_add_sprocedure(prog, Vdot_sos, t*(T-t),sproc_vars,4);
%     [prog, Vdot_sos] = spotless_add_sprocedure(prog, Vdot_sos, outer_radius-x'*x,x,4);
    for k=1:nU,
      [prog, Vdot_sos,bmult{j}{k},coeff] = spotless_add_sprocedure(prog, Vdot_sos, umat(j,k)*B(k),sproc_vars,2);
    end
    prog = prog.withSOS(Vdot_sos + gamma);
  end
  
  spot_options = spotprog.defaultOptions;
  spot_options.verbose = true;
  spot_options.do_fr = false;
  spot_options.sos_slack = -1e-8;
  solver = @spot_mosek;
  sol = prog.minimize(gamma,solver,spot_options);
  
  if sol.status ~= spotsolstatus.STATUS_PRIMAL_AND_DUAL_FEASIBLE
    keyboard
  end
  
  if sol.eval(gamma) < -1e-6 || sol.status ~= spotsolstatus.STATUS_PRIMAL_AND_DUAL_FEASIBLE
    rho_mult_min = rho_mult;
    for j=1:2^nU
      mult_opt{j} = sol.eval(mult{j});
      for k=1:nU,
        bmult_opt{j}{k} = sol.eval(bmult{j}{k});
      end
    end
    rho_opt = rho_mult*rho;
  else
    rho_mult_max = rho_mult;
    if i ==1,
      keyboard % uh oh
    end
  end
  
  if ~isinf(rho_mult_max)
    rho_mult = (rho_mult_max + rho_mult_min)/2;
  else
    rho_mult = 1.1*rho_mult;
  end
end

rho_opt = 1; %overwriting rho_opt, just using this search to get better multipliers
end

function [V_opt,B_opt,rho_opt] = binarySearchVandB(t,x,s,f,g,T,r,reset_constraint,Vprev,umat,V_init,Q_init,Q_init_T,scale_transform,mult,bmult,rho_init)
% cost_option
% 1 - determinant
% 2 - integral
cost_option = 2; 
max_iter = 5;
outer_radius = 2;

nX = length(x);
nU = size(g,2);

scale_mat = scale_transform';

for i=1:max_iter,
  %%
  prog = spotsosprog;
  prog = prog.withIndeterminate(x);
  prog = prog.withIndeterminate(t);
  if ~isempty(s)
    prog = prog.withIndeterminate(s);
  end
  sproc_vars = [t;x];
  
  %   [prog,gamma] = prog.newPos(1);
  
  rho_init_T = double(subs(rho_init,t,T));
  
  [prog,rho] = prog.newFreePoly(monomials(t,1:2));
  rho = rho  + 1; % set rho(0) = 1
  
  [prog,V] = prog.newFreePoly(reshape(monomials(x,2:2)*monomials(t,(0:2))',[],1));
  [prog,B] = prog.newFreePoly(reshape(monomials(x,1:2)*monomials(t,(0:2))',[],1),nU);
%   [prog,B] = prog.newFreePoly(reshape(monomials(x,1)*monomials(t,(0:1))',[],1),nU);
%   [prog,B] = prog.newFreePoly(monomials(x,1:2),nU);
  S0 = subs(diff(diff(subs(V,t,0),x)',x)/2,x,zeros(nX,1));
  ST = subs(diff(diff(subs(V,t,T),x)',x)/2,x,zeros(nX,1));
  
  [prog,gamma] = prog.newFree(1);
  prog = prog.withPos(gamma+10);
  
  if i==1,
    % initialize cost
    if cost_option == 1
      [~,cost_const] = calc_cost(nX,t,T,S0,ST,rho,Q_init,Q_init_T,rho_init,scale_mat);
      
      % based on the fact that cost(V0) = 0
      % gets around the fact that subs(Q) doesn't work with Q PSD
      cost_init = -cost_const;
    else
      [~,cost_const] = calc_integral_cost(prog,x,t,V_init,outer_radius,T,scale_mat);
      cost_init = cost_const;
    end
    
    if sign(cost_init) > 0
      cost_mult = 1/1.1;
    else
      cost_mult = 1.1;
    end
    
    cost_min = -inf;
    cost_max = cost_init;
    cost_val = cost_init;
  end
 
  
  prog = prog.withSOS(subs(V,t,0));
  prog = prog.withSOS(subs(V,t,T));
  
  rhodot = diff(rho,t);
  
  for j=1:2^nU
    Vdot = diff(V,x)*(f + g*umat(j,:)') + diff(V,t);
    Vdot_sos = rhodot-Vdot - mult{j}*(rho-V);
    for k=1:nU,
      Vdot_sos = Vdot_sos - bmult{j}{k}*umat(j,k)*B(k);
    end
    [prog, Vdot_sos] = spotless_add_sprocedure(prog, Vdot_sos, t*(T-t),sproc_vars,4);
%     [prog, Vdot_sos] = spotless_add_sprocedure(prog, Vdot_sos, 1000-x'*x,x,4);
    [prog, Vdot_sos] = spotless_add_sprocedure(prog, Vdot_sos, 5*rho_init - V_init,x,2);
    prog = prog.withSOS(Vdot_sos + gamma);
  end
  
  reset_sos = subs(V - rho,t,T)*(1 + [s;x]'*[s;x])^2;
  Vprev_reset = subs(Vprev,x,r);
  %   [prog, reset_sos] = spotless_add_sprocedure(prog, reset_sos,Vprev - 1,[s;x],4);
  [prog, reset_sos] = spotless_add_sprocedure(prog, reset_sos,Vprev_reset - 1,[s;x],4);
  [prog, reset_sos] = spotless_add_sprocedure(prog, reset_sos,2 - Vprev_reset,[s;x],4);
  if ~isempty(reset_constraint)
    [prog, reset_sos] = spotless_add_eq_sprocedure(prog, reset_sos,reset_constraint,[s;x],3);
  end
  prog = prog.withSOS(reset_sos);
  
  if cost_option == 1
    cost = calc_cost(nX,t,T,S0,ST,rho,Q_init,Q_init_T,rho_init,scale_mat);
  else
    cost = calc_integral_cost(prog,x,t,V,outer_radius,T,scale_mat);
  end
  
  
  prog = prog.withPos(cost_val - cost);  
  
  spot_options = spotprog.defaultOptions;
  spot_options.verbose = true;
  spot_options.do_fr = false;
  spot_options.sos_slack = -1e-6;
  solver = @spot_mosek;
  
  sol = prog.minimize(gamma,solver,spot_options);
  %%
  if sol.status ~= spotsolstatus.STATUS_PRIMAL_AND_DUAL_FEASIBLE
    keyboard
  end
  
  if sol.eval(gamma) < -1e-6 || sol.status ~= spotsolstatus.STATUS_PRIMAL_AND_DUAL_FEASIBLE
    cost_max = cost_val;
    V_opt = sol.eval(V);
    B_opt = sol.eval(B);
    rho_opt = sol.eval(rho);
  else
    cost_min = cost_val;
    if i ==1,
      keyboard % uh oh
    end
  end
  
  if ~isinf(cost_min)
    cost_val = (cost_max + cost_min)/2;
  else
    cost_val = cost_val*cost_mult;
  end
end

S0_trans = sol.eval(scale_mat*S0*scale_mat');
ST_trans = sol.eval(scale_mat*ST*scale_mat');
Q_init_trans = scale_mat*Q_init*scale_mat';
Q_init_T_trans = scale_mat*Q_init_T*scale_mat';

det_init = det(Q_init_trans);
det_init_T = det(Q_init_T_trans)*double(subs(rho_init,t,T))^-nX;
det_new = det(S0_trans);
det_new_T = det(ST_trans)*sol.eval(subs(rho,t,T))^-nX;

display(sprintf('Determinant from %f to %f, percent change %f',det_init,det_new,100-100*det_new/det_init));
display(sprintf('T-determinant from %f to %f, percent change %f',det_init_T,det_new_T,100-100*det_new_T/det_init_T));

if (det_new > det_init) && (det_new_T > det_init_T)
%   keyboard
end
end

function [cost,cost_const] = calc_cost(nX,t,T,S0,ST,rho,Q_init,Q_init_T,rho_init,scale_mat)
T_scale = 1;

rho_init_T = double(subs(rho_init,t,T));
S0_trans = scale_mat*S0*scale_mat;
ST_trans = scale_mat*ST*scale_mat';
Q_init_trans = scale_mat*Q_init*scale_mat';
Q_init_T_trans = scale_mat*Q_init_T*scale_mat';

det_init = det(Q_init_trans);
det_init_T = det(Q_init_T_trans);

cost_coeffs = det_init*inv(Q_init_trans);
cost_coeffs_T = T_scale*det_init_T*inv(Q_init_T_trans);

cost_rho = -nX/rho_init_T*det_init_T*subs(rho,t,T);

cost = sum(sum((S0_trans-Q_init_trans).*cost_coeffs));
cost = cost + rho_init_T^(-nX)*(sum(sum((ST_trans - Q_init_T_trans).*cost_coeffs_T)) + cost_rho);

cost = cost/norm([cost_coeffs(:);cost_coeffs_T(:)],inf);

coeffs = cost.coeff;
pows = cost.pow;
assert(pows(1) == 0)

cost_const = coeffs(1);
cost = cost - cost_const;
end



function [cost,cost_const] = calc_integral_cost(prog,x,t,V,radius,T,scale_mat)
nX = length(x);
A_diag = ones(1,nX)*radius;
cost = spotlessIntegral(prog,V,x,A_diag,t,[0 T]);

% add cost in x-direction
% the "1" isn't quite right, 
cost_line = spotlessIntegral(prog,subs(V,[t;x],[0;x(1)*scale_mat(:,1)]),x(1),radius/norm(scale_mat(:,1)),[],[]); 

cost = cost + 1000*cost_line;

% 
if isnumeric(cost)
  cost_const = cost;
  cost = 0;
else
  coeffs = cost.coeff;
  pows = cost.pow;
  
  if any(pows == 0)
    assert(pows(1) == 0)
    cost_const = coeffs(1);
    cost = cost - cost_const;
  else
    cost_const = 0;
  end
end
end

