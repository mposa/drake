function V = strictlyFeasibleResetTest(x,s,V0,r,reset_constraint,Vprev)
nX = length(x);
V = V0;


Q_init = double(diff(diff(V,x)',x))/2;


% scale problem data based on t=0
% y = T*x
% x = T^-1 y
y = msspoly('y',nX);
scale_transform = Q_init^(.5);
x_y = inv(scale_transform)*y;
r_y = scale_transform*subs(r,x,x_y);
reset_constraint_y = subs(reset_constraint,x,x_y);
V0_y = subs(V0,x,x_y);
Q_init_y = eye(nX);
Vprev_y = subs(Vprev,x,x_y);


V = binarySearchV(y,s,r_y,reset_constraint_y,Vprev_y,Q_init_y,scale_transform);
V = subs(V,y,scale_transform*x);
% [mult,bmult] = iter1(t,x,f,g,T,umat,V0,B0,rho0);
% [V,B,rho] = iter2(t,x,s,f,g,T,r,reset_constraint,Vprev,umat,Q_init,Q_init_T,mult,bmult,rho0);


end

function V_opt = binarySearchV(x,s,r,reset_constraint,Vprev,Q_init,scale_transform)
max_iter = 5;

scale_mat = scale_transform';

for i=1:max_iter,
  prog = spotsosprog;
  prog = prog.withIndeterminate(x);
  prog = prog.withIndeterminate(s);
  rho = 1;
  
  [prog,V] = prog.newFreePoly(reshape(monomials(x,2:2),[],1));
  S = diff(diff(V,x)',x)/2;
  
  [prog,gamma] = prog.newFree(1);
  
  if i == 1,
    % initialize cost    
    [~,cost_const] = calc_cost(S,Q_init,scale_mat);
    
    if cost_const > 0
      keyboard % assuming it's negative for some of the stuff below
    end
    
    cost_min = -inf;
    cost_max = -cost_const;
    cost_val = -cost_const;
  end
  
  
  prog = prog.withSOS(V);
  
  
  reset_sos = (V - rho)*(1 + [s;x]'*[s;x])^2;
  Vprev_reset = subs(Vprev,x,r);
  %   [prog, reset_sos] = spotless_add_sprocedure(prog, reset_sos,Vprev - 1,[s;x],4);
  [prog, reset_sos] = spotless_add_sprocedure(prog, reset_sos,Vprev_reset - 1,[s;x],4);
  [prog, reset_sos] = spotless_add_eq_sprocedure(prog, reset_sos,reset_constraint,[s;x],3);
  prog = prog.withSOS(reset_sos+gamma);
  
  
  cost = calc_cost(S,Q_init,scale_mat);
  prog = prog.withPos(cost_val - cost);  
  
  spot_options = spotprog.defaultOptions;
  spot_options.verbose = true;
  spot_options.do_fr = false;
  spot_options.sos_slack = -1e-5;
  solver = @spot_mosek;
  
  sol = prog.minimize(gamma,solver,spot_options);
  
  if sol.eval(gamma) < -1e-6
    cost_max = cost_val;
    V_opt = sol.eval(V);

  else
    cost_min = cost_val;
    if i ==1,
      keyboard % uh oh
    end
  end
  
  if ~isinf(cost_min)
    cost_val = (cost_max + cost_min)/2;
  else
    cost_val = cost_val/1.1;
  end  
end

S_trans = sol.eval(scale_mat*S*scale_mat');
Q_init_trans = scale_mat*Q_init*scale_mat';

det_init = det(Q_init_trans);
det_new = det(S_trans);

display(sprintf('Determinant from %f to %f, percent change %f',det_init,det_new,100-100*det_new/det_init));

if (det_new > det_init)
  keyboard
end
end

function [cost,cost_const] = calc_cost(S,Q_init,scale_mat)

S0_trans = scale_mat*S*scale_mat;
Q_init_trans = scale_mat*Q_init*scale_mat';

det_init = det(Q_init_trans);

cost_coeffs = det_init*inv(Q_init_trans);

cost = sum(sum((S0_trans-Q_init_trans).*cost_coeffs));

cost = cost/norm(cost_coeffs(:),inf);

coeffs = cost.coeff;
pows = cost.pow;
assert(pows(1) == 0)

cost_const = coeffs(1);
cost = cost - cost_const;
end

