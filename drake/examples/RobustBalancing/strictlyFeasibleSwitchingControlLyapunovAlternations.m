function [ V,B ] = strictlyFeasibleSwitchingControlLyapunovAlternations(x,f,g,V0,B0,R)
if nargin < 6
  R = [];
end

degree = 2;

N = 1;
nX = length(x);
nU = size(g,2);

Q_init = double(subs(diff(diff(V0,x)',x)/2,x,zeros(nX,1)));

V = V0;
B = B0;
ndgrid_arg = mat2cell(repmat([-1;1],1,nU),2,ones(1,nU)');
[ugrid{1:nU}] = ndgrid(ndgrid_arg{:});
umat = zeros(2^nU,nU);
for i=1:nU,
  umat(:,i) = ugrid{i}(:);
end


if nargin < 5
  R = [];
end

% scale problem data
% y = T*x
% x = T^-1 y
y = msspoly('y',nX);
T = Q_init^(.5);
x_y = inv(T)*y;
f_y = T*subs(f,x,x_y);
g_y = T*subs(g,x,x_y);
V_y = subs(V,x,x_y);
B_y = subs(B,x,x_y);
R_y = cell(length(R),1);
for i=1:length(R),
  R_y{i} = inv(T)'*R{i}*inv(T);
end
Q_init_y = eye(nX);

%% iter 1
outer_radius = 2;
inner_radius = .1;
rho0 = 1;
% [rho,mult,bmult] = binarySearchRho(x,f,g,umat,V,B,outer_radius,inner_radius,rho0);
[rho_y,mult_y,bmult_y] = binarySearchRho(y,f_y,g_y,umat,V_y,B_y,outer_radius,inner_radius,rho0,T);
V_y = V_y/rho_y;

%% iter 2
% [V,B] = binarySearchVandB(x,f,g,umat,mult,bmult,outer_radius,inner_radius,Q_init,R,eye(nX));
[V_y,B_y] = binarySearchVandB(y,f_y,g_y,umat,mult_y,bmult_y,outer_radius,inner_radius,V_y,Q_init_y,R_y,T,degree);
V = subs(V_y,y,T*x);
B = subs(B_y,y,T*x);
end

function [rho_opt,mult_opt,bmult_opt] = binarySearchRho(x,f,g,umat,V,B,outer_radius,inner_radius,rho0,T)
max_iter = 4;
rho = rho0;
rho_min = -inf;
rho_max = inf;
nU = size(g,2);

for i=1:max_iter
  %%
  prog = spotsosprog;
  prog = prog.withIndeterminate(x);
  [prog,gamma] = prog.newFree(1);
  
  for j=1:2^nU
    Vdot = diff(V,x)*(f + g*umat(j,:)');
    [prog, Vdot_sos,mult{j},coeff] = spotless_add_eq_sprocedure(prog, -Vdot, rho-V,x,4);
    for k=1:nU,
      [prog, Vdot_sos,bmult{j}{k},coeff] = spotless_add_sprocedure(prog, Vdot_sos, umat(j,k)*B(k),x,4);
    end
%     [prog, Vdot_sos] = spotless_add_sprocedure(prog, Vdot_sos, outer_radius-x'*x,x,4);
%     [prog, Vdot_sos] = spotless_add_sprocedure(prog, Vdot_sos, -inner_radius+x'*x,x,4);
    prog = prog.withSOS(Vdot_sos+gamma);
  end
  
  spot_options = spotprog.defaultOptions;
  spot_options.verbose = true;
  spot_options.sos_slack = -1e-6;
  spot_options.clean_primal = false;
  solver = @spot_mosek;
  sol = prog.minimize(gamma,solver,spot_options);
  
  if sol.eval(gamma) < -1e-6
    rho_min = rho;
    for j=1:2^nU
      mult_opt{j} = sol.eval(mult{j});
      for k=1:nU,
        bmult_opt{j}{k} = sol.eval(bmult{j}{k});
      end
    end
  else
    rho_max = rho;
  end
  
  if isinf(rho_max)
    rho = 1.01*rho;
  elseif isinf(rho_min)
    rho = .98*rho;
  else
    rho = (rho_max + rho_min)/2;
  end
end
if isinf(rho_min)
  keyboard % uh oh
else
  rho_opt = rho_min;
end
end


function [V,B] = binarySearchVandB(x,f,g,umat,mult,bmult,outer_radius,inner_radius,V_init,Q_init,R,T,degree)
%%
% cost_option
% 1 - determinant
% 2 - integral
cost_option = 2; 
max_iter = 4;

nX = length(x);
nU = size(g,2);

% scale_mat = eye(nX);
scale_mat = T';

% c(Q) = f(Q) + g
% c'(Q) = f(Q) = c(Q) - g
% c(Q) < 0 <==> c'(Q) < -g

for i=1:max_iter
  
  prog = spotsosprog;
  prog = prog.withIndeterminate(x);
  rho = 1;
  [prog,gamma] = prog.newFree(1);
  
  if degree == 2
    [prog,Q] = prog.newPSD(nX);
    V = x'*Q*x;
  else
    [prog,V] = prog.newFreePoly(monomials(x,2:degree));    
    % could add s-procedure here, if needed
    prog = prog.withSOS(V);
    
    Q = subs(diff(diff(V,x)',x)/2,x,zeros(nX,1));
  end
  
  [prog,B] = prog.newFreePoly(monomials(x,1:2),nU);
  %   B = B0; 
  R = [];
  if ~isempty(R)
    if iscell(R)
      for j=1:length(R),
        prog = prog.withPSD(Q-R{j});
      end
    else
      prog = prog.withPSD(Q-R);
    end
  end
  
  for j=1:2^nU
    Vdot = diff(V,x)*(f + g*umat(j,:)');
    Vdot_sos = -Vdot - mult{j}*(rho-V);
    for k=1:nU,
      Vdot_sos = Vdot_sos - bmult{j}{k}*umat(j,k)*B(k);
    end
    [prog, Vdot_sos] = spotless_add_sprocedure(prog, Vdot_sos, outer_radius-x'*x,x,4);
    [prog, Vdot_sos] = spotless_add_sprocedure(prog, Vdot_sos, -inner_radius+x'*x,x,4);
    prog = prog.withSOS(Vdot_sos+gamma);
  end
  
  if i==1,
    % initialize cost
    if cost_option == 1
      [~,cost_const] = calc_cost(Q,Q_init,scale_mat);
      
      % based on the fact that cost(V0) = 0
      % gets around the fact that subs(Q) doesn't work with Q PSD
      cost_init = -cost_const;      
    else
      [~,cost_const] = calc_integral_cost(prog,x,V_init,outer_radius,scale_mat);      
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

  if cost_option == 1
    cost = calc_cost(Q,Q_init,scale_mat);
  else
    cost = calc_integral_cost(prog,x,V,outer_radius,scale_mat);
  end
  
  prog = prog.withPos(cost_val - cost);
  
  spot_options = spotprog.defaultOptions;
  spot_options.verbose = false;
  spot_options.sos_slack = -1e-6;
  spot_options.clean_primal = false;
  solver = @spot_mosek;
  sol = prog.minimize(gamma,solver,spot_options);
  
  if sol.eval(gamma) < -1e-6
    cost_max = cost_val;
    V_opt = sol.eval(V);
    B_opt = sol.eval(B);
    Q_opt = sol.eval(Q);
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

Q_init_det = det(scale_mat*Q_init*scale_mat');
Q_det = det(scale_mat*Q_opt*scale_mat');
display(sprintf('Determinant from %f to %f, percent change %f',Q_init_det,Q_det,100 - 100*Q_det/Q_init_det));
V = V_opt;
B = B_opt;
end

function [cost,cost_const] = calc_cost(Q,Q_init,scale_mat)
Q_trans = scale_mat*Q*scale_mat';
Q_init_trans = scale_mat*Q_init*scale_mat';

cost_coeffs = det(Q_init_trans)*inv(Q_init_trans);

% add cost on Q(1,1)
cost_coeffs = cost_coeffs + 5*scale_mat(:,1)*scale_mat(1,:);

cost_coeffs = cost_coeffs/norm(cost_coeffs(:),inf);
cost = sum(sum((Q_trans - Q_init_trans).*cost_coeffs));


coeffs = cost.coeff;
pows = cost.pow;
assert(pows(1) == 0)

cost_const = coeffs(1);
cost = cost - cost_const;
end

function [cost,cost_const] = calc_integral_cost(prog,x,V,radius,scale_mat)
nX = length(x);
A_diag = ones(1,nX)*radius;
cost = spotlessIntegral(prog,V,x,A_diag,[],[]);

% add cost in x-direction
% the "1" isn't quite right, 
cost_line = spotlessIntegral(prog,subs(V,x,x(1)*scale_mat(:,1)),x(1),radius/norm(scale_mat(:,1)),[],[]); 

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