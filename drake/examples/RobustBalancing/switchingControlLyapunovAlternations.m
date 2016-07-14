function [ V,B ] = switchingControlLyapunovAlternations(x,f,g,V0,B0,R)
if nargin < 6
  R = [];
end

N = 1;
nX = length(x);
nU = size(g,2);

ndgrid_arg = mat2cell(repmat([-1;1],1,nU),2,ones(1,nU)');
[ugrid{1:nU}] = ndgrid(ndgrid_arg{:});
umat = zeros(2^nU,nU);
for i=1:nU,
  umat(:,i) = ugrid{i}(:);
end


if nargin < 5
  R = [];
end

Q_init = double(diff(diff(V0,x)',x))/2;

% scale problem data
% y = T*x
% x = T^-1 y
y = msspoly('y',nX);
T = Q_init^(.5);
x_y = inv(T)*y;
f_y = T*subs(f,x,x_y);
g_y = T*subs(g,x,x_y);
% V0_y = subs(V0,x,x_y);
V0_y = y'*y;
B0_y = subs(B0,x,x_y);
R_y = cell(length(R),1);
for i=1:length(R),
  R_y{i} = inv(T)'*R{i}*inv(T);
end
Q_init_y = eye(nX);

[mult,bmult] = iter1(y,f_y,g_y,umat,V0_y,B0_y);
[V_y,B_y] = iter2(y,f_y,g_y,umat,Q_init_y,R_y,T,mult,bmult);
V = subs(V_y,y,T*x);
B = subs(B_y,y,T*x);

% [mult,bmult] = iter1(x_mss,f,g,umat,V0,B0);
% [V,B] = iter2(x_mss,f,g,umat,Q_init,R,mult,bmult);

end
%% iter 1
function [mult,bmult] = iter1(x,f,g,umat,V,B)
nX = length(x);
nU = size(g,2);
prog = spotsosprog;
prog = prog.withIndeterminate(x);
[prog,gamma] = prog.newPos(1);
prog = prog.withPos(1 - gamma);

for j=1:2^nU
  Vdot = diff(V,x)*(f + g*umat(j,:)')*(1+x'*x)^0;
  [prog, Vdot_sos,mult{j},coeff] = spotless_add_eq_sprocedure(prog, -Vdot, 1-V,x,4);
  for k=1:nU,
    [prog, Vdot_sos,bmult{j}{k},coeff] = spotless_add_sprocedure(prog, Vdot_sos, umat(j,k)*B(k),x,4);
  end
  prog = prog.withSOS(Vdot_sos - gamma*(x'*x)^2);
end

spot_options = spotprog.defaultOptions;
spot_options.verbose = true;
%   spot_options.do_fr = true;
spot_options.clean_primal = true;
solver = @spot_mosek;
sol = prog.minimize(-gamma,solver,spot_options);

for j=1:2^nU
  mult{j} = sol.eval(mult{j});
  for k=1:nU,
    bmult{j}{k} = sol.eval(bmult{j}{k});
  end
end
end


function [V,B] = iter2(x,f,g,umat,Q_init,R,T,mult,bmult)
%% iter 2
nX = length(x);
nU = size(g,2);
prog = spotsosprog;
prog = prog.withIndeterminate(x);
rho = 1;
[prog,Q] = prog.newPSD(nX);
V = x'*Q*x;
%   V = V0;

[prog,B] = prog.newFreePoly(monomials(x,1:2),nU);
%   B = B0;

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
  Vdot = diff(V,x)*(f + g*umat(j,:)')*(1+x'*x)^0;
  Vdot_sos = -Vdot - mult{j}*(rho-V);
  for k=1:nU,
    Vdot_sos = Vdot_sos - bmult{j}{k}*umat(j,k)*B(k);
  end
  prog = prog.withSOS(Vdot_sos + 1e-8*(x'*x));
end

scale_mat = eye(length(x));
% scale_mat(1,1) = 20;
scale_mat = T'*scale_mat;

cost_coeffs = det(scale_mat*Q_init*scale_mat')*inv(scale_mat*Q_init*scale_mat');
cost_coeffs = cost_coeffs/norm(cost_coeffs(:),inf);
cost = sum(sum(scale_mat*(Q-Q_init)*scale_mat'.*cost_coeffs));

%
Q_full = T'*Q*T;
cost = cost + 10*Q_full(1,1);  % should do this before scaling. see strictly feasible script

spot_options = spotprog.defaultOptions;
spot_options.verbose = true;
%   spot_options.do_fr = true;
spot_options.clean_primal = true;
solver = @spot_mosek;


sol = prog.minimize(cost,solver,spot_options);

Q_init_det = det(scale_mat*Q_init*scale_mat');
Q_det = det(sol.eval(scale_mat*Q*scale_mat'));
display(sprintf('Determinant from %f to %f, percent change %f',Q_init_det,Q_det,100 - 100*Q_det/Q_init_det));
V = sol.eval(V);
B = sol.eval(B);
%   u = sol.eval(u)
end

