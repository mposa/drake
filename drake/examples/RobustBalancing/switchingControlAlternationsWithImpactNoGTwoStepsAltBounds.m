function [ V,B,rho] = switchingControlAlternationsWithImpactNoGTwoStepsAltBounds(x,s,f,g,V0,rho0,B0,T,r,reset_constraint,Vprev)
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
x_y = inv(scale_transform)*y;
f_y = scale_transform*subs(f,x,x_y);
g_y = scale_transform*subs(g,x,x_y);
r_y = scale_transform*subs(r,x,x_y);
reset_constraint_y = subs(reset_constraint,x,x_y);
V0_y = subs(V0,x,x_y);
B0_y = subs(B0,x,x_y);
Q_init_y = eye(nX);
Q_init_T_y = inv(scale_transform)'*Q_init_T*inv(scale_transform);
Vprev_y = subs(Vprev,x,x_y);


[mult,bmult] = iter1(t,y,f_y,g_y,T,umat,V0_y,B0_y,rho0);
[V0_y,B0_y,rho] = iter2(t,y,s,f_y,g_y,T,r_y,reset_constraint_y,Vprev_y,umat,Q_init_y,Q_init_T_y,scale_transform,mult,bmult,rho0);
V = subs(V0_y,y,scale_transform*x);
B = subs(B0_y,y,scale_transform*x);

% [mult,bmult] = iter1(t,x,f,g,T,umat,V0,B0,rho0);
% [V,B,rho] = iter2(t,x,s,f,g,T,r,reset_constraint,Vprev,umat,Q_init,Q_init_T,mult,bmult,rho0);


end


function [mult,bmult] = iter1(t,x,f,g,T,umat,V,B,rho)
%% step 1
nX = length(x);
nU = size(g,2);
prog = spotsosprog;
prog = prog.withIndeterminate(x);
prog = prog.withIndeterminate(t);
sproc_vars = [t;x];

[prog,gamma] = prog.newPos(1);
prog = prog.withPos(1 - gamma);

rhodot = diff(rho,t);

for j=1:2^nU
  Vdot = diff(V,x)*(f + g*umat(j,:)') + diff(V,t);
  
  [prog, Vdot_sos,mult{j},coeff] = spotless_add_eq_sprocedure(prog, rhodot-Vdot, rho-V,x,2);
  [prog, Vdot_sos] = spotless_add_sprocedure(prog, Vdot_sos, t*(T-t),sproc_vars,4);
  for k=1:nU,
    [prog, Vdot_sos,bmult{j}{k},coeff] = spotless_add_sprocedure(prog, Vdot_sos, umat(j,k)*B(k),x,4);
  end
  Vdot_degree = even_degree(Vdot_sos,sproc_vars)-2;
  prog = prog.withSOS(Vdot_sos - gamma*(x'*x)^(Vdot_degree/2));
end

%   rad_sos = subs(V-rho,t,T)*(1+x'*x);
%   [prog, rad_sos] = spotless_add_sprocedure(prog, rad_sos, 4*a*d-b^2,x,2);
%   [prog, rad_sos] = spotless_add_sprocedure(prog, rad_sos, b,x,2);
%   prog = prog.withSOS(rad_sos);
%
%   ab_sos = subs(V-rho,t,T)*(1+x'*x);
%   [prog, ab_sos] = spotless_add_sprocedure(prog, ab_sos, 2*a-b,x,2);
%   [prog, ab_sos] = spotless_add_sprocedure(prog, ab_sos, -a+b-d,x,2);
%   [prog, ab_sos] = spotless_add_sprocedure(prog, ab_sos, -b,x,2);
%   prog = prog.withSOS(ab_sos);

spot_options = spotprog.defaultOptions;
spot_options.verbose = false;
spot_options.do_fr = false;
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

function [V,B,rho] = iter2(t,x,s,f,g,T,r,reset_constraint,Vprev,umat,Q_init,Q_init_T,scale_transform,mult,bmult,rho0)
%% step 2 of 2
nX = length(x);
nU = size(g,2);
prog = spotsosprog;
prog = prog.withIndeterminate(x);
prog = prog.withIndeterminate(t);
prog = prog.withIndeterminate(s);
sproc_vars = [t;x];

%   [prog,gamma] = prog.newPos(1);

rho_T_nom = double(subs(rho0,t,T));

[prog,rho] = prog.newFreePoly(monomials(t,1:2));
rho = rho  + 1; % set rho(0) = 1

[prog,V] = prog.newFreePoly(reshape(monomials(x,2:2)*monomials(t,(0:2))',[],1));
[prog,B] = prog.newFreePoly(monomials(x,1:2),nU);
S0 = diff(diff(subs(V,t,0),x)',x)/2;
ST = diff(diff(subs(V,t,T),x)',x)/2;
%   [prog,Q] = prog.newPSD(nX);
%   V = x'*Q*x;
%   S0 = Q;
%   ST = Q;


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
  Vdot_degree = even_degree(Vdot_sos,sproc_vars);
  prog = prog.withSOS(Vdot_sos + 1e-6*(x'*x)^(Vdot_degree/2));
end

reset_sos = subs(V - rho,t,T)*(1 + [s;x]'*[s;x])^2;
Vprev_reset = subs(Vprev,x,r);
%   [prog, reset_sos] = spotless_add_sprocedure(prog, reset_sos,Vprev - 1,[s;x],4);
[prog, reset_sos] = spotless_add_sprocedure(prog, reset_sos,Vprev_reset - 1,[s;x],4);
[prog, reset_sos] = spotless_add_eq_sprocedure(prog, reset_sos,reset_constraint,[s;x],3);
prog = prog.withSOS(reset_sos);

spot_options = spotprog.defaultOptions;
spot_options.verbose = true;
spot_options.do_fr = false;
spot_options.sos_slack = 0e-6;
solver = @spot_mosek;

%   cost = 0.01*trace(S0)+5*S0(2,2);

scale_mat = eye(length(x));
scale_mat(1) = 1;
scale_mat = scale_transform'*scale_mat;
%   scale_mat(5) = 0;

det_init = det(scale_mat*Q_init*scale_mat');
det_init_T = det(scale_mat*Q_init_T*scale_mat');

% linearization of determinant
T_scale = 1;

cost_coeffs = det_init*inv(scale_mat*Q_init*scale_mat');
cost_coeffs_T = T_scale*det_init_T*inv(scale_mat*Q_init_T*scale_mat');

cost_rho = -length(x)/rho_T_nom*det_init_T*subs(rho,t,T);

cost = 1*sum(sum(scale_mat*(S0-Q_init)*scale_mat'.*cost_coeffs));
cost = cost + rho_T_nom^(-length(x))*(sum(sum(scale_mat*(ST-Q_init_T)*scale_mat'.*cost_coeffs_T)) + cost_rho);

cost = cost/norm([cost_coeffs(:);cost_coeffs_T(:)],inf);

%   xstar = [sqrt(10);0;0;-1;zeros(2,1)];
%   cost = cost + .01*xstar'*S0*xstar;

%   cost = cost + 1*S0(1,1);

sol = prog.minimize(cost,solver,spot_options);

det_new = det(scale_mat*double(sol.eval(S0))*scale_mat');
det_new_T = det(scale_mat*double(sol.eval(ST))*scale_mat');
display(sprintf('Determinant from %f to %f, percent change %f',det_init,det_new,100-100*det_new/det_init));
display(sprintf('T-determinant from %f to %f, percent change %f',det_init_T,det_new_T,100-100*det_new_T/det_init_T));

V = sol.eval(V);
rho = sol.eval(rho);
B = sol.eval(B);

end
