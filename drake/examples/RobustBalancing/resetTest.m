function V = resetTest(x,s,V0,r,reset_constraint,Vprev)
nX = length(x);


Q_init = double(diff(diff(V0,x)',x))/2;


% scale problem data based on t=0
% y = T*x
% x = T^-1 y
y = msspoly('y',nX);
scale_transform = eye(nX);

V = searchV(x,s,r,reset_constraint,Vprev,Q_init,scale_transform);
% V = subs(V,y,scale_transform*x);
% [mult,bmult] = iter1(t,x,f,g,T,umat,V0,B0,rho0);
% [V,B,rho] = iter2(t,x,s,f,g,T,r,reset_constraint,Vprev,umat,Q_init,Q_init_T,mult,bmult,rho0);


end

function V_opt = searchV(x,s,r,reset_constraint,Vprev,Q_init,scale_transform)

scale_mat = scale_transform';

prog = spotsosprog;
prog = prog.withIndeterminate(x);
prog = prog.withIndeterminate(s);
rho = 1;

[prog,V] = prog.newFreePoly(reshape(monomials(x,2:2),[],1));
S = diff(diff(V,x)',x)/2;

prog = prog.withSOS(V);


reset_sos = (V - rho)*(1 + [s;x]'*[s;x])^2;
Vprev_reset = subs(Vprev,x,r);
%   [prog, reset_sos] = spotless_add_sprocedure(prog, reset_sos,Vprev - 1,[s;x],4);
[prog, reset_sos] = spotless_add_sprocedure(prog, reset_sos,Vprev_reset - 1,[s;x],4);
[prog, reset_sos] = spotless_add_eq_sprocedure(prog, reset_sos,reset_constraint,[s;x],3);
prog = prog.withSOS(reset_sos);


cost = calc_cost(S,Q_init,scale_mat);

spot_options = spotprog.defaultOptions;
spot_options.verbose = true;
spot_options.do_fr = false;
spot_options.sos_slack = 0;
solver = @spot_mosek;

sol = prog.minimize(cost,solver,spot_options);


S_trans = sol.eval(scale_mat*S*scale_mat');
Q_init_trans = scale_mat*Q_init*scale_mat';

det_init = det(Q_init_trans);
det_new = det(S_trans);

display(sprintf('Determinant from %f to %f, percent change %f',det_init,det_new,100-100*det_new/det_init));

if (det_new > det_init)
  keyboard
end
V_opt = sol.eval(V);
end

function [cost,cost_const] = calc_cost(S,Q_init,scale_mat)


cost = trace(S);
end

