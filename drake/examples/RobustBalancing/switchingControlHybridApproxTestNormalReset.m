function [V,rho]=switchingControlHybridApproxTestNormalReset(x_mss,s_mss,V0,rho0,T,r,reset_constraint,Vprev)
t = msspoly('t',1);

Q_init = double(diff(diff(subs(V0,t,0),x_mss)',x_mss))/2;
Q_init_T = double(diff(diff(subs(V0,t,T),x_mss)',x_mss))/2;
rho_T_nom = double(subs(rho0,t,T));


prog = spotsosprog;
prog = prog.withIndeterminate(x_mss);
% prog = prog.withIndeterminate(t);
prog = prog.withIndeterminate(s_mss);


% [prog,rho] = prog.newFreePoly(monomials(t,1:2));
% rho = rho  + 1; % set rho(0) = 1

rho = 1;

[prog,V] = prog.newFreePoly(monomials(x_mss,2:2));
S0 = diff(diff(subs(V,t,0),x_mss)',x_mss)/2;
ST = diff(diff(subs(V,t,T),x_mss)',x_mss)/2;

V0_reset = subs(Vprev,x_mss,r);
reset_sos = subs(V - rho,t,t)*(1 + [s_mss;x_mss]'*[s_mss;x_mss])^2;
[prog, reset_sos] = spotless_add_sprocedure(prog, reset_sos,V0_reset - 1,[s_mss;x_mss],4);
[prog, reset_sos] = spotless_add_eq_sprocedure(prog, reset_sos,reset_constraint,[s_mss;x_mss],4);
prog = prog.withSOS(reset_sos);

% [prog, V_sos] = spotless_add_sprocedure(prog, V*(1+[t;s_mss;x_mss]'*[t;s_mss;x_mss])^0,t*(T-t),[t;s_mss;x_mss],2);
% prog = prog.withSOS(V_sos);

scale_diag = ones(length(x_mss),1);
scale_diag(1) = 10;
scale_diag(3) = 10;
scale_mat = diag(scale_diag);



det_init = det(scale_mat*Q_init*scale_mat');
det_init_T = det(scale_mat*Q_init_T*scale_mat');

% linearization of determinant
cost_coeffs = det_init*inv(scale_mat*Q_init*scale_mat');
cost_coeffs_T = det_init_T*inv(scale_mat*Q_init_T*scale_mat');


cost_rho = -0*length(x_mss)/rho_T_nom*det_init_T*subs(rho,t,T);

cost = 1*sum(sum(scale_mat*(S0-Q_init)*scale_mat'.*cost_coeffs));
cost = cost + 0*rho_T_nom^(-length(x_mss))*(sum(sum(scale_mat*(ST-Q_init_T)*scale_mat'.*cost_coeffs_T)) + cost_rho);

cost = cost/norm(cost_coeffs(:),inf);

% cost = trace(S0) + trace(ST);

% cost = S0(1) + ST(1);

spot_options = spotprog.defaultOptions;
spot_options.verbose = true;
spot_options.do_fr = false;
spot_options.clean_primal = false;
solver = @spot_mosek;
sol = prog.minimize(cost,solver,spot_options);

V = sol.eval(V);
rho = sol.eval(rho);
det(sol.eval(scale_mat*S0*scale_mat))

end