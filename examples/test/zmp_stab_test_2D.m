megaclear
degree = 6;
T = 1;

R = 1;
R_diag = [2 1 1 1 1];
RT_diag = R_diag;

RT = .01;
A = 1/R*diag(1./(R_diag.^2));
AT = 1/RT*diag(1./(RT_diag.^2));

q_cm = msspoly('q',2);
v_cm = msspoly('v',2);
x_sup = msspoly('x',2);
u = msspoly('u',6);
t = msspoly('t',1);

f_sup1 = u(1:2);
f_sup2 = u(3:4);
sup1dot = u(5);
sup2dot = u(6);

W_vars = [q_cm;v_cm;x_sup(2)];
V_vars = [q_cm;v_cm;x_sup(2);t];

prog = spotsosprog();

prog = prog.withIndeterminate(q_cm);
prog = prog.withIndeterminate(v_cm);
prog = prog.withIndeterminate(x_sup(2));
prog = prog.withIndeterminate(u);
prog = prog.withIndeterminate(t);


%% Functions V, W
[prog,W,coefw]= prog.newFreePoly(monomials(W_vars,0:degree));
[prog,V,coefv]= prog.newFreePoly(monomials(V_vars,0:degree));

%% Dynamics
% state = [q_cm;v_cm;x_sup]
% f = [v_cm; f_sup1 + f_sup2; sup1dot;sup2dot];
% offset by x_sup(1)

f = [v_cm(1) + sup1dot; v_cm(2); f_sup1 + f_sup2; sup1dot + sup2dot];

Vdot = T*diff(V,[q_cm;v_cm;x_sup(2)])*f + diff(V,t);

T = 1; %rescale time

%% SOS functions
% (1) -Vdot >= 0 for x admissable and x in B
% (2) W(x) - V(0,x) - 1 >= 0 for (x) admissable and in B
% (3) V(x,T) >= 0 for x in BT and admissable
% (4) W >= 0 for x in B

sos_1 = -Vdot;
sos_2 = W - subs(V,t,0) - 1;
sos_3 = subs(V,t,T);  % fixed end time
sos_4 = W;

%% Add in constraints
const_deg = degree;
sig = {};
coefsig = {};

% constraints on U
% f_z >= 0
mu = 1;
[prog, sos_1, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_1, f_sup1(2), [V_vars], const_deg);
[prog, sos_1, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_1, f_sup2(2), [V_vars], const_deg);
% friction cone
[prog, sos_1, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_1, f_sup1(2)^2 - mu^2*f_sup1(1)^2, [V_vars], const_deg);
[prog, sos_1, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_1, f_sup2(2)^2 - mu^2*f_sup2(1)^2, [V_vars], const_deg);
% speed limits
[prog, sos_1, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_1, 1 - sup1dot^2, [V_vars], const_deg);
[prog, sos_1, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_1, 1 - sup2dot^2, [V_vars], const_deg);

%no moving and force at the same time
[prog, sos_1, sig{end+1}, coefsig{end+1}] = spotless_add_eq_sprocedure(prog, sos_1, f_sup1(2)*sup1dot, [V_vars], const_deg);
[prog, sos_1, sig{end+1}, coefsig{end+1}] = spotless_add_eq_sprocedure(prog, sos_1, f_sup2(2)*sup2dot, [V_vars], const_deg);

% Ball constraints
h_B = 1 - W_vars'*A*W_vars;
h_BT = 1 - W_vars'*AT*W_vars;

[prog, sos_1, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_1, h_B, V_vars, const_deg);
[prog, sos_2, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_2, h_B, V_vars, const_deg);
[prog, sos_3, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_3, h_BT, W_vars, const_deg);
[prog, sos_4, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_4, h_B, W_vars, const_deg);

%% Setup cost function
[vars,alphas,coeff] = decomp(W, prog.freeVar);

assert(isequal(W_vars,vars([1;3;2;4;5])))

alphas = alphas(:,[1 3 2 4 5]);

nX = 5;

% [~,alphas] = monomials(x_vars,0:degree);
betas = 0.5*(alphas + 1);
Ra = (R.^(sum(alphas,2) + nX))./(sum(alphas,2) + nX);
IS = 2*prod(gamma(betas),2)./(gamma(sum(betas,2)));
l = Ra.*IS;
alphaszero = (mod(alphas,2) ~= 0);
alphaszero = any(alphaszero,2);
l(alphaszero) = 0;
l = l.*prod(repmat(R_diag,size(alphas,1),1).^(alphas+1),2);

%% Solve program
prog = prog.withSOS(sos_1);
prog = prog.withSOS(sos_2);
prog = prog.withSOS(sos_3);
prog = prog.withSOS(sos_4);

options = spotprog.defaultOptions;
options.verbose = 1;
options.do_fr = false;
sol = prog.minimize(coeff*l,@spot_mosek,options);

%% plotting
close all
Vsol = sol.eval(V)

[X_CM,Z_CM] = meshgrid(linspace(-R_diag(1),R_diag(1),100),linspace(-R_diag(2),R_diag(2),100));

Vtmp = sol.eval(subs(V,[t;v_cm;x_sup(2)],[1;0;0;0]));
Vval = msubs(Vtmp,q_cm,[X_CM(:)'; Z_CM(:)']);
Vval = reshape(Vval,size(X_CM,1),[]);
figure(1)
[cl,h] = contour(X_CM,Z_CM,Vval);
clabel(cl,h);
