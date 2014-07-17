megaclear
degree = 4;
T = .5;
impact_time_scale = .1;
R = 1; % don't change
R_diag = [1 1 2 2];
RT_diag = R_diag;

RT = .01;
A = 1/R*diag(1./(R_diag.^2));
AT = 1/RT*diag(1./(RT_diag.^2));

g = 2;

prog = spotsosprog();

q = msspoly('q',2);
v = msspoly('v',2);
lx = msspoly('lx',1);
t = msspoly('t',1);

lzsq = 1;

W_vars = [q;v];
V_vars = [q;v;t];

prog = prog.withIndeterminate(q);
prog = prog.withIndeterminate(v);
prog = prog.withIndeterminate(lx);
prog = prog.withIndeterminate(t);

%% Functions V, W
[prog,W,coefw]= prog.newFreePoly(monomials(W_vars,0:degree));
[prog,V,coefv]= prog.newFreePoly(monomials(V_vars,0:degree));

%% Dynamics
f_free = [v;0;-g] + [0;0;-20*q(1);-0*q(2)];
f_impact = [0;0;lx;lzsq^2];

phi = q(2);
phidot = v(2);
psi = v(1);

Vdot_free = T*diff(V,[q;v])*f_free + diff(V,t);
Vdot_impact = diff(V,[q;v])*f_impact + impact_time_scale*diff(V,t);

T = 1; % rescale time

%% SOS functions
% (1) -Vdot_free(x) >= 0 for x admissable and x in B
% (2) -Vdot_impact_1(x,l) >= 0 for (x,l) admissable and x in B
% (3) W(x) - V(0,x) - 1 >= 0 for (x) admissable and in B
% (4) V(x,T) >= 0 for x in BT and admissable
% (5) w >= 0 for x in B

sos_1 = -Vdot_free;
sos_2 = -Vdot_impact;
sos_3 = W - subs(V,t,0) - 1;
sos_4 = subs(V,t,T);  % fixed end time
sos_5 = W;


%% Add in constraints
const_deg = degree;
sig = {};
coefsig = {};

% non-penetration (admissability of x)
[prog, sos_1, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_1, phi, [V_vars], const_deg);
[prog, sos_3, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_3, phi, W_vars, const_deg);
[prog, sos_4, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_4, phi, W_vars, const_deg);

% Ball constraints
h_B = 1 - [q;v]'*A*[q;v];
h_BT = 1 - [q;v]'*AT*[q;v];

[prog, sos_1, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_1, h_B, [V_vars], const_deg);
[prog, sos_2, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_2, h_B, [V_vars;lx(1)], const_deg);
[prog, sos_3, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_3, h_B, W_vars, const_deg);
[prog, sos_4, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_4, h_BT, W_vars, const_deg);
[prog, sos_5, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_5, h_B, W_vars, const_deg);


% Contact constraints (admissability of lambda)
[prog, sos_2, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_2, -phidot, [V_vars;lx(1)], const_deg);
[prog, sos_2, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_2, -lx(1)*psi, [V_vars;lx(1)], const_deg-1);
[prog, sos_2, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_2, lzsq(1)^2 - lx(1)^2, [V_vars;lx(1)], const_deg);
[prog, sos_2, sig{end+1}, coefsig{end+1}] = spotless_add_eq_sprocedure(prog, sos_2, phi, [V_vars;lx(1)], const_deg+1);
[prog, sos_2, sig{end+1}, coefsig{end+1}] = spotless_add_eq_sprocedure(prog, sos_2, (lzsq(1)^2 - lx(1)^2)*psi, [V_vars;lx(1)], const_deg-2);  %should this be psi^2?

%% Setup cost function
[vars,alphas,coeff] = decomp(W, prog.freeVar);
assert(isequal(vars([1;3]),q))
assert(isequal(vars([2 4]),v))

alphas = alphas(:,[1 3 2 4]);

nX = 4;

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
prog = prog.withSOS(sos_5);

options = spotprog.defaultOptions;
options.verbose = 1;
options.do_fr = false;
sol = prog.minimize(coeff*l,@spot_mosek,options);

%% plotting
close all
Vsol = sol.eval(V)
[X,Z] = meshgrid(linspace(-R_diag(1),R_diag(1),100),linspace(-.1*R_diag(2),R_diag(2),100));


Vtmp = sol.eval(subs(V,[t;v],[0;0;0]));
Vval = msubs(Vtmp,q,[X(:)'; Z(:)']);
Vval = reshape(Vval,size(X,1),[]);
figure(1)
[cl,h] = contour(X,Z,Vval);
clabel(cl,h);

Wtmp = sol.eval(subs(W,[t;v],[0;0;0]));
Wval = msubs(Wtmp,q,[X(:)'; Z(:)']);
Wval = reshape(Wval,size(X,1),[]);
figure(2)
[cl,h] = contour(X,Z,Wval);
clabel(cl,h);

Vtmp = subs(Vsol,[q(2);v(2);t],[0;0;0]);
[X,XD] = meshgrid(linspace(-R_diag(1),R_diag(1),100),linspace(-R_diag(3),R_diag(3),100));
Vval = msubs(Vtmp,[q(1);v(1)],[X(:)'; XD(:)']);
Vval = reshape(Vval,size(X,1),[]);

figure(3)
[cl,h] = contour(X,XD,Vval);
clabel(cl,h);
xlabel('x')
ylabel('xd')