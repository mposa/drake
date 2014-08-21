% clear all
% things to document:
%   good lyapunov function
%   mosek providing bogus results on feasible problem
megaclear
g = 9.81;
% g = 10;
sos_option = 1;
switch sos_option
  case 1
    withSOS_fun = @withSOS;
    PSD_fun = @newPSD;
  case 2
    withSOS_fun = @withSDSOS;
    PSD_fun = @newSDD;
  case 3
    withSOS_fun = @withDSOS;
    PSD_fun = @newDD;
end
sos_fun = @spot_mosek;

z = msspoly('z',1);
zd = msspoly('zd',1);

y = msspoly('y',1);
yd = msspoly('yd',1);

prog = spotsosprog();
prog = prog.withIndeterminate(z);
prog = prog.withIndeterminate(zd);
prog = prog.withIndeterminate(y);
prog = prog.withIndeterminate(yd);

degree = 4;
x_vars = [z;zd;y;yd];
V_monoms = monomials(x_vars,1:degree);
% V_monoms = V_monoms(2:end-4);
[prog,V,coefv]= prog.newFreePoly(V_monoms);


% SOS equations
% (1) V >= 0 for q admissable
% (2) Vdot_free <= 0 for q admissable
% (3) Vdot_impact <= 0 for q,lx admissable, phi=0, phidot<=0

u = -10*y - yd + g;

f_free = [zd;-g;yd;u];
f_impact_1 = [0;1;0;0];
f_impact_2 = [0;1;0;-1];

Vdot_free = diff(V,x_vars)*f_free;
Vdot_impact_1 = diff(V,x_vars)*f_impact_1;
Vdot_impact_2 = diff(V,x_vars)*f_impact_2;

alph = .1*x_vars'*x_vars;

sos_1 = V - alph;
sos_2 = -Vdot_free - alph;
sos_3 = -Vdot_impact_1 - alph;
sos_4 = -Vdot_impact_2 - alph;

%% Add in constraints
phi_1 = 1 + z;
phidot_1 = zd;

phi_2 = z - y;
phidot_2 = zd - yd;

const_deg = degree;
sig = {};
coefsig = {};

const1_monoms = monomials(x_vars,0:degree);
% const1_monoms = const1_monoms(1:end-5);
% const1_monoms = const1_monoms([1:6 10:end]);


const2_monoms = monomials(x_vars,0:degree);
% const2_monoms = const2_monoms(1:end-1);
% const2_monoms = const2_monoms([1:7 10:end]);

const3_monoms = monomials(x_vars,0:degree);
% const3_monoms = const3_monoms(1:end-1);

% non-penetration (admissability of x)
[prog, sos_1, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_1, phi_1, [x_vars], const1_monoms, sos_option);
[prog, sos_2, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_2, phi_1, [x_vars], const2_monoms, sos_option);
[prog, sos_3, sig{end+1}, coefsig{end+1}] = spotless_add_eq_sprocedure(prog, sos_3, phi_1, [x_vars], const_deg);
[prog, sos_4, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_4, phi_1, [x_vars], const2_monoms, sos_option);

[prog, sos_1, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_1, phi_2, [x_vars], const1_monoms, sos_option);
[prog, sos_2, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_2, phi_2, [x_vars], const2_monoms, sos_option);
[prog, sos_3, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_3, phi_2, [x_vars], const2_monoms, sos_option);
[prog, sos_4, sig{end+1}, coefsig{end+1}] = spotless_add_eq_sprocedure(prog, sos_4, phi_2, [x_vars], const_deg);

% phidot <= 0 contact constraint
[prog, sos_3, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_3, -phidot_1, [x_vars], const3_monoms, sos_option);
[prog, sos_4, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_4, -phidot_2, [x_vars], const3_monoms, sos_option);

%% Solve program
prog = withSOS_fun(prog,sos_1);
prog = withSOS_fun(prog,sos_2);
prog = withSOS_fun(prog,sos_3);
prog = withSOS_fun(prog,sos_4);

options = spotprog.defaultOptions;
options.verbose = 1;
options.do_fr = false;
sol = prog.minimize(0,sos_fun,options);