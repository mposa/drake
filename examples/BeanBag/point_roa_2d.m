% clear all
% things to document:
%   good lyapunov function
%   mosek providing bogus results on feasible problem
megaclear
% g = 9.81;
g = 10;
sos_option = 1;
switch sos_option
  case 1
%     sos_fun = @spot_sedumi;
    sos_fun = @spot_mosek_sos;
%     sos_fun = @spot_sdpnal;
    withSOS_fun = @withSOS;
    PSD_fun = @newPSD;
  case 2
    sos_fun = @spot_mosek_sdsos;
    withSOS_fun = @withSDSOS;
    PSD_fun = @newSDD;
  case 3
    sos_fun = @spot_mosek_dsos;
%     sos_fun = @spot_gurobi_dsos;7
    withSOS_fun = @withDSOS;
    PSD_fun = @newDD;
end

z = msspoly('z',1);
zd = msspoly('zd',1);
lx = msspoly('lx',1);
xd = msspoly('xd',1);

prog = spotsosprog();
prog = prog.withIndeterminate(z);
prog = prog.withIndeterminate(zd);
prog = prog.withIndeterminate(xd);
prog = prog.withIndeterminate(lx);

degree = 4;
x_vars = [z;xd;zd];
V_monoms = monomials([x_vars],1:degree);
% V_monoms = V_monoms(2:end-4);
[prog,V,coefv]= prog.newFreePoly(V_monoms);

E = .5*zd^2 + .5*xd^2 + g*z;

% V = (E + E^2) + .5*(zd^3 + zd*z) + xd^2*zd;


% SOS equations
% (1) V >= 0 for q admissable
% (2) Vdot_free <= 0 for q admissable
% (3) Vdot_impact <= 0 for q,lx admissable, phi=0, phidot<=0

f_free = [zd;0;-g];
f_impact = [0;lx;1];

Vdot_free = diff(V,[z;xd;zd])*f_free;
Vdot_impact = diff(V,[z;xd;zd])*f_impact;

alph = .1;

sos_1 = V - alph*(z + xd^2+ zd^2);
sos_2 = -Vdot_free - alph*(z + xd^2 + zd^2);
sos_3 = -Vdot_impact - alph*(z + xd^2 + zd^2);

%% Add in constraints
phi = z;
phidot = zd;
psi = xd;

const_deg = degree;
sig = {};
coefsig = {};

const1_monoms = monomials(x_vars,0:degree);
% const1_monoms = const1_monoms(1:end-5);
% const1_monoms = const1_monoms([1:6 10:end]);


const2_monoms = monomials(x_vars,0:degree);
% const2_monoms = const2_monoms(1:end-1);
% const2_monoms = const2_monoms([1:7 10:end]);

const3_monoms = monomials([x_vars;lx],0:degree);
% const3_monoms = const3_monoms(1:end-1);

% non-penetration (admissability of x)
[prog, sos_1, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_1, phi, [x_vars], const1_monoms, sos_option);
[prog, sos_2, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_2, phi, [x_vars], const2_monoms, sos_option);
[prog, sos_3, sig{end+1}, coefsig{end+1}] = spotless_add_eq_sprocedure(prog, sos_3, phi, [x_vars;lx], const_deg);

% phidot <= 0 contact constraint
[prog, sos_3, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_3, -phidot, [x_vars;lx], const3_monoms, sos_option);
[prog, sos_3, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_3, -psi*lx, [x_vars;lx], const3_monoms, sos_option);
[prog, sos_3, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_3, 1-lx^2, [x_vars;lx], const3_monoms, sos_option);
[prog, sos_3, sig{end+1}, coefsig{end+1}] = spotless_add_eq_sprocedure(prog, sos_3, psi*(1-lx^2), [x_vars;lx], const_deg);

%% Solve program
prog = withSOS_fun(prog,sos_1);
prog = withSOS_fun(prog,sos_2);
prog = withSOS_fun(prog,sos_3);

options = spotprog.defaultOptions;
options.verbose = 1;
options.do_fr = true;
sol = prog.minimize(0,@spot_mosek,options);


% %%
% prog2 = spotsosprog();
% prog2 = prog2.withIndeterminate(z);
% prog2 = prog2.withIndeterminate(zd);
% 
% [prog2,tmp] = prog2.newFree(1);
% 
% prog2 = withSOS_fun(prog2,-z^2 + 1e-10*z);
% % prog2 = withSOS_fun(prog2,sol.eval(sos_1));
% % prog2 = withSOS_fun(prog2,sol.eval(sos_2));
% % prog2 = withSOS_fun(prog2,sol.eval(sos_3));
% 
% options2 = options;
% options2.do_fr = false;
% sol2 = prog2.minimize(0,@spot_mosek,options2);