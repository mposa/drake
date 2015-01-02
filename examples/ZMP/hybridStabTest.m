megaclear
prog = spotsosprog();

mu = 1;

x = msspoly('x',1);
xd = msspoly('xd',1);

% z = msspoly('z',1);
% zd = msspoly('zd',1);

% xf = msspoly('xf',1);

prog = prog.withIndeterminate(x);
prog = prog.withIndeterminate(xd);
% prog = prog.withIndeterminate(z);
% prog = prog.withIndeterminate(zd);
% prog = prog.withIndeterminate(xf);

% vars = [x;xd;z;zd;xf];
vars = [x;xd];

degree = 1;

% [prog,lx1,coeflx1] = prog.newFreePoly(vars,degree);
% [prog,lz1,coeflz1] = prog.newFreePoly(vars,degree);
% [prog,lx2,coeflx2] = prog.newFreePoly(vars,degree);
% [prog,lzd,coeflz2] = prog.newFreePoly(vars,degree);

[prog,u,coefu] = prog.newFreePoly(monomials(vars,0:degree));

xf1 = .1;
xf2 = -.1;

f = [xd; xf1*u + xf2*(1-u)];

V = xd^2 + x^2;

Vdot = diff(V,vars)*f;

rho = .01;

sos_1 = -Vdot;
sos_2 = u + 1;
sos_3 = 1 - u;

const_deg = degree;
sig = {};
coefsig = {};


% non-penetration (admissability of x)
[prog, sos_1, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_1, rho - V,vars, const_deg);
[prog, sos_2, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_2, rho - V,vars, const_deg);
[prog, sos_3, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_3, rho - V,vars, const_deg);

%% Solve program
prog = prog.withSOS(sos_1);
prog = prog.withSOS(sos_2);
prog = prog.withSOS(sos_3);
options = spotprog.defaultOptions;
options.verbose = true;
sol = prog.minimize(0,@spot_mosek,options);
