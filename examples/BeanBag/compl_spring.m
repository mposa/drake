% two masses connected by two springs (m1 to the sky, m2 to m1) with
% contact between m2 and the ground, via complementarity
megaclear
degree = 8;

prog = spotsosprog();

q = msspoly('q',2);
v = msspoly('v',2);
prog = prog.withIndeterminate(q);
prog = prog.withIndeterminate(v);

m1 = 1;
m2 = 1;
k1 = 1;
k2 = 1;
b1 = 1;
b2 = 1;
ground_height = -.5;

x = [q;v];
f_free = [-1/m1*k1*q(1) - 1/m1*b1*v(1) - 1/m1*k2*(q(1) - q(2)) - 1/m1*b2*(v(1) - v(2));
          1/m2*k2*(q(1) - q(2))+1/m2*b2*(v(1) - v(2));v];
        
f_impact = [0;1/m2;0;0];
        
guard_1 = q(2) - ground_height;

[prog,V,coefv]= prog.newFreePoly(monomials(x,1:degree));

prog = prog.withSOS(V);

const_deg = degree;
sig = {};
coefsig = {};

% Vdot <= 0 in both modes
sos_vdot = -diff(V,x)*f_free;
[prog, sos_vdot, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_vdot, guard_1, x, const_deg);

sos_vdot_impact = -diff(V,x)*f_impact;
[prog, sos_vdot_impact, sig{end+1}, coefsig{end+1}] = spotless_add_eq_sprocedure(prog, sos_vdot_impact, guard_1, x, const_deg);
[prog, sos_vdot_impact, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_vdot_impact, -v(2), x, const_deg);

prog = prog.withSOS(sos_vdot);
prog = prog.withSOS(sos_vdot_impact);

options = spotprog.defaultOptions;
options.verbose = 1;
options.do_fr = false;
sol = prog.minimize(0,@spot_mosek,options);