% two masses connected by two springs (m1 to the sky, m2 to m1) with
% contact between m2 and the ground, via a hybrid method
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

x_1 = [q;v];
x_2 = [q(1);v(1)];
% x = [q;v]
f_free = [-1/m1*k1*q(1) - 1/m1*b1*v(1) - 1/m1*k2*(q(1) - q(2)) - 1/m1*b2*(v(1) - v(2));
          1/m2*k2*(q(1) - q(2))+1/m2*b2*(v(1) - v(2));v];
        
% x= [q(1);v(1)]
% where q(2) = ground_height, v(2) = 0;
f_ground = [-1/m1*k1*q(1) - 1/m1*k2*(q(1) - ground_height);v(1)];
        
guard_1 = q(2) - ground_height;
guard_2 = ground_height - q(1);

jump_1 = [q(1);v(1)];
jump_2 = [q(1);ground_height;v(1);0];

[prog,V1,coefv1]= prog.newFreePoly(monomials(x_1,1:degree));
[prog,V2,coefv2]= prog.newFreePoly(monomials(x_2,0:degree));

prog = prog.withSOS(V1);
prog = prog.withSOS(V2);

const_deg = degree;
sig = {};
coefsig = {};

% Vdot <= 0 in both modes
sos_vdot_1 = -diff(V1,x_1)*f_free;
[prog, sos_vdot_1, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_vdot_1, guard_1, x_1, const_deg);

sos_vdot_2 = -diff(V2,x_2)*f_ground;
[prog, sos_vdot_2, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_vdot_2, guard_2, x_2, const_deg);

prog = prog.withSOS(sos_vdot_1);
prog = prog.withSOS(sos_vdot_2);

% V decreases across the jumpes

sos_jump_1 = V1 - subs(V2,x_2,jump_1);
[prog, sos_jump_1, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_jump_1, -guard_1, x_1, const_deg);

sos_jump_2 = V2 - subs(V1,x_1,jump_2);
[prog, sos_jump_2, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_jump_2, -guard_2, x_2, const_deg);

prog = prog.withSOS(sos_jump_1);
prog = prog.withSOS(sos_jump_2);

options = spotprog.defaultOptions;
options.verbose = 1;
options.do_fr = false;
sol = prog.minimize(0,@spot_mosek,options);