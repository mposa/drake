megaclear
sos_option = 1;

switch sos_option
  case 1
    sos_fun = @spot_mosek_sos;
    withSOS_fun = @withSOS;
  case 2
    sos_fun = @spot_mosek_sdsos;
    withSOS_fun = @withSDSOS;
  case 3
    sos_fun = @spot_mosek_dsos;
%     sos_fun = @spot_gurobi_dsos;
    withSOS_fun = @withDSOS;
end

V_degree = 4;

g = 9.81;

prog = spotsosprog();

%% Add indeterminate variables
q = msspoly('q',4);
qd = msspoly('v',4);
lx = msspoly('lx',2);
lz = msspoly('lz',2);
lzsq = [1;1];

u = msspoly('u',1);

x = q(1);
z = q(2);
pitch = q(3);
theta = q(4);

xd = qd(1);
zd = qd(2);
pitchd = qd(3);
thetad = qd(4);

v_vars = [q(2:4);qd];
x_vars = v_vars;

prog = prog.withIndeterminate(q(2:4));
prog = prog.withIndeterminate(qd);
prog = prog.withIndeterminate(lx);

%% Dynamics
% [H,C,B,phi,phidot,psi,J,J_f,K,S,U] = torsoEOM_mss(q,qd,s_vec,c_vec);

% [phi_poly,phidot_poly,psi_poly,f_free_poly,f_impact_poly,E_poly,q,qd] = torsoPolyEOM(taylor_deg);
load torso_eom_poly

K = [10 1];
% u = -K*[s_th;thetad];
% U = U + K(1)*(1-c_th);
phi_poly = clean(phi_poly);
f_free_poly = clean(f_free_poly);
f_impact_poly = clean(f_impact_poly);

f_free_poly = subs(f_free_poly,u, -K*[q(4);qd(4)]);

phi = phi_poly;
psi = psi_poly;
phidot = phidot_poly;

%% Lyapunov function
[prog,V,coefv] = prog.newFreePoly(monomials(v_vars,1:4));
% V = E_poly + .5*K(1)*q(4)^2;

% [prog, equil_eqn] = prog.withEqs(subs(V,[z;s;c;s_th;c_th;qd],[0;0;1;0;1;0;0;0;0]));

Vdot_free = diff(V,[q;qd])*[qd;f_free_poly];
Vdot_impact_1 = diff(V,qd)*subs(f_impact_poly,[lx;lz],[lx(1);0;1;0]);
Vdot_impact_2 = diff(V,qd)*subs(f_impact_poly,[lx;lz],[0;lx(2);0;1]);


%% SOS functions
% (1) -Vdot_free(x) >= 0 for x admissable in B_o
% (2) -Vdot_impact_1(x,l) >= 0 for (x,l) admissable and x in B_o
% (3) -Vdot_impact_2(x,l) >= 0 for (x,l) admissable and x in B_o
% (4) V(x) >= 0 for x admissable and in B_o
% (5) V(x) - 1 >= 0 for x on bdry(B_o)
% (6) 1 - V(x) >= 0 for x in B_i

sos_1 = -Vdot_free;
sos_2 = -Vdot_impact_1;
sos_3 = -Vdot_impact_2;
sos_4 = V + .1;
sos_5 = V - 1;
% sos_5 = .5*vm*qd'*H*qd - 1;
% sos_5 = 10*U - .01;
sos_6 = 1 - V;


prog_bkp = prog;
%% Add in constraints
prog = prog_bkp;


const_deg = 4;
sig = {};
coefsig = {};


% Ball constraints
ball_vec = x_vars;

rho_i = .01;
rho_o = .05;

% Ao2 = Ao;
Ao2 = zeros(7)*z;
Ao2(1,1) = 100;
Ao2(2,2) = 5;
Ao2(3,3) = 2;

% Ao2(4:7,4:7) = .25*diff(diff(E_poly,qd)',qd);
Ao2(4:7,4:7) = eye(4);

h_Bo2 = rho_o - ball_vec'*Ao2*ball_vec;

h_Bi = rho_i - ball_vec'*Ao2*ball_vec;


doSOS = [1 1 1 1 1 1];

if doSOS(1)
  % non-penetration (admissability of x)
  [prog, sos_1, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_1, phi(1), [x_vars], const_deg, sos_option);
  [prog, sos_1, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_1, phi(2), [x_vars], const_deg, sos_option);
  [prog, sos_1, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_1, h_Bo2, [x_vars], const_deg, sos_option);
  prog = withSOS_fun(prog,sos_1);
end

if doSOS(2)
  [prog, sos_2, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_2, phi(2), [x_vars;lx(1)], const_deg, sos_option);
  % Contact constraints (admissability of lambda)
  [prog, sos_2, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_2, -phidot(1), [x_vars;lx(1)], const_deg, sos_option);
  [prog, sos_2, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_2, -lx(1)*psi(1), [x_vars;lx(1)], const_deg, sos_option);
  [prog, sos_2, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_2, lzsq(1)^2 - lx(1)^2, [x_vars;lx(1)], const_deg, sos_option);
  [prog, sos_2, sig{end+1}, coefsig{end+1}] = spotless_add_eq_sprocedure(prog, sos_2, phi(1), [x_vars;lx(1)], const_deg);
  [prog, sos_2, sig{end+1}, coefsig{end+1}] = spotless_add_eq_sprocedure(prog, sos_2, (lzsq(1)^2 - lx(1)^2)*psi(1), [x_vars;lx(1)], const_deg);  %should this be psi^2?
  [prog, sos_2, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_2, h_Bo2, [x_vars;lx(1)], const_deg, sos_option);
  prog = withSOS_fun(prog,sos_2);
end

if doSOS(3)
  [prog, sos_3, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_3, phi(1), [x_vars;lx(2)], const_deg, sos_option);
  [prog, sos_3, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_3, -phidot(2), [x_vars;lx(2)], const_deg, sos_option);
  [prog, sos_3, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_3, -lx(2)*psi(2), [x_vars;lx(2)], const_deg, sos_option);
  [prog, sos_3, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_3, lzsq(2)^2 - lx(2)^2, [x_vars;lx(2)], const_deg, sos_option);
  [prog, sos_3, sig{end+1}, coefsig{end+1}] = spotless_add_eq_sprocedure(prog, sos_3, phi(2), [x_vars;lx(2)], const_deg);
  [prog, sos_3, sig{end+1}, coefsig{end+1}] = spotless_add_eq_sprocedure(prog, sos_3, (lzsq(2)^2 - lx(2)^2)*psi(2), [x_vars;lx(2)], const_deg);  %should this be psi^2?
  [prog, sos_3, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_3, h_Bo2, [x_vars;lx(2)], const_deg, sos_option);
  prog = withSOS_fun(prog,sos_3);
end

if doSOS(4)
  [prog, sos_4, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_4, phi(1), x_vars, const_deg, sos_option);
  [prog, sos_4, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_4, phi(2), x_vars, const_deg, sos_option);
  [prog, sos_4, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_4, h_Bo2, x_vars, const_deg, sos_option);
  prog = withSOS_fun(prog,sos_4);
end

if doSOS(5)
  [prog, sos_5, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_5, phi(1), x_vars, const_deg, sos_option);
  [prog, sos_5, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_5, phi(2), x_vars, const_deg, sos_option);
  [prog, sos_5, sig{end+1}, coefsig{end+1}] = spotless_add_eq_sprocedure(prog, sos_5, h_Bo2, x_vars, const_deg);
  prog = withSOS_fun(prog,sos_5);
end

if doSOS(6)
  [prog, sos_6, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_6, phi(1), x_vars, const_deg, sos_option);
  [prog, sos_6, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_6, phi(2), x_vars, const_deg, sos_option);
  [prog, sos_6, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_6, h_Bi, x_vars, const_deg, sos_option);
  prog = withSOS_fun(prog,sos_6);
end


%% Solve program


options = spotprog.defaultOptions;
options.verbose = 1;
options.verbose = 1;

sol = prog.minimize(0,sos_fun,options);

Vsol = sol.eval(V);
