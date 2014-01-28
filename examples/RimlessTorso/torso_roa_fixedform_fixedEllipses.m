clear all
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

% degree = 4;
b_degree = 1;
u_degree = 4;

Ao = 1000*eye(9);

rho = .1;
Ai = 100*eye(9);
Ai = Ai/rho;

g = 9.81;

prog = spotsosprog();

%% Add indeterminate variables
q = msspoly('q',4);
qd = msspoly('qd',4);
s_vec = msspoly('s',4);
c_vec = msspoly('c',4);
lx = msspoly('lx',2);
lzsq = [1;1];

x = q(1);
z = q(2);

s = s_vec(3);
s_th = s_vec(4);

c = c_vec(3);
c_th = c_vec(4);

xd = qd(1);
zd = qd(2);
pitchd = qd(3);
thetad = qd(4);

v_vars = [q(2);s_vec(3:4);c_vec(3:4);qd];
x_vars = v_vars;

prog = prog.withIndeterminate(q(2));
prog = prog.withIndeterminate(s_vec(3:4));
prog = prog.withIndeterminate(c_vec(3:4));
prog = prog.withIndeterminate(qd);
prog = prog.withIndeterminate(lx);




%% Dynamics
[H,C,B,phi,phidot,psi,J,J_f,K,S,U] = torsoEOM_mss(q,qd,s_vec,c_vec);

K = [10 1];
u = -K*[s_th;thetad];
U = U + K(1)*(1-c_th);

%% Lyapunov function
[prog,b,coefb] = prog.newFreePoly(monomials([z;s;c;s_th;c_th],0:b_degree),4);
[prog,Uq,coefu] = prog.newFreePoly(monomials([z;s;c;s_th;c_th],0:u_degree));
[prog,vm]=prog.newFree(1);
% b = 0;
% Uq = U;
% load torso_data_new0_125
% vm = vmsol;
% Uq = Usol;
% b = bsol;
V = .5*vm*qd'*H*qd + b'*H*qd + Uq;
% V = .5*vm*qd'*H*qd + vm*U;
E = .5*qd'*H*qd + U;
% vm = 1; V = .5*vm*qd'*H*qd + vm*U; b = 0;


[prog, equil_eqn] = prog.withEqs(subs(V,[z;s;c;s_th;c_th;qd],[0;0;1;0;1;0;0;0;0]));

Vdot_free = diff(V,[x;z;s;c;s_th;c_th])*[qd(1:2);c*pitchd;-s*pitchd;c_th*thetad;-s_th*thetad] + (vm*qd + b)'*(-C + B*u);
Vdot_impact_1 = (vm*qd + b)'*(J(1,:)'*lzsq(1) + J_f(1,:)'*lx(1));
Vdot_impact_2 = (vm*qd + b)'*(J(2,:)'*lzsq(2) + J_f(2,:)'*lx(2));


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
sos_4 = V;
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
ball_vec = [z;s;1-c;s_th;1-c_th;qd];
% h_Bo = 1 - ball_vec'*Ao*ball_vec;
% h_Bi = 1 - ball_vec'*Ai*ball_vec;

rho_i = .4;
rho_o = 1.5;

% Ao2 = Ao;
Ao2 = zeros(9)*z;
% Ao2(1,1) = 10;
% Ao2(2,2) = .5;
% Ao2(3,3) = .5;
% Ao2(4,4) = 2.5/5;
% Ao2(5,5) = 2.5/5;

Ao2(1,1) = 100;
Ao2(2,2) = 5;
Ao2(3,3) = 5;
Ao2(4,4) = 2;
Ao2(5,5) = 2;

Ao2(6:9,6:9) = .25*H;
% Ao2(6,6) = 0;
% Ao2(7,7) = 0;
% Ao2(8,8) = 0;
% Ao2(9,9) = 0;
h_Bo2 = rho_o - ball_vec'*Ao2*ball_vec;
%START COMMENTS, 1/27
% Set Ao2 to 100,5,5,2,2,.25H
% K = 10
%  worked .4/1.5
%  failed .6/1.5
%  failed .6/1
%  failed .6/2
%  failed .5/2
%  failed .5/1.5

% Set Ao2 to 100,5,5,1,1,.25H
% K = 10
%  failed .5/1.5
%  worked .4/1.5

% IGNORE COMMENTS BELOW THIS LINE

% h_Bo2 = .1 - (2-c-c_th) - z^2 - .05*.5*qd'*H*qd;

% divided th by 2, worked with .03, .045, unk .06, .055

% set K=10, worked .05, unk .07, .06

% set K=10, divided A_th by 5, changed Bo2 to .15 (from .1)
% worked .05, .1, .125, failed .2, .15

% as above, but K=50, failed .15, .1 trying and .05
% as above, divided A_th by 2 only, trying .1, failed

%K=20, A/5, worked .05

%with SDSOS, worked .02, unk .04, .03
%searched for V with .01,.02,.03,.05,.06, .0625 worked unk with .07,
%.065,.6375
h_Bi = rho_i - ball_vec'*Ao2*ball_vec; %worked with .01 and E, but failed sdsos


% changed hbo to hbo2
% ugh, been using the wrong hbo all along, for V>=1 and vdot <= 0
% FAILED K=20, A/5, bi=.05,bo=.15

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


% Solve program


options = spotprog.defaultOptions;
options.verbose = 1;
options.verbose = 1;
options.trig.enable = true;
options.trig.sin = [s;s_th];
options.trig.cos = [c;c_th];
sol = prog.minimize(0,sos_fun,options);

Vsol = sol.eval(V);
bsol = sol.eval(b);
Usol = sol.eval(Uq);
vmsol = sol.eval(vm);