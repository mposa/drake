megaclear
iter = 0;
sos_option = 1;
do_backoff = 0;
do_clean = 0;

switch sos_option
  case 1
    sos_fun = @spot_mosek_sos;
    withSOS_fun = @withSOS;
    PSD_fun = @newPSD;
  case 2
    sos_fun = @spot_mosek_sdsos;
    withSOS_fun = @withSDSOS;
    PSD_fun = @newSDD;
  case 3
    sos_fun = @spot_mosek_dsos;
    withSOS_fun = @withDSOS;
    PSD_fun = @newDD;
end

V_degree = 2;

g = 9.81;

prog = spotsosprog();

%% Add indeterminate variables
q = msspoly('q',9);
qd = msspoly('qd',9);

x = q(1:3:end);
z = q(2:3:end);
theta = q(3:3:end);

xd = qd(1:3:end);
zd = qd(2:3:end);
thetad = qd(3:3:end);

lx = msspoly('lx',3);
u = msspoly('u',6);
lzsq = [1;1;1];

v_vars = [q(1:2);q(4:end);qd];
x_vars = v_vars;

prog = prog.withIndeterminate([q;qd;lx]);

m = [1;1;1];
R = [.1;.1;.1];
I = 2/4*m.*R.^2;

f_free = [0;-g;0;0;-g;0;0;-g;0] + [zeros(3,1);u];

% xnom =  [0;.14;0;-.05;0;0;.05;0;0;zeros(9,1)];
% xnom =  [0;.12;0;-.09;0;0;.09;0;0;zeros(9,1)];
xnom =  [0;.2;0;-.12;0;0;.12;0;0;zeros(9,1)];

K = zeros(6,18);
K(1:3,4:6) = diag([100;100;10]);
K(1:3,13:15) = diag([30;30;3]);
K(4:6,7:9) = diag([100;100;10]);
K(4:6,16:18) = diag([30;30;3]);
% K = K*0;
u_control = -K*([q;qd] - xnom) + [0;1.5*g;0;0;1.5*g;0];
% u_control = u_control*0;
f_free = subs(f_free,u,u_control.*(1./[m(2);m(2);I(2);m(3);m(3);I(2)]));

r_perp_1_2 = [x(2) - x(1); z(2) - z(1)];
r_perp_2_3 = [x(3) - x(2); z(3) - z(2)];
r_perp_3_1 = [x(1) - x(3); z(1) - z(3)];

r_fric_1_2 = [-r_perp_1_2(2);r_perp_1_2(1)];
r_fric_2_3 = [-r_perp_2_3(2);r_perp_2_3(1)];
r_fric_3_1 = [-r_perp_3_1(2);r_perp_3_1(1)];

f_impact_1_2 = [(-r_perp_1_2 - lx(1)*r_fric_1_2)/m(1);...
  -R(1)/I(1)*lx(1);...
  (r_perp_1_2 + lx(1)*r_fric_1_2)/m(2);...
  -R(2)/I(2)*lx(1);...
  zeros(3,1)];

f_impact_2_3 = [zeros(3,1);...
  (-r_perp_2_3 - lx(2)*r_fric_2_3)/m(2);...
  -R(2)/I(2)*lx(2);...
  (r_perp_2_3 + lx(2)*r_fric_2_3)/m(3);...
  -R(3)/I(3)*lx(2)];

f_impact_3_1 = [
    (r_perp_3_1 + lx(3)*r_fric_3_1)/m(1);...
  -R(1)/I(1)*lx(3);...
  zeros(3,1);...
  (-r_perp_3_1 - lx(3)*r_fric_3_1)/m(3);...
  -R(3)/I(3)*lx(3)];

phi = 100*[r_perp_1_2'*r_perp_1_2 - (R(1) + R(2))^2;...
  r_perp_2_3'*r_perp_2_3 - (R(2) + R(3))^2;...
  r_perp_3_1'*r_perp_3_1 - (R(3) + R(1))^2];

phidot = diff(phi,q)*qd;

% psi = [(qd(4:5) - qd(1:2) + (R(1)*qd(3) + R(2)*qd(6))*r_fric_1_2)'*r_fric_1_2;...
%   (qd(7:8) - qd(4:5) + (R(2)*qd(6) + R(3)*qd(9))*r_fric_2_3)'*r_fric_2_3;...
%   (qd(1:2) - qd(7:8) + (R(3)*qd(9) + R(1)*qd(3))*r_fric_3_1)'*r_fric_3_1];

psi = 100*[(qd(4:5) - qd(1:2))'*r_fric_1_2 - R(1)*qd(3) - R(2)*qd(6);...
  (qd(7:8) - qd(4:5))'*r_fric_2_3 - R(2)*qd(6) - R(3)*qd(9);...
  (qd(1:2) - qd(7:8))'*r_fric_3_1 - R(3)*qd(9) - R(1)*qd(3)];

E = sum(m(1)*g.*z(1)) + sum(.5*m.*(xd.^2 + zd.^2) + .5*I.*thetad.^2) + sum(.5*K(:,1:9)*(q-xnom(1:9)).^2);
E = E - subs(E,x_vars,xnom([1 2 4:end]));

%% Lyapunov function
% V = 10*E;
if iter==0,
  [prog,V,coefv] = prog.newFreePoly(monomials(v_vars,1:V_degree));
elseif iter==1,
  load iter_0
  V = Vsol;
elseif ~even(iter),
  load iter_2
  V = Vsol;
else
  load iter_1
  [prog,V,coefv] = prog.newFreePoly(monomials(v_vars,1:V_degree));
end

Vdot_free = diff(V,q)*qd + diff(V,qd)*f_free;
Vdot_impact_1_2 = diff(V,qd)*f_impact_1_2;
Vdot_impact_2_3 = diff(V,qd)*f_impact_2_3;
Vdot_impact_3_1 = diff(V,qd)*f_impact_3_1;


%% Ball constraints
ball_vec = x_vars;
% Ao2 = 20*eye(18);

Ao2 = zeros(17);
Ao2(1:2,1:2) = diag([1000;1000]);
Ao2(3:end,3:end) = diag([  100.0000
  100.0000
   10.0000
  100.0000
  100.0000
   10.0000
    1.0000
    1.0000
    0.0050
    1.0000
    1.0000
    0.0050
    1.0000
    1.0000
    0.0050]);
% Ao2(3:8,3:8) = K(:,4:9);
% Ao2(9:17,9:17) = diag([4;4;.02;1;1;.0013;1;1;.0013]);

if iter==0,
  cost = 0;
  rho = .01;
  Ai = Ao2;
else
%   [prog,rho] = prog.newFree(1);
%   cost = -rho;

%   Ai = Ao2;
%   [prog, Ai] = PSD_fun(prog,6);

  [prog,Ai_diag] = prog.newPos(17);
  Ai = diag(Ai_diag);
  cost = trace(Ai);

  rho = .1;
end

rho_i = rho;
rho_o = .1;



h_Bo2 = rho_o - (ball_vec-xnom([1 2 4:end]))'*Ao2*(ball_vec-xnom([1 2 4:end]));

h_Bi = rho_i - (ball_vec-xnom([1 2 4:end]))'*Ai*(ball_vec-xnom([1 2 4:end]));

%% SOS functions
% (1) -Vdot_free(x) >= 0 for x admissable in B_o
% (2) -Vdot_impact_1_2(x,l) >= 0 for (x,l) admissable and x in B_o
% (3) -Vdot_impact_2_3(x,l) >= 0 for (x,l) admissable and x in B_o
% (4) -Vdot_impact_3_1(x,l) >= 0 for (x,l) admissable and x in B_o
% (5) V(x) - 1 >= 0 for x in bnd(B_o);
% (6) 1 - V(x) >= 0 for x in B_i

%1) try energy, pick a better x0

% [prog,d] = prog.newFree(1);
% cost = d;

sos_1 = -Vdot_free;
sos_2 = -Vdot_impact_1_2;
sos_3 = -Vdot_impact_2_3;
sos_4 = -Vdot_impact_3_1;
sos_5 = V - 1;

if iter==0,
  sos_6 = (1 - V);
else
  sos_6 = -h_Bi;%*(1 + z^2 + qd'*qd + 1 - c);
end

prog_bkp = prog;

%% Add in constraints
prog = prog_bkp;


const_deg = 2;
sig = {};
coefsig = {};


doSOS = [1 1 1 1 1 1];
% doSOS = [0 0 0 0 1 0];

if doSOS(1)
  % non-penetration (admissability of x)
  [prog, sos_1, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_1, phi(1), [x_vars], const_deg, sos_option);
  [prog, sos_1, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_1, phi(2), [x_vars], const_deg, sos_option);
  [prog, sos_1, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_1, phi(3), [x_vars], const_deg, sos_option);
  [prog, sos_1, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_1, h_Bo2, [x_vars], const_deg, sos_option);
  
  prog = withSOS_fun(prog,sos_1);
end

if doSOS(2)
  [prog, sos_2, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_2, phi(2), [x_vars;lx(1)], const_deg, sos_option);
  [prog, sos_2, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_2, phi(3), [x_vars;lx(1)], const_deg, sos_option);
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
  [prog, sos_3, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_3, phi(3), [x_vars;lx(2)], const_deg, sos_option);
  [prog, sos_3, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_3, -phidot(2), [x_vars;lx(2)], const_deg, sos_option);
  [prog, sos_3, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_3, -lx(2)*psi(2), [x_vars;lx(2)], const_deg, sos_option);
  [prog, sos_3, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_3, lzsq(2)^2 - lx(2)^2, [x_vars;lx(2)], const_deg, sos_option);
  [prog, sos_3, sig{end+1}, coefsig{end+1}] = spotless_add_eq_sprocedure(prog, sos_3, phi(2), [x_vars;lx(2)], const_deg);
  [prog, sos_3, sig{end+1}, coefsig{end+1}] = spotless_add_eq_sprocedure(prog, sos_3, (lzsq(2)^2 - lx(2)^2)*psi(2), [x_vars;lx(2)], const_deg);  %should this be psi^2?
  [prog, sos_3, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_3, h_Bo2, [x_vars;lx(2)], const_deg, sos_option);
  
  prog = withSOS_fun(prog,sos_3);
  
end

if doSOS(4)
  [prog, sos_4, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_4, phi(1), [x_vars;lx(3)], const_deg, sos_option);
  [prog, sos_4, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_4, phi(2), [x_vars;lx(3)], const_deg, sos_option);
  [prog, sos_4, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_4, -phidot(3), [x_vars;lx(3)], const_deg, sos_option);
  [prog, sos_4, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_4, -lx(3)*psi(3), [x_vars;lx(3)], const_deg, sos_option);
  [prog, sos_4, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_4, lzsq(3)^2 - lx(3)^2, [x_vars;lx(3)], const_deg, sos_option);
  [prog, sos_4, sig{end+1}, coefsig{end+1}] = spotless_add_eq_sprocedure(prog, sos_4, phi(3), [x_vars;lx(3)], const_deg);
  [prog, sos_4, sig{end+1}, coefsig{end+1}] = spotless_add_eq_sprocedure(prog, sos_4, (lzsq(3)^2 - lx(3)^2)*psi(3), [x_vars;lx(3)], const_deg);  %should this be psi^2?
  [prog, sos_4, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_4, h_Bo2, [x_vars;lx(3)], const_deg, sos_option);
  
  prog = withSOS_fun(prog,sos_3);
  
end

if doSOS(5)
  [prog, sos_5, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_5, phi(1), x_vars, const_deg, sos_option);
  [prog, sos_5, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_5, phi(2), x_vars, const_deg, sos_option);
  [prog, sos_5, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_5, phi(3), x_vars, const_deg, sos_option);
  [prog, sos_5, sig{end+1}, coefsig{end+1}] = spotless_add_eq_sprocedure(prog, sos_5, h_Bo2, x_vars, const_deg);
  prog = withSOS_fun(prog,sos_5);
end

if doSOS(6)
  [prog, sos_6, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_6, phi(1), x_vars, const_deg, sos_option);
  [prog, sos_6, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_6, phi(2), x_vars, const_deg, sos_option);
  [prog, sos_6, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_6, phi(3), x_vars, const_deg, sos_option);
  
  if iter==0,
    [prog, sos_6, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_6, h_Bi, x_vars, const_deg, sos_option);
  elseif ~even(iter),

    [prog, sos_6, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_6, V - 1, x_vars, const_deg, sos_option);
  else
    load iter_1
    sos_6 = sos_6 - (V-1)*sos6_mult;
  end
  
  prog = withSOS_fun(prog,sos_6);
end


% Solve program


options = spotprog.defaultOptions;
options.verbose = 1;
sol = prog.minimize(cost,sos_fun,options);


%%
Vsol = sol.eval(V);
R = double(sol.eval(rho));
AI = double(sol.eval(Ai));
if iter==0,
  save iter_0 Vsol R Ao2 AI
elseif ~even(iter),
  sos6_mult = sol.eval(sig{end});
%   sos6_mult_2 = sol.eval(sig{end-1});
%   sos6_mult_3 = sol.eval(sig{end-2});
  save iter_1 Vsol sos6_mult R Ao2 AI
%   save iter_1 Vsol sos6_mult sos6_mult_2 sos6_mult_3 R Ao2 AI
else
  save iter_2 Vsol R Ao2 AI
end
