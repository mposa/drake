megaclear
iter = 2;
sos_option = 1;

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
    %     sos_fun = @spot_gurobi_dsos;
    withSOS_fun = @withDSOS;
    PSD_fun = @newDD;
    
end

% degree = 4;
b_degree = 1;
u_degree = 4;

change_variables = true;

% Ao = 1000*eye(9);


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


x0 = double(subs(x_vars,[q;qd;s_vec;c_vec],[zeros(12,1);ones(4,1)]));

options = spotprog.defaultOptions;
options.verbose = 1;
options.verbose = 1;
options.trig.enable = true;
options.trig.sin = [s;s_th];
options.trig.cos = [c;c_th];
options.clean_primal = true;


%% Dynamics
[H,C,B,phi,phidot,psi,J,J_f,K,S,U] = torsoEOM_mss(q,qd,s_vec,c_vec);

% do variable change
if change_variables
  T = [1 0 0 0;0 1 0 0; 0 0 1 0; 0 0 1 1];
  H = inv(T)'*H*inv(T);
  C = inv(T)'*C;
  B = inv(T)'*B;
  H=prog.trigExprReduction(reshape(subs(H(:),[s_th;c_th],[s_th*c - c_th*s; c_th*c+s_th*s]),4,[]),[s;s_th],[c;c_th]);
  C=prog.trigExprReduction(reshape(subs(C(:),[s_th;c_th;thetad],[s_th*c - c_th*s; c_th*c+s_th*s;thetad - pitchd]),4,[]),[s;s_th],[c;c_th]);
  U=prog.trigExprReduction(subs(U,[s_th;c_th],[s_th*c - c_th*s; c_th*c+s_th*s]),[s;s_th],[c;c_th]);
  
  K = [10 1];
  u = -K(1)*(s_th*c - c_th*s) - K(2)*(thetad - pitchd);  %TODO: ADD POTENTIAL ENERGY HERE
  U = U + K(1)*(1-s_th*s-c_th*c);
  
%   u = -K*[s_th;thetad];
%   U = U + K(1)*(1-c_th);
else
  K = [10 1];
  u = -K*[s_th;thetad];
  U = U + K(1)*(1-c_th);
end

H = clean(H);
C = clean(C);

%% Lyapunov function
if iter > 0
  load(sprintf('iter_%d',iter-1))
end
%   load(sprintf('iter_%d',iter-2))

if even(iter)% && 0
  [prog,b,coefb] = prog.newFreePoly(monomials([z;s;c;s_th;c_th],0:b_degree),4);
  [prog,Uq,coefu] = prog.newFreePoly(monomials([z;s;c;s_th;c_th],0:u_degree));
  [prog,vm]=prog.newFree(1);
elseif iter==1% || 1,
  b = bsol;
  Uq = Usol;
  vm = vmsol;
else
  b = bsol;
  Uq = Usol;
  vm = vmsol;
end


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

% rho_i = .4;
rho_o = 1;

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
Ao2(6:9,6:9) = .25*double(subs(H,x_vars,x0));
% Ao2(6:9,6:9) = .25*H;
%Ao2(6:9,6:9) = diag([3;3;.5;.25]);

% Ao2(6,6) = 0;
% Ao2(7,7) = 0;
% Ao2(8,8) = 0;
% Ao2(9,9) = 0;



cost_option = 3;

if iter==0,
  cost = 0;
  rho = .1;
  Ai = Ao2;
else
  switch cost_option
    case 1 % fix rho, diagonal Ai
      rho = .1;
      [prog,Ai_diag] = prog.newPos(9);
      [prog,AHmult] = prog.newPos(1);
      cost = sum(Ai_diag);
      %   Ai_diag(1) = 100*Ai_diag(1);
      Ai = diag(Ai_diag);
      Ai(6:9,6:9) = Ai(6:9,6:9) + AHmult*H;
    case 2 % minimize rho, fix trace, diagonal Ai
      [prog,rho] = prog.newFree(1);
      cost = -rho;
      [prog,Ai_diag] = prog.newPos(9);
      Ai_diag(1) = 100*Ai_diag(1);
      prog = prog.withEqs(sum(Ai_diag) - 6);
    case 3 % fix rho, PSD Ai
      [prog, Ai] = PSD_fun(prog,9);
      cost = trace(Ai);
      Ai(1,1) = Ai(1,1)*100;
      rho = .1;
  end
end
Ao = Ao2;
rho_i = rho;
rho_o = 1;

h_Bo2 = rho_o - ball_vec'*Ao2*ball_vec;


% if iter==0
%   cost=0;
%   rho_i = .3;
%   Ai = Ao2;
% else 
%   [prog,Ai_diag] = prog.newPos(9);
%   Ai = Ao2*diag(Ai_diag);
%   cost = Ai_diag(1) + .1*sum(Ai_diag);
% %   [prog,rho] = prog.newFree(1);
% %   cost = rho;
%   
% %   Ai = Ao2;
% %   cost = 0;
%   rho_i = .4;
% end

h_Bi = rho_i - ball_vec'*Ai*ball_vec; %worked with .01 and E, but failed sdsos


% changed hbo to hbo2
% ugh, been using the wrong hbo all along, for V>=1 and vdot <= 0
% FAILED K=20, A/5, bi=.05,bo=.15


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
if iter==0
  sos_6 = 1 - V;
else
  sos_6 = -h_Bi*(1 + z^2 + qd'*qd + 2 - c - c_th);
end

use_additional_eqs = true;

% if ~even(iter)
%   doSOS = [0 0 0 0 0 1];
% else
doSOS = [1 1 1 1 1 1];
% doSOS = [0 0 0 0 0 1];  %worked 1,2,3,4,5
% end

if doSOS(1)
  % non-penetration (admissability of x)
  [prog, sos_1, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_1, phi(1), [x_vars], const_deg, sos_option, options);
  [prog, sos_1, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_1, phi(2), [x_vars], const_deg, sos_option, options);
  [prog, sos_1, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_1, h_Bo2, [x_vars], const_deg, sos_option, options);
  if use_additional_eqs
    prog = prog.withEqs(subs(sig{end},x_vars,x0));
    prog = prog.withEqs(subs(diff(sig{end},qd(1)),x_vars,x0));
    prog = prog.withEqs(subs(diff(sig{end},qd(2)),x_vars,x0));
    prog = prog.withEqs(subs(diff(sig{end},qd(3)),x_vars,x0));
    prog = prog.withEqs(subs(diff(sig{end},qd(4)),x_vars,x0));
    prog = prog.withEqs(subs(diff(sig{end},s)*c + diff(sig{end},c)*-s,x_vars,x0));
    prog = prog.withEqs(subs(diff(sig{end},s_th)*c_th + diff(sig{end},c_th)*-s_th,x_vars,x0));
  end
  prog = withSOS_fun(prog,sos_1);
end

if doSOS(2)
  [prog, sos_2, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_2, phi(2), [x_vars;lx(1)], const_deg, sos_option, options);
  % Contact constraints (admissability of lambda)
  [prog, sos_2, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_2, -phidot(1), [x_vars;lx(1)], const_deg, sos_option, options);
  [prog, sos_2, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_2, -lx(1)*psi(1), [x_vars;lx(1)], const_deg, sos_option, options);
  [prog, sos_2, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_2, lzsq(1)^2 - lx(1)^2, [x_vars;lx(1)], const_deg, sos_option, options);
  [prog, sos_2, sig{end+1}, coefsig{end+1}] = spotless_add_eq_sprocedure(prog, sos_2, phi(1), [x_vars;lx(1)], const_deg);
  [prog, sos_2, sig{end+1}, coefsig{end+1}] = spotless_add_eq_sprocedure(prog, sos_2, (lzsq(1)^2 - lx(1)^2)*psi(1), [x_vars;lx(1)], const_deg);  %should this be psi^2?
  [prog, sos_2, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_2, h_Bo2, [x_vars;lx(1)], const_deg, sos_option, options);
  if use_additional_eqs
    prog = prog.withEqs(subs(sig{end},[x_vars;lx(1)],[x0;0]));
    prog = prog.withEqs(subs(diff(sig{end},qd(1)),[x_vars;lx(1)],[x0;0]));
    prog = prog.withEqs(subs(diff(sig{end},qd(2)),[x_vars;lx(1)],[x0;0]));
    prog = prog.withEqs(subs(diff(sig{end},qd(3)),[x_vars;lx(1)],[x0;0]));
    prog = prog.withEqs(subs(diff(sig{end},qd(4)),[x_vars;lx(1)],[x0;0]));
    prog = prog.withEqs(subs(diff(sig{end},s)*c + diff(sig{end},c)*-s,[x_vars;lx(1)],[x0;0]));
    prog = prog.withEqs(subs(diff(sig{end},s_th)*c_th + diff(sig{end},c_th)*-s_th,[x_vars;lx(1)],[x0;0]));
  end
  prog = withSOS_fun(prog,sos_2);
end

if doSOS(3)
  [prog, sos_3, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_3, phi(1), [x_vars;lx(2)], const_deg, sos_option, options);
  [prog, sos_3, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_3, -phidot(2), [x_vars;lx(2)], const_deg, sos_option, options);
  [prog, sos_3, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_3, -lx(2)*psi(2), [x_vars;lx(2)], const_deg, sos_option, options);
  [prog, sos_3, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_3, lzsq(2)^2 - lx(2)^2, [x_vars;lx(2)], const_deg, sos_option, options);
  [prog, sos_3, sig{end+1}, coefsig{end+1}] = spotless_add_eq_sprocedure(prog, sos_3, phi(2), [x_vars;lx(2)], const_deg);
  [prog, sos_3, sig{end+1}, coefsig{end+1}] = spotless_add_eq_sprocedure(prog, sos_3, (lzsq(2)^2 - lx(2)^2)*psi(2), [x_vars;lx(2)], const_deg);  %should this be psi^2?
  [prog, sos_3, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_3, h_Bo2, [x_vars;lx(2)], const_deg, sos_option, options);
  
  if use_additional_eqs
    prog = prog.withEqs(subs(sig{end},[x_vars;lx(2)],[x0;0]));
    prog = prog.withEqs(subs(diff(sig{end},qd(1)),[x_vars;lx(2)],[x0;0]));
    prog = prog.withEqs(subs(diff(sig{end},qd(2)),[x_vars;lx(2)],[x0;0]));
    prog = prog.withEqs(subs(diff(sig{end},qd(3)),[x_vars;lx(2)],[x0;0]));
    prog = prog.withEqs(subs(diff(sig{end},qd(4)),[x_vars;lx(2)],[x0;0]));
    prog = prog.withEqs(subs(diff(sig{end},s)*c + diff(sig{end},c)*-s,[x_vars;lx(2)],[x0;0]));
    prog = prog.withEqs(subs(diff(sig{end},s_th)*c_th + diff(sig{end},c_th)*-s_th,[x_vars;lx(2)],[x0;0]));
  end
  
  prog = withSOS_fun(prog,sos_3);
end

if doSOS(4)
  [prog, sos_4, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_4, phi(1), x_vars, const_deg, sos_option, options);
  [prog, sos_4, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_4, phi(2), x_vars, const_deg, sos_option, options);
  [prog, sos_4, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_4, h_Bo2, x_vars, const_deg, sos_option, options);
  if use_additional_eqs
    prog = prog.withEqs(subs(sig{end},x_vars,x0));
    prog = prog.withEqs(subs(diff(sig{end},qd(1)),x_vars,x0));
    prog = prog.withEqs(subs(diff(sig{end},qd(2)),x_vars,x0));
    prog = prog.withEqs(subs(diff(sig{end},qd(3)),x_vars,x0));
    prog = prog.withEqs(subs(diff(sig{end},qd(4)),x_vars,x0));
    prog = prog.withEqs(subs(diff(sig{end},s)*c + diff(sig{end},c)*-s,x_vars,x0));
    prog = prog.withEqs(subs(diff(sig{end},s_th)*c_th + diff(sig{end},c_th)*-s_th,x_vars,x0));
  end
  prog = withSOS_fun(prog,sos_4);
end

if doSOS(5)
  [prog, sos_5, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_5, phi(1), x_vars, const_deg, sos_option, options);
  [prog, sos_5, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_5, phi(2), x_vars, const_deg, sos_option, options);
  [prog, sos_5, sig{end+1}, coefsig{end+1}] = spotless_add_eq_sprocedure(prog, sos_5, h_Bo2, x_vars, const_deg);
  prog = withSOS_fun(prog,sos_5);
end

if doSOS(6)
  [prog, sos_6, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_6, phi(1), x_vars, const_deg, sos_option, options);
  [prog, sos_6, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_6, phi(2), x_vars, const_deg, sos_option, options);
  if iter==0,
    [prog, sos_6, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_6, h_Bi, x_vars, const_deg, sos_option, options);
  elseif ~even(iter)
    [prog, sos_6, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_6, V - 1, x_vars, const_deg, sos_option, options);
    sos6_mult = sig{end};
  else
    sos_6 = sos_6 - (V-1)*sos6_mult;
  end
  prog = withSOS_fun(prog,sos_6);
end



% Solve program




sol = prog.minimize(cost,sos_fun,options);

sqrt(1/double(sol.eval(Ai(1,1)/double(sol.eval(rho_i)))))

Vsol = sol.eval(V);
bsol = sol.eval(b);
Usol = sol.eval(Uq);
vmsol = sol.eval(vm);
AI = sol.eval(Ai);
AO = double(sol.eval(Ao));
R = double(sol.eval(rho_i));
D = double(sol.eval(rho_o));
if even(iter)
  save iter_0 Vsol bsol vmsol Usol Ao2 AI R
  save(strcat(sprintf('iter_%d',iter)),'Vsol','bsol','vmsol','Usol','R','AO','AI','D')
else
  sos6_mult = sol.eval(sos6_mult);
  save(strcat(sprintf('iter_%d',iter)),'Vsol','bsol','vmsol','Usol','R','AO','AI','D','sos6_mult')
end

