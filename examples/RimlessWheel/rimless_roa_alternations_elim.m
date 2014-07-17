% clear all
% from flex_mult
% eliminating a variable (z) for contact conditions.
%
% merged in phi calculations from rimless_roa_moment_inf_time
megaclear
for iter=31:50,
display(sprintf('Starting iter %d',iter))
sos_option = 1;
do_backoff = 0;
do_clean = 0;
do_eig = 0;

switch sos_option
  case 1
%     sos_fun = @spot_sedumi;
%     sos_fun = @spot_mosek_sos;
%     sos_fun = @spot_sdpnal;
    withSOS_fun = @withSOS;
    PSD_fun = @newPSD;
  case 2
%     sos_fun = @spot_mosek_sdsos;
    withSOS_fun = @withSDSOS;
    PSD_fun = @newSDD;
  case 3
%     sos_fun = @spot_mosek_dsos;
%     sos_fun = @spot_gurobi_dsos;7
    withSOS_fun = @withDSOS;
    PSD_fun = @newDD;
end
sos_fun= @spot_mosek;

% degree = 4;
V_degree = 4;

g = 9.81;

prog = spotsosprog();

file_prefix = 'reg8_sq_lowmu_elim_iter_%d';

%% Add indeterminate variables
q = msspoly('q',3);
qd = msspoly('qd',3);
s_vec = msspoly('s',3);
c_vec = msspoly('c',3);
lx = msspoly('lx',2);
lzsq = [1;1];

x = q(1);
z = q(2);

s = s_vec(3);
c = c_vec(3);

xd = qd(1);
zd = qd(2);
pitchd = qd(3);

v_vars = [q(2);s_vec(3);c_vec(3);qd];
x_vars = v_vars;

prog = prog.withIndeterminate(q(2));
prog = prog.withIndeterminate(s_vec(3));
prog = prog.withIndeterminate(c_vec(3));
prog = prog.withIndeterminate(qd);
prog = prog.withIndeterminate(lx);

x0 = [0;0;1;0;0;0];

options = spotprog.defaultOptions;
options.verbose = 1;
options.trig.enable = true;
options.trig.sin = s;
options.trig.cos = c;
options.scale_monomials = true;
options.backoff = false;
% options.clean_primal = true;
% options.regularize = true;
options.clean_primal = false;
options.regularize = true;
options.regularize_eps = 1e-8;

%% Dynamics
% [H,C,B,phi,phidot,psi,J,J_f,K,S,U] = torsoEOM_mss(q,qd,s_vec,c_vec);
% H = clean(H);
% C = clean(C);
% K = [10 1];
% u = -K*[s_th;thetad];
% U = U + K(1)*(1-c_th);


cpi8 = cos(pi/8);
spi8 = sin(pi/8);
rt2 = sqrt(2);


% f_impact=[0, 0;
%   0, 0;
%   0, 0;
%   0, 0;
%   lx(1), lx(2);
%   lzsq(1)^2, lzsq(2)^2;
%   4*lzsq(1)^2*(spi8*c + s*cpi8)- 4*(-lx(1))*(c*cpi8 - s*spi8), 4*lzsq(2)^2*(-spi8*c + s*cpi8)- 4*(-lx(2))*(c*cpi8 + s*spi8)];
% 
% phi = [z + s*spi8 - c*cpi8 + cpi8; z - s*spi8 - c*cpi8 + cpi8];
% 
% phidot = [zd + c*pitchd*spi8 + s*pitchd*cpi8; zd - c*pitchd*spi8 + s*pitchd*cpi8];
% 
% psi = [xd + c*pitchd*cpi8 - s*pitchd*spi8; xd + c*pitchd*cpi8 + s*pitchd*spi8];

f_free=[xd;
  zd;
  c*pitchd;
  -s*pitchd;
  0;
  -g
  0];

ground_angle = 0;

fk_foot{1} = [x + s*cpi8 - c*spi8;z + s*spi8 - c*cpi8 + cpi8];
fk_foot{2} = [x + s*cpi8 + c*spi8;z - s*spi8 - c*cpi8 + cpi8];

T_contact = [cos(ground_angle) -sin(ground_angle); sin(ground_angle) cos(ground_angle)];

phi = [T_contact(2,:)*fk_foot{1}; T_contact(2,:)*fk_foot{2}];

J_fk{1} = diff(fk_foot{1},[x;z;s;c]);
J_fk{2} = diff(fk_foot{2},[x;z;s;c]);

vk_foot{1} = T_contact*J_fk{1}*f_free(1:4);
vk_foot{2} = T_contact*J_fk{2}*f_free(1:4);

phidot = [vk_foot{1}(2); vk_foot{2}(2)];

psi = [vk_foot{1}(1); vk_foot{2}(1)];

f_impact = [zeros(4,2); diag([1;1;4])*[1 0 0 0;0 1 0 0; 0 0 c -s]*[J_fk{1}'*[lx(1);lzsq(1)^2] J_fk{2}'*[lx(2);lzsq(2)^2]]]; 

E = .5*xd^2 + .5*zd^2 + 1/8*pitchd^2 + g*z;

%% Lyapunov function
if iter > 0
  load(datapath(sprintf(file_prefix,iter-1)))
end

if iter==0,
  [prog,V,coefv] = prog.newFreePoly(monomials(v_vars,0:V_degree));
  [prog, equil_eqn] = prog.withEqs(subs(V,[z;s;c;qd],[0;0;1;0;0;0]));
  
  prog = prog.withEqs(subs(diff(V,qd(1)),x_vars,x0));
  prog = prog.withEqs(subs(diff(V,qd(2)),x_vars,x0));
  prog = prog.withEqs(subs(diff(V,qd(3)),x_vars,x0));
  prog = prog.withEqs(subs(diff(V,s)*c + diff(V,c)*-s,x_vars,x0));
elseif ~even(iter),
  V = Vsol;
else
%   load iter_1
  [prog,V,coefv] = prog.newFreePoly(monomials(v_vars,0:V_degree));
  [prog, equil_eqn] = prog.withEqs(subs(V,[z;s;c;qd],[0;0;1;0;0;0]));
  
  prog = prog.withEqs(subs(diff(V,qd(1)),x_vars,x0));
  prog = prog.withEqs(subs(diff(V,qd(2)),x_vars,x0));
  prog = prog.withEqs(subs(diff(V,qd(3)),x_vars,x0));
  prog = prog.withEqs(subs(diff(V,s)*c + diff(V,c)*-s,x_vars,x0));
end

% V = E;

Vdot_free = diff(V,[x;z;s;c;xd;zd;pitchd])*f_free;
Vdot_impact_1 = diff(V,[x;z;s;c;xd;zd;pitchd])*f_impact(:,1);
Vdot_impact_2 = diff(V,[x;z;s;c;xd;zd;pitchd])*f_impact(:,2);


%% Ball constraints
ball_vec = [z;s;1-c;qd];
% Ao2 = zeros(6);
% Ao2(1,1) = 100;
% Ao2(2,2) = 2;
% Ao2(3,3) = 2;
% Ao2(4,4) = .5;
% Ao2(5,5) = .5;
% Ao2(6,6) = .1;

Ao2 = zeros(6);
Ao2(1,1) = 60;
Ao2(2,2) = 2;
Ao2(3,3) = 2;
Ao2(4,4) = .5;
Ao2(5,5) = .5;
Ao2(6,6) = .1;

Ao2 = Ao2*1.5;

if iter==0,
  cost = 0;
  rho_i = .2;
  Ai = Ao2;
  Ao = Ao2;
  rho_o = 1;
  rho_Vo = 1;
%   Ai = diag([50;6;10;10;10;10]);
%   Ai = diag([50;50;50;50;50;50]);
% Ai = diag([10;.25;.005;.1;.1;.05]);
% rho = .03;
%   rho = .1;
elseif ~even(iter)
  [prog, Ao] = PSD_fun(prog,6);
  [prog, rho_Vo] = prog.newFree(1);
  cost = -rho_Vo;
  Ai = AI;
  rho_i = R;
  rho_o = 1;
else
  cost_option = 4;
  
  switch cost_option
    case 1 % fix rho, diagonal Ai
      rho_i = .1;
      [prog,Ai_diag] = prog.newPos(6);
      cost = sum(Ai_diag);
      %   Ai_diag(1) = 100*Ai_diag(1);
      Ai = diag(Ai_diag);
    case 2 % minimize rho, fix trace, diagonal Ai
      [prog,rho_i] = prog.newFree(1);
      cost = -rho_i;
      [prog,Ai_diag] = prog.newPos(6);
      prog = prog.withEqs(sum(Ai_diag) - 6);
      Ai_diag(1) = 100*Ai_diag(1);
    case 3 % fix rho, PSD Ai
      [prog, Ai] = PSD_fun(prog,6);
      cost = trace(Ai);
%      Ai(1,1) = Ai(1,1)*100;
      rho_i = .1;
    case 4 % minimize rho, fix trace, psd Ai
      [prog,rho_i] = prog.newFree(1);
      cost = -rho_i;
%       Ai11mult = .1;
      [prog, Ai] = PSD_fun(prog,6);
%       Ai(1,1) = 100*Ai(1,1);
      prog = prog.withEqs(trace(Ai) - 100);
%       prog = prog.withEqs(Ai(1,1) + trace(Ai(2:end,2:end)) - 100);
%       Ai(1,1) = Ai11mult*Ai(1,1);
%       Ai_diag(1) = 100*Ai_diag(1);
  end
  Ao = AO;
  rho_o = D;
  rho_Vo = 1;
end

h_Bo2 = rho_o - ball_vec'*Ao*ball_vec;

h_Bi = rho_i - ball_vec'*Ai*ball_vec;

if even(iter)
  rho_V = 1;
else
  rho_V = 1.0;
end

prog_bkp = prog;
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
sos_5 = (V - rho_Vo);%*(1 + z^2 + qd'*qd + 1 - c)^2;

if iter==0,
  sos_6 = (rho_V - V);
else
%   sos_6 = -h_Bi*tmp_mult;
  sos_6 = -h_Bi*(1 + z^2 + qd'*qd + 1 - c)^2;
end


%% Add in constraints
prog = prog_bkp;

mu = .2;

const_deg = 4;
sig = {};
coefsig = {};


doSOS = [1 1 1 1 1 1];
% cost = 0;
% doSOS = [1 1 1 1 0 0 0];

if ~even(iter) && 0
  doSOS = [0 0 0 0 0 1 0];
end

if doSOS(1)
  % non-penetration (admissability of x)
  [prog, sos_1, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_1, phi(1), [x_vars], const_deg, sos_option, options);
  [prog, sos_1, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_1, phi(2), [x_vars], const_deg, sos_option, options);
  
  if even(iter) || iter==0
    [prog, sos_1, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_1, h_Bo2, [x_vars], const_deg, sos_option, options);
    sos1_mult = sig{end};
    
    prog = prog.withEqs(subs(sig{end},x_vars,x0));
    prog = prog.withEqs(subs(diff(sig{end},qd(1)),x_vars,x0));
    prog = prog.withEqs(subs(diff(sig{end},qd(2)),x_vars,x0));
    prog = prog.withEqs(subs(diff(sig{end},qd(3)),x_vars,x0));
    prog = prog.withEqs(subs(diff(sig{end},s)*c + diff(sig{end},c)*-s,x_vars,x0));
  else
    sos_1 = sos_1 - h_Bo2*sos1_mult;
  end
  
  prog = withSOS_fun(prog,sos_1);
end

if doSOS(2)
  [prog, sos_2, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_2, phi(2), [x_vars(2:end);lx(1)], const_deg, sos_option, options);
  % Contact constraints (admissability of lambda)
  [prog, sos_2, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_2, -phidot(1), [x_vars(2:end);lx(1)], const_deg, sos_option, options);
  [prog, sos_2, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_2, -lx(1)*psi(1), [x_vars(2:end);lx(1)], const_deg-2, sos_option, options);
  [prog, sos_2, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_2, mu^2*lzsq(1)^2 - lx(1)^2, [x_vars(2:end);lx(1)], const_deg, sos_option, options);
%   [prog, sos_2, sig{end+1}, coefsig{end+1}] = spotless_add_eq_sprocedure(prog, sos_2, phi(1), [x_vars;lx(1)], const_deg+1);
  [prog, sos_2, sig{end+1}, coefsig{end+1}] = spotless_add_eq_sprocedure(prog, sos_2, (mu^2*lzsq(1)^2 - lx(1)^2)*psi(1), [x_vars(2:end);lx(1)], const_deg-2);  %should this be psi^2?
  
  if even(iter) || iter==0
    [prog, sos_2, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_2, h_Bo2, [x_vars(2:end);lx(1)], const_deg, sos_option, options);
    sos2_mult = sig{end};
    
    prog = prog.withEqs(subs(sig{end},[x_vars;lx(1)],[x0;0]));
    prog = prog.withEqs(subs(diff(sig{end},qd(1)),[x_vars;lx(1)],[x0;0]));
    prog = prog.withEqs(subs(diff(sig{end},qd(2)),[x_vars;lx(1)],[x0;0]));
    prog = prog.withEqs(subs(diff(sig{end},qd(3)),[x_vars;lx(1)],[x0;0]));
    prog = prog.withEqs(subs(diff(sig{end},s)*c + diff(sig{end},c)*-s,[x_vars;lx(1)],[x0;0]));
  else
    sos_2 = sos_2 - h_Bo2*sos2_mult;
  end
  sos_2 = subs(sos_2,z,z-phi(1));
  prog = withSOS_fun(prog,sos_2);
end

if doSOS(3)
  [prog, sos_3, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_3, phi(1), [x_vars(2:end);lx(2)], const_deg, sos_option, options);
  [prog, sos_3, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_3, -phidot(2), [x_vars(2:end);lx(2)], const_deg, sos_option, options);
  [prog, sos_3, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_3, -lx(2)*psi(2), [x_vars(2:end);lx(2)], const_deg-2, sos_option, options);
  [prog, sos_3, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_3, mu^2*lzsq(2)^2 - lx(2)^2, [x_vars(2:end);lx(2)], const_deg, sos_option, options);
%   [prog, sos_3, sig{end+1}, coefsig{end+1}] = spotless_add_eq_sprocedure(prog, sos_3, phi(2), [x_vars;lx(2)], const_deg+1);
  [prog, sos_3, sig{end+1}, coefsig{end+1}] = spotless_add_eq_sprocedure(prog, sos_3, (mu^2*lzsq(2)^2 - lx(2)^2)*psi(2), [x_vars(2:end);lx(2)], const_deg-2);  %should this be psi^2?
  
  if even(iter) || iter==0
    [prog, sos_3, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_3, h_Bo2, [x_vars(2:end);lx(2)], const_deg, sos_option, options);
    sos3_mult = sig{end};
    
    prog = prog.withEqs(subs(sig{end},[x_vars;lx(2)],[x0;0]));
    prog = prog.withEqs(subs(diff(sig{end},qd(1)),[x_vars;lx(2)],[x0;0]));
    prog = prog.withEqs(subs(diff(sig{end},qd(2)),[x_vars;lx(2)],[x0;0]));
    prog = prog.withEqs(subs(diff(sig{end},qd(3)),[x_vars;lx(2)],[x0;0]));
    prog = prog.withEqs(subs(diff(sig{end},s)*c + diff(sig{end},c)*-s,[x_vars;lx(2)],[x0;0]));
    
  else
    sos_3 = sos_3 - h_Bo2*sos3_mult;
  end
  sos_3 = subs(sos_3,z,z-phi(2));
  prog = withSOS_fun(prog,sos_3);
  
end

if doSOS(4)
  [prog, sos_4, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_4, phi(1), x_vars, const_deg, sos_option, options);
  [prog, sos_4, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_4, phi(2), x_vars, const_deg, sos_option, options);
  
  if even(iter) || iter==0
    [prog, sos_4, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_4, h_Bo2, x_vars, const_deg, sos_option, options);
    sos4_mult = sig{end};
    
    prog = prog.withEqs(subs(sig{end},x_vars,x0));
    prog = prog.withEqs(subs(diff(sig{end},qd(1)),x_vars,x0));
    prog = prog.withEqs(subs(diff(sig{end},qd(2)),x_vars,x0));
    prog = prog.withEqs(subs(diff(sig{end},qd(3)),x_vars,x0));
    prog = prog.withEqs(subs(diff(sig{end},s)*c + diff(V,c)*-s,x_vars,x0));
    
  else
    sos_4 = sos_4 - h_Bo2*sos4_mult;
  end

  prog = withSOS_fun(prog,sos_4);
end

if doSOS(5)
  [prog, sos_5, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_5, phi(1), x_vars, const_deg, sos_option, options);
  [prog, sos_5, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_5, phi(2), x_vars, const_deg, sos_option, options);
  
  if even(iter) || iter==0 
    [prog, sos_5, sig{end+1}, coefsig{end+1}] = spotless_add_eq_sprocedure(prog, sos_5, h_Bo2, x_vars, const_deg);
    sos5_mult = sig{end};
  else
    sos_5 = sos_5 - h_Bo2*sos5_mult;
  end
  
  prog = withSOS_fun(prog,sos_5);
end

if doSOS(6)
  [prog, sos_6, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_6, phi(1), x_vars, const_deg, sos_option, options);
  [prog, sos_6, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_6, phi(2), x_vars, const_deg, sos_option, options);
  %     [prog, sos_6, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_6, phi(1)*phi(2), x_vars, const_deg, sos_option, options);
  
  if iter==0,
    [prog, sos_6, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_6, h_Bi, x_vars, const_deg, sos_option, options);
    sos6_mult = sig{end};
  elseif ~even(iter),
    %     [prog, sos_6, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_6, (V - 1)*phi(1), x_vars, const_deg, sos_option, options);
    %     [prog, sos_6, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_6, (V - 1)*phi(2), x_vars, const_deg, sos_option, options);
    [prog, sos_6, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_6, V - rho_V, x_vars, const_deg-2, sos_option, options);
    sos6_mult = sig{end};
  else
    %     sos_6 = sos_6 - (V-1)*phi(1)*sos6_mult_3;
    %     sos_6 = sos_6 - (V-1)*phi(2)*sos6_mult_2;
    sos_6 = sos_6 - (V-rho_V)*sos6_mult;
  end
  
  prog = withSOS_fun(prog,sos_6);
end

% if doSOS(7)
%   if iter==0
%     sos_7 = subs(1-V,[s;c;qd],[0;1;zeros(3,1)]);
%     [prog, sos_7, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_7, phi(1), x_vars, const_deg, sos_option, options);
%     [prog, sos_7, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_7, phi(2), x_vars, const_deg, sos_option, options);
%     [prog, sos_7, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_7, .15^2-z^2, z, const_deg, sos_option, options);
%     prog = withSOS_fun(prog,sos_7);
%   end
% end


% Solve program

% options.do_fr = true;
%options.solveroptions.MSK_IPAR_INTPNT_BASIS = 'MSK_BI_NEVER';
% options.regularize = false;
sol = prog.minimize(cost,sos_fun,options);

display(sprintf('Finished iter %d with cost=%f',iter,double(sol.eval(cost))))

sqrt(1/double(sol.eval(Ai(1,1)/double(sol.eval(rho_i)))))

do_backoff = do_backoff && ~isequal(cost,0);

% for i=1:length(prog.sosExpr)
%   [v,d] = eig(double(sol.eval(sol.gramMatrices{i})));
%   if any(diag(d) <= 1e-6)
%     options.gram_reduction{i} = v(diag(d) > 1e-6,:);
%   else
%     options.gram_reduction{i} = eye(length(d));
%   end
% end

if do_eig
  [prog,gamma] = prog.newPos(1);
  cost = cost + gamma;
  for i=1:length(prog.sosExpr),
    prog.sosExpr(i) = prog.sosExpr(i) + 1e-6*gamma*sol.gramMonomials{i}'*sol.gramMonomials{i};
  end
  sol = prog.minimize(cost,sos_fun,options);
end

if do_clean && do_backoff
  display(sprintf('Eliminating %d variables by cleaning',sum(abs(double(sol.eval(prog.freeVar))) < 1e-6)))
  prog = prog.withEqs(prog.freeVar(abs(double(sol.eval(prog.freeVar))) < 1e-6));
  costval = double(sol.eval(cost));
  prog = prog.withPos(costval + .1*abs(costval) - cost);
  cost = 0;
  sol = prog.minimize(cost,sos_fun,options);
else
  
  if do_clean
    display(sprintf('Eliminating %d variables by cleaning',sum(abs(double(sol.eval(prog.freeVar))) < 1e-6)))
    elim_vars_ind = find(abs(double(sol.eval(prog.freeVar))) < 1e-6);
    elim_vars = prog.freeVar(elim_vars_ind);
    keep_ind = setdiff(1:length(prog.freeVar),elim_vars_ind); 
    prog.A = prog.A(:,keep_ind);
    prog.freeVar = prog.freeVar(keep_ind);
%     prog = prog.withEqs(prog.freeVar(abs(double(sol.eval(prog.freeVar))) < 1e-6));

    prog.sosExpr = subs(prog.sosExpr,elim_vars,zeros(length(elim_vars_ind),1));

    sol = prog.minimize(cost,sos_fun,options);
  end
  
  if do_backoff
    costval = double(sol.eval(cost));

    prog = prog.withPos(costval + .1*abs(costval) - cost);
%     cost = 0;
    
%     [prog,gamma] = prog.newFree(1);
%     cost = gamma;
%     for i=1:length(prog.sosExpr),
%       prog.sosExpr(i) = prog.sosExpr(i) + gamma*sol.gramMonomials{i}'*sol.gramMonomials{i};
%     end
    sol = prog.minimize(0,sos_fun,options);
  end
end

Vsol = sol.eval(V)/double(sol.eval(rho_Vo));
R = double(sol.eval(rho_i));
AI = double(sol.eval(Ai));
AO = double(sol.eval(Ao));
D = double(sol.eval(rho_o));
if ~even(iter),
%   sos6_mult_2 = sol.eval(sig{end-1});
%   sos6_mult_3 = sol.eval(sig{end-2});
sos6_mult = sol.eval(sos6_mult)*double(sol.eval(rho_Vo));
  save(datapath(strcat(sprintf(file_prefix,iter))),'Vsol','R','AO','AI','D','sos6_mult')
%   save iter_1 Vsol sos6_mult sos6_mult_2 sos6_mult_3 R Ao2 AI
else

  sos1_mult = sol.eval(sos1_mult);
  sos2_mult = sol.eval(sos2_mult);
  sos3_mult = sol.eval(sos3_mult);
  sos4_mult = sol.eval(sos4_mult);
  sos5_mult = sol.eval(sos5_mult);
  save(datapath(strcat(sprintf(file_prefix,iter))),'Vsol','sos1_mult','sos2_mult','sos3_mult','sos4_mult','sos5_mult','R','AO','AI','D')  
end
end
