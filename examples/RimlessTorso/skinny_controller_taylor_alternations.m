%very good results with z_scale = .1, SOS, degree 4, controller 3, reduced
%degree in eqn 5, 5e-6,1e-6, 1e-5

megaclear
for iter = 48:80,
display(sprintf('Starting iter %d',iter))
do_backoff = true;
sos_option = 1;
z_scale = .1;

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

V_degree = 4;

change_variables = true;

% Ao = 1000*eye(9);


g = 9.81;

prog = spotsosprog();

% file_prefix = 'lowdeg_skinny_cubic_controller_taylor_iter_%d';
file_prefix = 'zscale_skinny_cubic_controller_taylor_iter_%d';
% file_prefix = 'sdsos_zscale_skinny_cubic_controller_taylor_iter_%d';
if iter > 0
  load(datapath(sprintf(file_prefix,iter-1)))
end

controller_degree = 3;

%% Add indeterminate variables
q = msspoly('q',4);
qd = msspoly('qd',4);
u = msspoly('u',1);
lx = msspoly('lx',2);
lzsq = [1;1];

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


x0 = zeros(7,1);

options = spotprog.defaultOptions;
options.verbose = 1;
options.verbose = 1;
options.trig.enable = false;
options.clean_primal = true;
options.scale_monomials = true;
options.regularize = true;

if iter == 0,
  options.regularize_eps = 5e-6;
elseif even(iter)
  options.regularize_eps = 1e-6; %cubic, have been using 1e-6. linear, iter 2, needed 1e-5
else
  options.regularize_eps = 1e-5; 
end


%% Taylor Dynamics
taylor_degree = 3;
taylor_vars = [q;qd;u;lx];
load skinny_taylor_eom
f_impact_1 = clean(getmsspoly(f_impact_1,taylor_vars,taylor_degree));
f_impact_2 = clean(getmsspoly(f_impact_2,taylor_vars,taylor_degree));
f_free = clean(getmsspoly(f_free,taylor_vars,taylor_degree));
phi = clean(getmsspoly(phi,taylor_vars,taylor_degree-1));
phidot = clean(getmsspoly(phidot,taylor_vars,taylor_degree-1));
psi = clean(getmsspoly(psi,taylor_vars,taylor_degree));
E = clean(getmsspoly(E,taylor_vars,taylor_degree));

% invH = clean(getmsspoly(invH,taylor_vars,taylor_degree));
K = [10 1];
E = E + .5*K(1)*theta^2;


f_impact_1(2) = subs(f_impact_1(2)/z_scale,[z;zd],[z;zd]*z_scale);
f_impact_2(2) = subs(f_impact_2(2)/z_scale,[z;zd],[z;zd]*z_scale);
f_free(2) = subs(f_free(2)/z_scale,[z;zd],[z;zd]*z_scale);
phi = subs(phi,[z;zd],[z;zd]*z_scale);
phidot = subs(phidot,[z;zd],[z;zd]*z_scale);
psi = subs(psi,[z;zd],[z;zd]*z_scale);
E = subs(E,[z;zd],[z;zd]*z_scale);

prog_bkp = prog;
%% Add in constraints
prog = prog_bkp;


const_deg = 4;
sig = {};
coefsig = {};


% Ball constraints
ball_vec = v_vars;

% Ao2 = Ao;
Ao2 = zeros(7);

% Ao2(1,1) = 200;
% Ao2(2,2) = 5;
Ao2(1,1) = 1000*z_scale^2;
Ao2(2,2) = 25;
Ao2(3,3) = 5;
Ao2(4:7,4:7) = .5*double(subs(diff(diff(E,qd)',qd),x_vars,x0));



cost_option = 4;
odd_cost_option = 1;

if iter==0,
  cost = 0;
  if sos_option == 1
    rho_i = .02;
    rho_o = 1;
  else
    rho_i = .01;
    rho_o = 1;
  end
  Ai = Ao2;
  Ao = Ao2;
  
  rho_Vo = 1;
  
  controller = -K*[theta; thetad];
elseif ~even(iter)
  [prog,controller,coefu] = prog.newFreePoly(monomials(v_vars,1:controller_degree));  
  
  if iter == 1,
    Ao = Ao2;% on iter 1, numerical issues with arbitrary Ao
  else
    [prog, Ao] = PSD_fun(prog,7);
  end  
  
  Ai = AI;
  rho_o = 1;
  
  switch odd_cost_option
    case 1
      [prog, rho_Vo] = prog.newFree(1);
      rho_i = R;
      cost = -rho_Vo;
    case 2
      [prog, rho_Vo] = prog.newFree(1);
      [prog,rho_i] = prog.newFree(1);
      cost = -rho_Vo - rho_i;
    case 3
      rho_Vo = 1;
      [prog,rho_i] = prog.newFree(1);
      cost = -rho_i;
  end
else
  controller = controllersol;
  cost_option = 4;
  
  switch cost_option
    case 1 % fix rho, diagonal Ai
      rho_i = .1;
      [prog,Ai_diag] = prog.newPos(9);
      cost = sum(Ai_diag);
      %   Ai_diag(1) = 100*Ai_diag(1);
      Ai = diag(Ai_diag);
    case 2 % minimize rho, fix trace, diagonal Ai
      [prog,rho_i] = prog.newFree(1);
      cost = -rho_i;
      [prog,Ai_diag] = prog.newPos(9);
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
      [prog, Ai] = PSD_fun(prog,7);
      prog = prog.withEqs(trace(Ai) - 7);
      TA = diag([10*z_scale;ones(6,1)]);
      Ai = TA'*Ai*TA;
  end
  Ao = AO;
  rho_o = D;
  rho_Vo = 1;
end

f_free = subs(f_free,u,controller);

h_Bo2 = rho_o - ball_vec'*Ao*ball_vec;




h_Bi = rho_i - ball_vec'*Ai*ball_vec; 

%% Lyapunov function
mu = .5;

if even(iter)% && 0
  [prog,V,coefv] = prog.newFreePoly(monomials(v_vars,1:V_degree));  
  
%   prog = prog.withEqs(subs(diff(V,qd(1)),x_vars,x0));
%   prog = prog.withEqs(subs(diff(V,qd(2)),x_vars,x0));
%   prog = prog.withEqs(subs(diff(V,qd(3)),x_vars,x0));
%   prog = prog.withEqs(subs(diff(V,qd(4)),x_vars,x0));
%   prog = prog.withEqs(subs(diff(V,theta),x_vars,x0));
%   prog = prog.withEqs(subs(diff(V,pitch),x_vars,x0));
else
  V = Vsol;
end

Vdot_free = diff(V,q)*qd + diff(V,qd)*f_free;
Vdot_impact_1 = diff(V,qd)*f_impact_1;
Vdot_impact_2 = diff(V,qd)*f_impact_2;


% SOS functions
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
sos_5 = V - rho_Vo;
if iter==0
  sos_6 = 1 - V;
else
  sos_6 = -h_Bi*(1 + x_vars'*x_vars)^2;  %failed without this ^2? 
end

use_additional_eqs = false;

% doSOS = [1 1 1 1 1 1];
doSOS = [1 1 1 0 1 1];  % not using V >= 0 constraint for taylor expansion

% doSOS = [1 1 1 0 0 1];
% cost = 0;

if doSOS(1)
  % non-penetration (admissability of x)
  [prog, sos_1, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_1, phi(1), [x_vars], const_deg + controller_degree - 1, sos_option, options);
  [prog, sos_1, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_1, phi(2), [x_vars], const_deg + controller_degree - 1, sos_option, options);
  
  if even(iter)
    [prog, sos_1, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_1, h_Bo2, [x_vars], const_deg + controller_degree - 1, sos_option, options);
    sos1_mult = sig{end};
    if use_additional_eqs
      prog = prog.withEqs(subs(diff(sig{end},qd(1)),x_vars,x0));
      prog = prog.withEqs(subs(diff(sig{end},qd(2)),x_vars,x0));
      prog = prog.withEqs(subs(diff(sig{end},qd(3)),x_vars,x0));
      prog = prog.withEqs(subs(diff(sig{end},qd(4)),x_vars,x0));
      prog = prog.withEqs(subs(diff(sig{end},theta),x_vars,x0));
      prog = prog.withEqs(subs(diff(sig{end},pitch),x_vars,x0));
    end
  else
    sos_1 = sos_1 - h_Bo2*sos1_mult;
  end
  prog = withSOS_fun(prog,sos_1);
end

if doSOS(2)
  [prog, sos_2, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_2, phi(2), [x_vars;lx(1)], const_deg, sos_option, options);
  % Contact constraints (admissability of lambda)
  [prog, sos_2, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_2, -phidot(1), [x_vars;lx(1)], const_deg, sos_option, options);
  [prog, sos_2, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_2, -lx(1)*psi(1), [x_vars;lx(1)], const_deg-2, sos_option, options);
  [prog, sos_2, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_2, mu^2*lzsq(1)^2 - lx(1)^2, [x_vars;lx(1)], const_deg, sos_option, options);
  [prog, sos_2, sig{end+1}, coefsig{end+1}] = spotless_add_eq_sprocedure(prog, sos_2, phi(1), [x_vars;lx(1)], const_deg);
  [prog, sos_2, sig{end+1}, coefsig{end+1}] = spotless_add_eq_sprocedure(prog, sos_2, (mu^2*lzsq(1)^2 - lx(1)^2)*psi(1), [x_vars;lx(1)], const_deg-3);  %should this be psi^2?
  
  if even(iter)
    [prog, sos_2, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_2, h_Bo2, [x_vars;lx(1)], const_deg, sos_option, options);
    sos2_mult = sig{end};
    if use_additional_eqs
      prog = prog.withEqs(subs(diff(sig{end},qd(1)),[x_vars;lx(1)],[x0;0]));
      prog = prog.withEqs(subs(diff(sig{end},qd(2)),[x_vars;lx(1)],[x0;0]));
      prog = prog.withEqs(subs(diff(sig{end},qd(3)),[x_vars;lx(1)],[x0;0]));
      prog = prog.withEqs(subs(diff(sig{end},qd(4)),[x_vars;lx(1)],[x0;0]));
      prog = prog.withEqs(subs(diff(sig{end},theta),[x_vars;lx(1)],[x0;0]));
      prog = prog.withEqs(subs(diff(sig{end},pitch),[x_vars;lx(1)],[x0;0]));
    end
  else
    sos_2 = sos_2 - h_Bo2*sos2_mult;
  end
  prog = withSOS_fun(prog,sos_2);
end

if doSOS(3)
  [prog, sos_3, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_3, phi(1), [x_vars;lx(2)], const_deg, sos_option, options);
  [prog, sos_3, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_3, -phidot(2), [x_vars;lx(2)], const_deg, sos_option, options);
  [prog, sos_3, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_3, -lx(2)*psi(2), [x_vars;lx(2)], const_deg-2, sos_option, options);
  [prog, sos_3, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_3, mu^2*lzsq(2)^2 - lx(2)^2, [x_vars;lx(2)], const_deg, sos_option, options);
  [prog, sos_3, sig{end+1}, coefsig{end+1}] = spotless_add_eq_sprocedure(prog, sos_3, phi(2), [x_vars;lx(2)], const_deg);
  [prog, sos_3, sig{end+1}, coefsig{end+1}] = spotless_add_eq_sprocedure(prog, sos_3, (mu^2*lzsq(2)^2 - lx(2)^2)*psi(2), [x_vars;lx(2)], const_deg-3);  %should this be psi^2?
  
  if even(iter)
    [prog, sos_3, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_3, h_Bo2, [x_vars;lx(2)], const_deg, sos_option, options);
    sos3_mult = sig{end};
    
    if use_additional_eqs
      prog = prog.withEqs(subs(diff(sig{end},qd(1)),[x_vars;lx(2)],[x0;0]));
      prog = prog.withEqs(subs(diff(sig{end},qd(2)),[x_vars;lx(2)],[x0;0]));
      prog = prog.withEqs(subs(diff(sig{end},qd(3)),[x_vars;lx(2)],[x0;0]));
      prog = prog.withEqs(subs(diff(sig{end},qd(4)),[x_vars;lx(2)],[x0;0]));
      prog = prog.withEqs(subs(diff(sig{end},theta),[x_vars;lx(2)],[x0;0]));
      prog = prog.withEqs(subs(diff(sig{end},pitch),[x_vars;lx(2)],[x0;0]));
    end
  else
    sos_3 = sos_3 - h_Bo2*sos3_mult;
  end
  prog = withSOS_fun(prog,sos_3);
end

if doSOS(4)
  [prog, sos_4, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_4, phi(1), x_vars, const_deg, sos_option, options);
  [prog, sos_4, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_4, phi(2), x_vars, const_deg, sos_option, options);
  
  
  if even(iter)
    [prog, sos_4, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_4, h_Bo2, x_vars, const_deg, sos_option, options);
    sos4_mult = sig{end};
    if use_additional_eqs
      prog = prog.withEqs(subs(diff(sig{end},qd(1)),x_vars,x0));
      prog = prog.withEqs(subs(diff(sig{end},qd(2)),x_vars,x0));
      prog = prog.withEqs(subs(diff(sig{end},qd(3)),x_vars,x0));
      prog = prog.withEqs(subs(diff(sig{end},qd(4)),x_vars,x0));
      prog = prog.withEqs(subs(diff(sig{end},theta),x_vars,x0));
      prog = prog.withEqs(subs(diff(sig{end},pitch),x_vars,x0));
    end
  else
    sos_4 = sos_4 - h_Bo2*sos4_mult;
  end
  prog = withSOS_fun(prog,sos_4);
else
  sos4_mult = 0;
end

if doSOS(5)
  [prog, sos_5, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_5, phi(1), x_vars, const_deg-2, sos_option, options);
  [prog, sos_5, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_5, phi(2), x_vars, const_deg-2, sos_option, options);
  
  if even(iter)
    [prog, sos_5, sig{end+1}, coefsig{end+1}] = spotless_add_eq_sprocedure(prog, sos_5, h_Bo2, x_vars, const_deg-2);
    sos5_mult = sig{end};
  else
    sos_5 = sos_5 - h_Bo2*sos5_mult;
  end
   
  prog = withSOS_fun(prog,sos_5);
end

if doSOS(6)
  [prog, sos_6, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_6, phi(1), x_vars, const_deg, sos_option, options);
  [prog, sos_6, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_6, phi(2), x_vars, const_deg, sos_option, options);
  if iter==0,
    [prog, sos_6, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_6, h_Bi, x_vars, const_deg, sos_option, options);
  elseif ~even(iter)
    [prog, sos_6, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_6, V - 1, x_vars, const_deg-2, sos_option, options);
    sos6_mult = sig{end};
  else
    sos_6 = sos_6 - (V-1)*sos6_mult;
  end
  prog = withSOS_fun(prog,sos_6);
end

% Solve program


sol = prog.minimize(cost,sos_fun,options);

if do_backoff && ~isnumeric(cost)
  display('Backing off and solving again...')
  sol_beforebackoff = sol;
  prog_beforebackoff = prog;
  costval = double(sol.eval(cost));
  prog = prog.withPos(costval + .01*abs(costval) - cost);
  
  [prog,bo_eps] = prog.newFree(1);
  switch sos_option
    case 1
      for i=1:length(prog.sosExpr),
        prog.sosExpr(i) = prog.sosExpr(i) + bo_eps*sol.gramMonomials{i}'*sol.gramMonomials{i};
      end
    case 2
      for i=1:length(prog.sdsosExpr),
        prog.sdsosExpr(i) = prog.sdsosExpr(i) + bo_eps*sol.gramMonomials{i}'*sol.gramMonomials{i};
      end
  end
  sol = prog.minimize(bo_eps,sos_fun,options);
else
  bo_eps = 0;
end

sqrt(1/double(sol.eval(Ai(1,1)/double(sol.eval(rho_i)))))*z_scale

Vsol = sol.eval(V)/double(sol.eval(rho_Vo));
controllersol = sol.eval(controller);

AI = sol.eval(Ai);
AO = double(sol.eval(Ao));
R = double(sol.eval(rho_i));
D = double(sol.eval(rho_o));
bo_eps = double(sol.eval(bo_eps));
if even(iter)
  sos1_mult = sol.eval(sos1_mult);
  sos2_mult = sol.eval(sos2_mult);
  sos3_mult = sol.eval(sos3_mult);
  sos4_mult = sol.eval(sos4_mult);
  sos5_mult = sol.eval(sos5_mult);
  
  save(datapath(strcat(sprintf(file_prefix,iter))),'Vsol','R','AO','AI','D','sos1_mult','sos2_mult','sos3_mult','sos4_mult','sos5_mult','controllersol','bo_eps')
else
  sos6_mult = sol.eval(sos6_mult)*double(sol.eval(rho_Vo));
  save(datapath(strcat(sprintf(file_prefix,iter))),'Vsol','R','AO','AI','D','sos6_mult','controllersol','bo_eps')
end

end