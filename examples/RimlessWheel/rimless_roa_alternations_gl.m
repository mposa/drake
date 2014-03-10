% clear all
megaclear
for iter=5:20,
display(sprintf('Starting iter %d',iter))
sos_option = 1;
do_backoff = 0;
do_clean = 0;
do_eig = 0;

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

% degree = 4;
V_degree = 4;

g = 9.81;

prog = spotsosprog();

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
options.clean_primal = true;
options.scale_monomials = true;

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


f_impact=[0, 0;
  0, 0;
  0, 0;
  0, 0;
  lx(1), lx(2);
  lzsq(1)^2, lzsq(2)^2;
  4*lzsq(1)^2*(spi8*c + s*cpi8)- 4*(-lx(1))*(c*cpi8 - s*spi8), 4*lzsq(2)^2*(-spi8*c + s*cpi8)- 4*(-lx(2))*(c*cpi8 + s*spi8)];

phi = [z + s*spi8 - c*cpi8 + cpi8; z - s*spi8 - c*cpi8 + cpi8];

phidot = [zd + c*pitchd*spi8 + s*pitchd*cpi8; zd - c*pitchd*spi8 + s*pitchd*cpi8];

psi = [xd + c*pitchd*cpi8 - s*pitchd*spi8; xd + c*pitchd*cpi8 + s*pitchd*spi8];

f_free=[xd;
  zd;
  c*pitchd;
  -s*pitchd;
  0;
  -g
  0];

E = .5*xd^2 + .5*zd^2 + 1/8*pitchd^2 + g*z;

%% Lyapunov function
if iter > 0
  load(sprintf('gl_iter_%d',iter-1))
end

if ~even(iter),
  V = Vsol;
else
%   load iter_1
  [prog,V,coefv] = prog.newFreePoly(monomials(v_vars,0:V_degree));
  [prog, equil_eqn] = prog.withEqs(subs(V,[z;s;c;qd],[0;0;1;0;0;0]));
  [prog, equil_eqn2] = prog.withEqs(1 - subs(V,[z;s;c;qd],[.1;0;1;0;0;0]));
  
  prog = prog.withEqs(subs(diff(V,qd(1)),x_vars,x0));
  prog = prog.withEqs(subs(diff(V,qd(2)),x_vars,x0));
  prog = prog.withEqs(subs(diff(V,qd(3)),x_vars,x0));
  prog = prog.withEqs(subs(diff(V,s)*c + diff(V,c)*-s,x_vars,x0));
end

% V = E;

Vdot_free = diff(V,[x;z;s;c;xd;zd;pitchd])*f_free;
Vdot_impact_1 = diff(V,[x;z;s;c;xd;zd;pitchd])*f_impact(:,1);
Vdot_impact_2 = diff(V,[x;z;s;c;xd;zd;pitchd])*f_impact(:,2);


%% cost setup
if iter==0,
  sos1_mult = (z^2 + 1 - c + qd'*qd)^2;
  sos2_mult = (z^2 + 1 - c + qd'*qd)^2;
  sos3_mult = (z^2 + 1 - c + qd'*qd)^2;
end
if ~even(iter)
  rho = 1.5*R;
  cost = 0;
else
  [prog,rho] = prog.newFree(1);
  cost = -rho;
end

prog_bkp = prog;
%% SOS functions
% (1) -Vdot_free(x) >= 0 for x admissable and V < rho
% (2) -Vdot_impact_1(x,l) >= 0 for (x,l) admissable and V < rho
% (3) -Vdot_impact_2(x,l) >= 0 for (x,l) admissable and V < rho
% (4) V(x) >= 0 for x admissable


sos_1 = -Vdot_free;
sos_2 = -Vdot_impact_1;
sos_3 = -Vdot_impact_2;
sos_4 = V;

%% Add in constraints
prog = prog_bkp;

mu = 1;

const_deg = 4;
sig = {};
coefsig = {};


doSOS = [1 1 1 1];

if doSOS(1)
  % non-penetration (admissability of x)
  [prog, sos_1, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_1, phi(1), [x_vars], const_deg, sos_option, options);
  [prog, sos_1, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_1, phi(2), [x_vars], const_deg, sos_option, options);
  
  if ~even(iter)
    [prog, sos_1, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_1, rho - V, [x_vars], const_deg, sos_option, options);
    sos1_mult = sig{end};
    
    prog = prog.withEqs(subs(sig{end},x_vars,x0));
    prog = prog.withEqs(subs(diff(sig{end},qd(1)),x_vars,x0));
    prog = prog.withEqs(subs(diff(sig{end},qd(2)),x_vars,x0));
    prog = prog.withEqs(subs(diff(sig{end},qd(3)),x_vars,x0));
    prog = prog.withEqs(subs(diff(sig{end},s)*c + diff(sig{end},c)*-s,x_vars,x0));
  else
    sos_1 = sos_1 - (rho - V)*sos1_mult;
  end
  
  prog = withSOS_fun(prog,sos_1);
end

if doSOS(2)
  [prog, sos_2, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_2, phi(2), [x_vars;lx(1)], const_deg, sos_option, options);
  % Contact constraints (admissability of lambda)
  [prog, sos_2, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_2, -phidot(1), [x_vars;lx(1)], const_deg, sos_option, options);
  [prog, sos_2, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_2, -lx(1)*psi(1), [x_vars;lx(1)], const_deg, sos_option, options);
  [prog, sos_2, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_2, mu^2*lzsq(1)^2 - lx(1)^2, [x_vars;lx(1)], const_deg, sos_option, options);
  [prog, sos_2, sig{end+1}, coefsig{end+1}] = spotless_add_eq_sprocedure(prog, sos_2, phi(1), [x_vars;lx(1)], const_deg);
  [prog, sos_2, sig{end+1}, coefsig{end+1}] = spotless_add_eq_sprocedure(prog, sos_2, (mu^2*lzsq(1)^2 - lx(1)^2)*psi(1), [x_vars;lx(1)], const_deg);  %should this be psi^2?
  
  if ~even(iter)
    [prog, sos_2, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_2, rho - V, [x_vars;lx(1)], const_deg, sos_option, options);
    sos2_mult = sig{end};
    
    prog = prog.withEqs(subs(sig{end},[x_vars;lx(1)],[x0;0]));
    prog = prog.withEqs(subs(diff(sig{end},qd(1)),[x_vars;lx(1)],[x0;0]));
    prog = prog.withEqs(subs(diff(sig{end},qd(2)),[x_vars;lx(1)],[x0;0]));
    prog = prog.withEqs(subs(diff(sig{end},qd(3)),[x_vars;lx(1)],[x0;0]));
    prog = prog.withEqs(subs(diff(sig{end},s)*c + diff(sig{end},c)*-s,[x_vars;lx(1)],[x0;0]));
  else
    sos_2 = sos_2 - (rho - V)*sos2_mult;
  end

  prog = withSOS_fun(prog,sos_2);
end

if doSOS(3)
  [prog, sos_3, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_3, phi(1), [x_vars;lx(2)], const_deg, sos_option, options);
  [prog, sos_3, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_3, -phidot(2), [x_vars;lx(2)], const_deg, sos_option, options);
  [prog, sos_3, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_3, -lx(2)*psi(2), [x_vars;lx(2)], const_deg, sos_option, options);
  [prog, sos_3, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_3, mu^2*lzsq(2)^2 - lx(2)^2, [x_vars;lx(2)], const_deg, sos_option, options);
  [prog, sos_3, sig{end+1}, coefsig{end+1}] = spotless_add_eq_sprocedure(prog, sos_3, phi(2), [x_vars;lx(2)], const_deg);
  [prog, sos_3, sig{end+1}, coefsig{end+1}] = spotless_add_eq_sprocedure(prog, sos_3, (mu^2*lzsq(2)^2 - lx(2)^2)*psi(2), [x_vars;lx(2)], const_deg);  %should this be psi^2?
  
  if ~even(iter)
    [prog, sos_3, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_3, rho - V, [x_vars;lx(2)], const_deg, sos_option, options);
    sos3_mult = sig{end};
    
    prog = prog.withEqs(subs(sig{end},[x_vars;lx(2)],[x0;0]));
    prog = prog.withEqs(subs(diff(sig{end},qd(1)),[x_vars;lx(2)],[x0;0]));
    prog = prog.withEqs(subs(diff(sig{end},qd(2)),[x_vars;lx(2)],[x0;0]));
    prog = prog.withEqs(subs(diff(sig{end},qd(3)),[x_vars;lx(2)],[x0;0]));
    prog = prog.withEqs(subs(diff(sig{end},s)*c + diff(sig{end},c)*-s,[x_vars;lx(2)],[x0;0]));
    
  else
    sos_3 = sos_3 - (rho - V)*sos3_mult;
  end

  prog = withSOS_fun(prog,sos_3);
  
end

if doSOS(4)
  [prog, sos_4, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_4, phi(1), x_vars, const_deg, sos_option, options);
  [prog, sos_4, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_4, phi(2), x_vars, const_deg, sos_option, options);

  prog = withSOS_fun(prog,sos_4);
end



% Solve program

% options.do_fr = true;
%options.solveroptions.MSK_IPAR_INTPNT_BASIS = 'MSK_BI_NEVER';
options.regularize = false;
sol = prog.minimize(cost,sos_fun,options);

display(sprintf('Finished iter %d with cost=%f',iter,double(sol.eval(cost))))

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

Vsol = sol.eval(V);
R = double(sol.eval(rho));
sos1_mult = sol.eval(sos1_mult);
sos2_mult = sol.eval(sos2_mult);
sos3_mult = sol.eval(sos3_mult);
if ~even(iter),
  
  save(strcat(sprintf('gl_iter_%d',iter)),'Vsol','R','sos1_mult','sos2_mult','sos3_mult')
%   save iter_1 Vsol sos6_mult sos6_mult_2 sos6_mult_3 R Ao2 AI
else
  save(strcat(sprintf('gl_iter_%d',iter)),'Vsol','R','sos1_mult','sos2_mult','sos3_mult')
end
end
