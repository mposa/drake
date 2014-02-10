% clear all
megaclear
iter = 2;
sos_option = 1;
do_backoff = 0;
do_clean = 0;
do_eig = 0;

switch sos_option
  case 1
%     sos_fun = @spot_sedumi;
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
V_degree = 4;

g = 9.81;

prog = spotsosprog();

%% Add indeterminate variables
q = msspoly('q',3);
qd = msspoly('qd',3);
% qdp = msspoly('qdp',3);
qd_hess = msspoly('qdh',3);
s_vec = msspoly('s',3);
c_vec = msspoly('c',3);
lx = msspoly('lx',2);
lzsq = msspoly('lzs',2);

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
prog = prog.withIndeterminate(qd_hess);
prog = prog.withIndeterminate(lx);
prog = prog.withIndeterminate(lzsq);

x0 = [0;0;1;0;0;0];




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

qdp1 = qd + f_impact(5:7,1);
qdp2 = qd + f_impact(5:7,2);

phi = [z + s*spi8 - c*cpi8 + cpi8; z - s*spi8 - c*cpi8 + cpi8];

phidot = [zd + c*pitchd*spi8 + s*pitchd*cpi8; zd - c*pitchd*spi8 + s*pitchd*cpi8];
phidotp = [subs(phidot(1),qd,qdp1);subs(phidot(2),qd,qdp2)];

psi = [xd + c*pitchd*cpi8 - s*pitchd*spi8; xd + c*pitchd*cpi8 + s*pitchd*spi8];
psip = [subs(psi(1),qd,qdp1);subs(psi(2),qd,qdp2)];


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
  load(sprintf('iter_%d',iter-1))
end
if iter==0,
  [prog,V,coefv] = prog.newFreePoly(monomials(v_vars,0:V_degree));
  [prog, equil_eqn] = prog.withEqs(subs(V,[z;s;c;qd],[0;0;1;0;0;0]));
  
  prog = prog.withEqs(subs(diff(V,qd(1)),x_vars,x0));
  prog = prog.withEqs(subs(diff(V,qd(2)),x_vars,x0));
  prog = prog.withEqs(subs(diff(V,qd(3)),x_vars,x0));
  prog = prog.withEqs(subs(diff(V,s)*c + diff(V,c)*-s,x_vars,x0));
  
  prog = withSOS_fun(prog, qd_hess'*diff(diff(V,qd)',qd)*qd_hess);
elseif iter==1,
%   load iter_0
  V = Vsol;
elseif ~even(iter),
%   load iter_2
  V = Vsol;
else
%   load iter_1
  [prog,V,coefv] = prog.newFreePoly(monomials(v_vars,0:V_degree));
  [prog, equil_eqn] = prog.withEqs(subs(V,[z;s;c;qd],[0;0;1;0;0;0]));
  
  prog = prog.withEqs(subs(diff(V,qd(1)),x_vars,x0));
  prog = prog.withEqs(subs(diff(V,qd(2)),x_vars,x0));
  prog = prog.withEqs(subs(diff(V,qd(3)),x_vars,x0));
  prog = prog.withEqs(subs(diff(V,s)*c + diff(V,c)*-s,x_vars,x0));
  
  prog = withSOS_fun(prog, qd_hess'*diff(diff(V,qd)',qd)*qd_hess);
end

Vdot_free = diff(V,[x;z;s;c;xd;zd;pitchd])*f_free;
Vdot_impact_1 = subs(V,qd,qdp1) - V;
Vdot_impact_2 = subs(V,qd,qdp2) - V;


%% Ball constraints
ball_vec = [z;s;1-c;qd];
Ao2 = zeros(6);
Ao2(1,1) = 100;
Ao2(2,2) = 2;
Ao2(3,3) = 2;
Ao2(4,4) = .5;
Ao2(5,5) = .5;
Ao2(6,6) = .1;

if iter==0,
  cost = 0;
  rho = .2;
  Ai = Ao2;
%   Ai = diag([50;6;10;10;10;10]);
%   Ai = diag([50;50;50;50;50;50]);
%   rho = .1;
else
%   [prog,rho] = prog.newFree(1);
%   cost = -rho;

%   Ai = Ao2;
%   [prog, Ai] = PSD_fun(prog,6);

  [prog,Ai_diag] = prog.newPos(6);
  Ai = diag(Ai_diag);
  cost = Ai(1,1) + .1*trace(Ai);

%   [prog,Ascale] = prog.newFree(1);
%   cost = Ascale;
%   Ai = Ascale*diag([10;1;1;2;2;2]);  
  
  rho = .1;
end

rho_i = rho;
rho_o = 1;




h_Bo2 = rho_o - ball_vec'*Ao2*ball_vec;

h_Bi = rho_i - ball_vec'*Ai*ball_vec;


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

if iter==0,
  sos_6 = (1 - V);
else
%   sos_6 = -h_Bi*tmp_mult;
  sos_6 = -h_Bi*(1 + z^2 + qd'*qd + 1 - c);
end

prog_bkp = prog;
%% Add in constraints
prog = prog_bkp;


const_deg = 4;
sig = {};
coefsig = {};


doSOS = [1 1 1 1 1 1 0];
% doSOS = [1 1 1 1 1 0 1];

if ~even(iter) && 0
  doSOS = [0 0 0 0 0 1 0];
end

if doSOS(1)
  % non-penetration (admissability of x)
  [prog, sos_1, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_1, phi(1), [x_vars], const_deg, sos_option);
  [prog, sos_1, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_1, phi(2), [x_vars], const_deg, sos_option);
  [prog, sos_1, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_1, h_Bo2, [x_vars], const_deg, sos_option);
  
  prog = prog.withEqs(subs(sig{end},x_vars,x0));
  prog = prog.withEqs(subs(diff(sig{end},qd(1)),x_vars,x0));
  prog = prog.withEqs(subs(diff(sig{end},qd(2)),x_vars,x0));
  prog = prog.withEqs(subs(diff(sig{end},qd(3)),x_vars,x0));
  prog = prog.withEqs(subs(diff(sig{end},s)*c + diff(sig{end},c)*-s,x_vars,x0));
  
  prog = withSOS_fun(prog,sos_1);
end

if doSOS(2)
  [prog, sos_2, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_2, phi(2), [x_vars;lzsq(1);lx(1)], const_deg, sos_option);
  % Contact constraints (admissability of lambda)
  [prog, sos_2, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_2, -phidot(1), [x_vars;lzsq(1);lx(1)], const_deg, sos_option);
  [prog, sos_2, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_2, -lx(1)*psip(1), [x_vars;lzsq(1);lx(1)], const_deg, sos_option);
  [prog, sos_2, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_2, lzsq(1)^4 - lx(1)^2, [x_vars;lzsq(1);lx(1)], const_deg, sos_option);
  [prog, sos_2, sig{end+1}, coefsig{end+1}] = spotless_add_eq_sprocedure(prog, sos_2, phi(1), [x_vars;lzsq(1);lx(1)], const_deg);
  [prog, sos_2, sig{end+1}, coefsig{end+1}] = spotless_add_eq_sprocedure(prog, sos_2, phidotp(1), [x_vars;lzsq(1);lx(1)], const_deg);
  [prog, sos_2, sig{end+1}, coefsig{end+1}] = spotless_add_eq_sprocedure(prog, sos_2, (lzsq(1)^4 - lx(1)^2)*psip(1), [x_vars;lzsq(1);lx(1)], const_deg);  %should this be psi^2?
  [prog, sos_2, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_2, h_Bo2, [x_vars;lzsq(1);lx(1)], const_deg, sos_option);
  
  prog = prog.withEqs(subs(sig{end},[x_vars;lzsq(1);lx(1)],[x0;0;0]));
  prog = prog.withEqs(subs(diff(sig{end},qd(1)),[x_vars;lzsq(1);lx(1)],[x0;0;0]));
  prog = prog.withEqs(subs(diff(sig{end},qd(2)),[x_vars;lzsq(1);lx(1)],[x0;0;0]));
  prog = prog.withEqs(subs(diff(sig{end},qd(3)),[x_vars;lzsq(1);lx(1)],[x0;0;0]));
  prog = prog.withEqs(subs(diff(sig{end},s)*c + diff(sig{end},c)*-s,[x_vars;lzsq(1);lx(1)],[x0;0;0]));
  
  prog = withSOS_fun(prog,sos_2);
end

if doSOS(3)
  [prog, sos_3, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_3, phi(1), [x_vars;lx(2)], const_deg, sos_option);
  [prog, sos_3, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_3, -phidot(2), [x_vars;lzsq(2);lx(2)], const_deg, sos_option);
  [prog, sos_3, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_3, -lx(2)*psip(2), [x_vars;lzsq(2);lx(2)], const_deg, sos_option);
  [prog, sos_3, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_3, lzsq(2)^4 - lx(2)^2, [x_vars;lzsq(2);lx(2)], const_deg, sos_option);
  [prog, sos_3, sig{end+1}, coefsig{end+1}] = spotless_add_eq_sprocedure(prog, sos_3, phi(2), [x_vars;lzsq(2);lx(2)], const_deg);
  [prog, sos_3, sig{end+1}, coefsig{end+1}] = spotless_add_eq_sprocedure(prog, sos_3, phidotp(2), [x_vars;lzsq(2);lx(2)], const_deg);
  [prog, sos_3, sig{end+1}, coefsig{end+1}] = spotless_add_eq_sprocedure(prog, sos_3, (lzsq(2)^4 - lx(2)^2)*psip(2), [x_vars;lzsq(2);lx(2)], const_deg);  %should this be psi^2?
  [prog, sos_3, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_3, h_Bo2, [x_vars;lzsq(2);lx(2)], const_deg, sos_option);
  
  prog = prog.withEqs(subs(sig{end},[x_vars;lzsq(2);lx(2)],[x0;0;0]));
  prog = prog.withEqs(subs(diff(sig{end},qd(1)),[x_vars;lzsq(2);lx(2)],[x0;0;0]));
  prog = prog.withEqs(subs(diff(sig{end},qd(2)),[x_vars;lzsq(2);lx(2)],[x0;0;0]));
  prog = prog.withEqs(subs(diff(sig{end},qd(3)),[x_vars;lzsq(2);lx(2)],[x0;0;0]));
  prog = prog.withEqs(subs(diff(sig{end},s)*c + diff(sig{end},c)*-s,[x_vars;lzsq(2);lx(2)],[x0;0;0]));
  
  prog = withSOS_fun(prog,sos_3);
  
end

if doSOS(4)
  [prog, sos_4, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_4, phi(1), x_vars, const_deg, sos_option);
  [prog, sos_4, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_4, phi(2), x_vars, const_deg, sos_option);
  [prog, sos_4, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_4, h_Bo2, x_vars, const_deg, sos_option);
  
  prog = prog.withEqs(subs(sig{end},x_vars,x0));
  prog = prog.withEqs(subs(diff(sig{end},qd(1)),x_vars,x0));
  prog = prog.withEqs(subs(diff(sig{end},qd(2)),x_vars,x0));
  prog = prog.withEqs(subs(diff(sig{end},qd(3)),x_vars,x0));
  prog = prog.withEqs(subs(diff(sig{end},s)*c + diff(V,c)*-s,x_vars,x0));
  
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
  %     [prog, sos_6, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_6, phi(1)*phi(2), x_vars, const_deg, sos_option);
  
  if iter==0,
    [prog, sos_6, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_6, h_Bi, x_vars, const_deg, sos_option);
  elseif ~even(iter),
    %     [prog, sos_6, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_6, (V - 1)*phi(1), x_vars, const_deg, sos_option);
    %     [prog, sos_6, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_6, (V - 1)*phi(2), x_vars, const_deg, sos_option);
    [prog, sos_6, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_6, V - 1, x_vars, const_deg, sos_option);
  else
    %     sos_6 = sos_6 - (V-1)*phi(1)*sos6_mult_3;
    %     sos_6 = sos_6 - (V-1)*phi(2)*sos6_mult_2;
    sos_6 = sos_6 - (V-1)*sos6_mult;
  end
  
  prog = withSOS_fun(prog,sos_6);
end

if doSOS(7)
  if iter==0
    sos_7 = subs(1-V,[s;c;qd],[0;1;zeros(3,1)]);
    [prog, sos_7, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_7, phi(1), x_vars, const_deg, sos_option);
    [prog, sos_7, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_7, phi(2), x_vars, const_deg, sos_option);
    [prog, sos_7, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_7, .15^2-z^2, z, const_deg, sos_option);
    prog = withSOS_fun(prog,sos_7);
  end
end


% Solve program


options = spotprog.defaultOptions;
options.verbose = 1;
options.trig.enable = true;
options.trig.sin = s;
options.trig.cos = c;
sol = prog.minimize(cost,sos_fun,options);

sqrt(1/double(sol.eval(Ai(1,1)*10)))

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
    prog = prog.withEqs(prog.freeVar(abs(double(sol.eval(prog.freeVar))) < 1e-6));
    sol = prog.minimize(cost,sos_fun,options);
  end
  
  if do_backoff
    costval = double(sol.eval(cost));

    prog = prog.withPos(costval + .1*abs(costval) - cost);
    cost = 0;
    
%     [prog,gamma] = prog.newFree(1);
%     cost = gamma;
%     for i=1:length(prog.sosExpr),
%       prog.sosExpr(i) = prog.sosExpr(i) + gamma*sol.gramMonomials{i}'*sol.gramMonomials{i};
%     end
    sol = prog.minimize(cost,sos_fun,options);
  end
end

Vsol = sol.eval(V);
R = double(sol.eval(rho));
AI = double(sol.eval(Ai));
if iter==0,
  save iter_0 Vsol R Ao2 AI
elseif ~even(iter),
  sos6_mult = sol.eval(sig{end});
%   sos6_mult_2 = sol.eval(sig{end-1});
%   sos6_mult_3 = sol.eval(sig{end-2});
  save(strcat(sprintf('iter_%d',iter)),'Vsol','sos6_mult','R','Ao2','AI')
%   save iter_1 Vsol sos6_mult sos6_mult_2 sos6_mult_3 R Ao2 AI
else
%   save iter_2 Vsol R Ao2 AI
save(strcat(sprintf('iter_%d',iter)),'Vsol','R','Ao2','AI')
end
