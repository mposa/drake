megaclear

sos_option = 2;

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

g = 10;

prog = spotsosprog();

file_name = 'ballinhole_iter_%d';


%% Add indeterminate variables
q = msspoly('q',2);
qd = msspoly('qd',3);
s = msspoly('s',1);
c = msspoly('c',1);

nC = 5;
lx = msspoly('lx',nC);
lzsq = ones(nC,1);

x = q(1);
z = q(2);

xd = qd(1);
zd = qd(2);
thetad = qd(3);

v_vars = [q;qd];
x_vars = [q;s;c;qd];

prog = prog.withIndeterminate(q);
% prog = prog.withIndeterminate(s);
% prog = prog.withIndeterminate(c);
prog = prog.withIndeterminate(qd);
prog = prog.withIndeterminate(lx);

x0 = [zeros(3,1);1;zeros(3,1)];

options = spotprog.defaultOptions;
options.verbose = 1;
options.trig.enable = true;
options.trig.sin = s;
options.trig.cos = c;
options.clean_primal = false;
options.scale_monomials = true;
options.regularize = false;
options.regularize_eps = 1e-6;


%% Dynamics
ball_radius = 1;
hole_width = .5;
hole_depth = .5;
hole_center_u = .1;
mass = 1;
ball_inertia = .4*ball_radius^2*mass;
H = diag([1;1;ball_inertia]);
Hinv = inv(H);
K = 10;
f_free = [xd;zd;c*thetad;-s*thetad;-K*x;-g;0];

% for each contact, add set of inequalities
phi = [z - hole_depth;
       z - hole_depth; 
       x + hole_width/2; 
       -x + hole_width/2; 
       z];
contact_ineq{1} = [-x - hole_width/2];
contact_ineq{2} = [x - hole_width/2];
contact_ineq{3} = [hole_depth-z; phi(3) + hole_center_u; hole_center_u - phi(3)];
contact_ineq{4} = [hole_depth-z; phi(4) + hole_center_u; hole_center_u - phi(4)];
contact_ineq{5} = [];

skip_phi = [false;false;true;true;false];

J_N = [0 1 0;
       0 1 0;
       1 0 0;
      -1 0 0;
       0 1 0];
     
J_f = [1 0 ball_radius;
       1 0 ball_radius;
       0 -1 ball_radius;
       0 1 ball_radius;
       1 0 ball_radius];

f_free_guards{1} = [z - hole_depth];
f_free_guards{2} = [-z + hole_depth; -x + hole_width/2; x + hole_width/2;z];
% f_free_guards{1} = z;

%% Lyapunov function

[prog,V,coefv] = prog.newFreePoly(monomials(v_vars,0:V_degree));
[prog, equil_eqn] = prog.withEqs(subs(V,x_vars,x0));


Vdot_free = diff(V,x_vars)*f_free;
Vdot_impact = zeros(nC,1)*x;
for i=1:nC,
  Vdot_impact(i) = diff(V,qd)*Hinv*(J_N(i,:)'*lzsq(i) + J_f(i,:)'*lx(i));
end

prog_bkp = prog;
%% Add in constraints
prog = prog_bkp;

const_deg = 4;
sig = {};
coefsig = {};

mu = .2;

include_ball = false;
ball_fun = .2 - ([x;z;qd]'*diag([1;.5;ones(3,1)])*[x;z;qd]);

% SOS functions
% V(x) >= 0 for x admissable
% -Vdot_free(x) >= 0 for x admissable
% -Vdot_impact_1(x,l) >= 0 for (x,l) admissable

nSOS = 2*length(f_free_guards) + nC;
sos = zeros(nSOS,1)*x;
for i=1:length(f_free_guards),
  sos(i) = V - .1*(z + qd'*qd + x^2);
end

for i=1:length(f_free_guards),
  sos(i+length(f_free_guards)) = -Vdot_free - .0*(z);
end

for i=1:nC,
  sos(i+2*length(f_free_guards)) = -Vdot_impact(i);
end

for i=1:length(f_free_guards),
  for j=1:length(f_free_guards{i}),
    [prog, sos(i), sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos(i), f_free_guards{i}(j), [v_vars], const_deg, sos_option, options);
    [prog, sos(i+length(f_free_guards)), sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos(i+length(f_free_guards)), f_free_guards{i}(j), [v_vars], const_deg, sos_option, options);
    
    if include_ball
      [prog, sos(i), sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos(i), ball_fun, [v_vars], const_deg, sos_option, options);
      [prog, sos(i+length(f_free_guards)), sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos(i+length(f_free_guards)), ball_fun, [v_vars], const_deg, sos_option, options);
    end
  end
end

for i=1:nC,
  sos_ind = 2*length(f_free_guards) + i;
  sos_vars = [v_vars;lx(i)];
  for j=1:length(contact_ineq{i}),
    [prog, sos(sos_ind), sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos(sos_ind), contact_ineq{i}(j), sos_vars, const_deg, sos_option, options);
  end
  
  [prog, sos(sos_ind), sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos(sos_ind), -J_N(i,:)*qd, sos_vars, const_deg, sos_option, options);
  [prog, sos(sos_ind), sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos(sos_ind), -J_f(i,:)*qd*lx(i), sos_vars, const_deg, sos_option, options);
  [prog, sos(sos_ind), sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos(sos_ind), mu^2*lzsq(i)^2 - lx(i)^2, sos_vars, const_deg, sos_option, options);
  
  if ~skip_phi(i)
    [prog, sos(sos_ind), sig{end+1}, coefsig{end+1}] = spotless_add_eq_sprocedure(prog, sos(sos_ind), phi(i), sos_vars, const_deg);
  end
  [prog, sos(sos_ind), sig{end+1}, coefsig{end+1}] = spotless_add_eq_sprocedure(prog, sos(sos_ind), (mu^2*lzsq(i)^2 - lx(i)^2)*J_f(i,:)*qd, sos_vars, const_deg-1);
  
  if include_ball
    [prog, sos(sos_ind), sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos(sos_ind), ball_fun, sos_vars, const_deg, sos_option, options);
  end
end

for i=1:nSOS,
  prog = withSOS_fun(prog,sos(i));
end

% Solve program

prog = prog.withPos(10 - subs(V,[z;x;qd],[1;zeros(4,1)]));
% cost = ;
cost = 0;
tic
sol = prog.minimize(cost,sos_fun,options);
toc
Vsol = sol.eval(V);
