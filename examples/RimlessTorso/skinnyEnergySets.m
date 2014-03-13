megaclear

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


prog = spotsosprog();

%% Add indeterminate variables
q = msspoly('q',4);
qd = msspoly('qd',4);
s_vec = msspoly('s',4);
c_vec = msspoly('c',4);

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


x0 = double(subs(x_vars,[q;qd;s_vec;c_vec],[zeros(12,1);ones(4,1)]));

options = spotprog.defaultOptions;
options.verbose = 1;
options.verbose = 1;
options.trig.enable = true;
options.trig.sin = [s;s_th];
options.trig.cos = [c;c_th];
options.clean_primal = true;
options.scale_monomials = true;
options.regularize = true;

change_variables = true;

%% Dynamics
[H,C,B,phi,phidot,psi,J,J_f,K,S,U] = skinnyEOM_mss(q,qd,s_vec,c_vec);

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
U = clean(U);

E = .5*qd'*H*qd + U;
V = E;

prog_bkp = prog;

M = [[1 0 0 0 0; 0 c -s 0 0; 0 0 0 c_th -s_th] zeros(3,4); zeros(4,5) eye(4)];
gradE = M*diff(E,[z;s;c;s_th;c_th;qd])';
hessE = M*diff(gradE,[z;s;c;s_th;c_th;qd])';

%% Add in constraints
prog = prog_bkp;


const_deg = 4;
sig = {};
coefsig = {};


% Ball constraints
ball_vec = [z;s;1-c;s_th;1-c_th;qd];

% V = subs(V,[z;qd;s;c],[zeros(5,1);0;1]);
% ball_vec = [0;0;0;s_th;1-c_th;zeros(4,1);];


Ao2 = zeros(9);
Ao2(1,1) = 400;
Ao2(2,2) = 5;
Ao2(3,3) = 5;
Ao2(4,4) = .1;
Ao2(5,5) = .1;
Ao2(6:9,6:9) = .25*double(subs(H,x_vars,x0));

% Ao2 = subs(hessE,x_vars,x0);

[prog, rho_V] = prog.newFree(1);

cost = -rho_V;
rho_i = .03;

Ai = Ao2;
Ao = Ao2;
rho_o = .1;


h_Bo = rho_o - ball_vec'*Ao*ball_vec;
h_Bi = rho_i - ball_vec'*Ai*ball_vec;


sos_1 = V - rho_V;
sos_2 = rho_V - V;

doSOS = [1 0];

if doSOS(1),
  [prog, sos_1, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_1, phi(1), x_vars, const_deg, sos_option, options);
  [prog, sos_1, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_1, phi(2), x_vars, const_deg, sos_option, options);
  [prog, sos_1, sig{end+1}, coefsig{end+1}] = spotless_add_eq_sprocedure(prog, sos_1, h_Bo, x_vars, const_deg);
  prog = withSOS_fun(prog,sos_1);
end


if doSOS(2)
  [prog, sos_2, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_2, phi(1), x_vars, const_deg, sos_option, options);
  [prog, sos_2, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_2, phi(2), x_vars, const_deg, sos_option, options);
  [prog, sos_2, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_2, h_Bi, x_vars, const_deg, sos_option, options);
  prog = withSOS_fun(prog,sos_2);
end

sol = prog.minimize(cost,sos_fun,options);