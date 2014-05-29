megaclear 
sos_option = 1;

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

prog = spotsosprog();

F = 1;

q = msspoly('q',1);
v = msspoly('v',1);
f = msspoly('f',1);
qh = msspoly('qh',1);
vh = msspoly('vh',1);

fe = msspoly('fe',1);
fh_scale = msspoly('fs',1);
fh = msspoly('fh',1);
% fh_sgn = msspoly('fg',1);

prog = prog.withIndeterminate(q);
prog = prog.withIndeterminate(v);
prog = prog.withIndeterminate(f);
prog = prog.withIndeterminate(qh);
prog = prog.withIndeterminate(vh);
% prog = prog.withIndeterminate(fh_scale);
% prog = prog.withIndeterminate(fh_sgn);
prog = prog.withIndeterminate(fe);
prog = prog.withIndeterminate(fh);

fh_sgn = -v;
%%
K = [10;10];
Kf = 10;

y = q;
yh = qh;

% fh = fh_sgn*fh_scale^2;

x_dyn = [v;f];
xh_dyn = [vh;fh] + K*(y - yh);
fe_dyn = Kf*(y-yh);

[prog,Q] = PSD_fun(prog,3);
% [prog, gamma] = prog.newFree(1);
prog = prog.withEqs(trace(Q) - 1);

x_err = [q-qh;v-vh;F-fe];

V = x_err'*Q*x_err;
Vdot = diff(V,[q;v])*x_dyn + diff(V,[qh;vh])*xh_dyn + diff(V,fe)*fe_dyn;

ball = .1 - x_err'*x_err;

%%
x_vars = [q;v;f;qh;vh;fh;fe];
sos = -Vdot;

const_deg = 4;
sig = {};
coefsig = {};
[prog, sos, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos, -f*v, x_vars, const_deg, sos_option);
[prog, sos, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos, F-f, x_vars, const_deg, sos_option);
[prog, sos, sig{end+1}, coefsig{end+1}] = spotless_add_eq_sprocedure(prog, sos, v*(F-f), x_vars, const_deg);

[prog, sos, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos, -fh*vh, x_vars, const_deg, sos_option);
[prog, sos, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos, fe-fh, x_vars, const_deg-2, sos_option);
[prog, sos, sig{end+1}, coefsig{end+1}] = spotless_add_eq_sprocedure(prog, sos, v*(fe-fh), x_vars, const_deg-2);

% [prog, sos, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos, ball, x_vars, const_deg, sos_option);
[prog, sos, sig{end+1}, coefsig{end+1}] = spotless_add_eq_sprocedure(prog, sos, ball, x_vars, const_deg);

prog = withSOS_fun(prog,sos);

cost = 0;
options = spotprog.defaultOptions;
options.verbose = 1;
sol = prog.minimize(cost,sos_fun,options);