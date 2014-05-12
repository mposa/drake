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

%%
N_obs = 2;
obs_pos = randn(2,N_obs);
obs_rad = ones(2,N_obs)+45;
u_degree = 3;
const_deg = 4;

x = msspoly('x',1);
y = msspoly('y',1);
s = msspoly('s',1);
c = msspoly('c',1);
% u = msspoly('u',1);
prog = spotsosprog();
prog = prog.withIndeterminate(x);
prog = prog.withIndeterminate(y);
prog = prog.withIndeterminate(s);
prog = prog.withIndeterminate(c);


[prog,u,coefu] = prog.newFreePoly(monomials([x;y;s;c],0:u_degree));

v = 1;
% [x;y;s;c]
f = [s*v;c*v;u*c;-u*s];
  

options = spotprog.defaultOptions;
options.verbose = 1;
options.trig.enable = true;
options.trig.sin = s;
options.trig.cos = c;
options.clean_primal = false;
options.scale_monomials = true;
options.regularize = false;
options.regularize_eps = 1e-6;

sig = {};
coefsig = {};

sos = zeros(N_obs,1)*x;
for i=1:N_obs,
  k = 10;
  V_i = (x - obs_pos(1,i))^2 + (y - obs_pos(2,i))^2 + (x - obs_pos(1,i) + k*s)^2 + (y - obs_pos(2,i) + k*c)^2;
  x_vec = [x - obs_pos(1,i);y - obs_pos(2,i)];
  sos(i) = diff(V_i,[x;y;s;c])*f - .01;
  [prog, sos(i), sig{end+1}, coefsig{end+1}] = spotless_add_eq_sprocedure(prog, sos(i), V_i - obs_rad(i) , [x;y;s;c], const_deg);
  prog = withSOS_fun(prog,sos(i));
end

sol = prog.minimize(0,sos_fun,options);