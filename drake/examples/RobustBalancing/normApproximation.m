function [p_sol,q_sol] = normApproximation(f,R_diag, options)
% Solve for a p which approximates 


if ~isfield(options,'time_varying')
  options.time_varying = false;
end

if ~isfield(options,'T')
  options.T = 1;
end

nX = length(R_diag);

%% Solution method settings
degree = options.degree; % degree of V,W

time_varying = options.time_varying;
T = options.T;

%% Create SOS program
prog = spotsosprog;

%% Create indeterminates
[prog,x] = prog.newIndeterminate('x',nX);
vars = x;
if time_varying
  [prog,t] = prog.newIndeterminate('t',1);
  vars = [vars;t];
end

%% Create polynomials p(t,x) and q(t,x)
pbasis = monomials(vars,0:degree);
[prog,p,pcoeff] = prog.newFreePoly(pbasis);

[prog, Q] = prog.newPSD(length(pbasis));
q = pbasis'*Q*pbasis;

% p^2 <= q
prog = prog.withPSD([Q pcoeff; pcoeff' 1]);

%% SOS constraints
% State constraint
A = diag(1./(R_diag.^2));
h_X = 1 - x'*A*x;

% p >= 0
[prog, p_sos] = spotless_add_sprocedure(prog, p, h_X,vars,degree-2);
if time_varying
  [prog, p_sos] = spotless_add_sprocedure(prog, p_sos, t*(T-t),vars,degree-2);
end
prog = prog.withSOS(p_sos);

% q*f <= 1
qf_degree = 2*degree+even_degree(f,vars);
[prog, q_sos] = spotless_add_sprocedure(prog, 1-q*f, h_X,vars,qf_degree-2);
[prog, q_sos] = spotless_add_sprocedure(prog, q_sos, h_X*f,vars,2*degree-2);
if time_varying
  [prog, q_sos] = spotless_add_sprocedure(prog, q_sos, t*(T-t),vars,qf_degree-2);
end
prog = prog.withSOS(q_sos);

%% Solve
if time_varying
  cost = spotlessIntegral(prog,-p,x,R_diag,t,[0 T]);
else
  cost = spotlessIntegral(prog,-p,x,R_diag,[],[]);
end


spot_options = spotprog.defaultOptions;
spot_options.verbose = true;
spot_options.do_fr = true;
solver = @spot_mosek;
% solver = @spot_sedumi;
sol = prog.minimize(cost,solver,spot_options);

p_sol = sol.eval(p);
q_sol = sol.eval(q);
end