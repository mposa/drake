function [p_sol] = normUpperBound(f,R_diag, options)
% Solve for a p which approximates 


if ~isfield(options,'time_varying')
  options.time_varying = false;
end

if ~isfield(options,'T')
  options.T = 1;
end

if ~isfield(options,'W')
  options.W = [];
end
W = options.W;

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


%% SOS constraints
% State constraint
A = diag(1./(R_diag.^2));
h_X = 1 - x'*A*x;

mult_degree = even_degree(p+f,vars);

% p >= f
[prog, p_sos] = spotless_add_sprocedure(prog, (p-f), h_X,vars,mult_degree-2);
if ~isempty(W)
  W_degree = even_degree(W,vars);
  [prog, p_sos] = spotless_add_sprocedure(prog, p_sos, W,vars,mult_degree-W_degree);
end
if time_varying
  [prog, p_sos] = spotless_add_sprocedure(prog, p_sos, t*(T-t),vars,mult_degree-2);
end
prog = prog.withSOS(p_sos);

% p >= -f
[prog, p_sos] = spotless_add_sprocedure(prog, (p+f), h_X,vars,mult_degree-2);
if time_varying
  [prog, p_sos] = spotless_add_sprocedure(prog, p_sos, t*(T-t),vars,mult_degree-2);
end
if ~isempty(W)
  [prog, p_sos] = spotless_add_sprocedure(prog, p_sos, W,vars,mult_degree-W_degree);
end
prog = prog.withSOS(p_sos);


if ~isempty(W)
  [prog, p_sos] = spotless_add_sprocedure(prog, p-1, h_X,vars,degree-2);
  [prog, p_sos] = spotless_add_sprocedure(prog, p_sos, -W,vars,degree-W_degree);
%   prog = prog.withSOS(p_sos);
end

%% Solve
if time_varying
  cost = spotlessIntegral(prog,p,x,R_diag,t,[0 T]);
else
  cost = spotlessIntegral(prog,p,x,R_diag,[],[]);
end


spot_options = spotprog.defaultOptions;
spot_options.verbose = true;
spot_options.do_fr = true;
solver = @spot_mosek;
% solver = @spot_sedumi;
sol = prog.minimize(cost,solver,spot_options);

p_sol = sol.eval(p);
end