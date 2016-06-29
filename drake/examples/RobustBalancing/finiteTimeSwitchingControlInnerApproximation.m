function [Vsol,Wsol] = finiteTimeSwitchingControlInnerApproximation(model,T,B,R_diag,target,options)

%scaling of state vector
if ~isfield(options,'scale')
  scale = ones(model.num_states,1);
elseif length(options.scale) == 1
  scale = ones(model.num_states,1)*options.scale;
else
  scale = options.scale;
end

if ~isfield(options,'beta')
  beta = 0;
else
  beta = options.beta;
end

scale_inv = 1./scale;
degree = options.degree; % degree of V,W

%% Create SOS program
prog = spotsosprog;


%% Create indeterminates

[prog,x]=prog.newIndeterminate('x', model.num_states); % state
[prog,t]=prog.newIndeterminate('t', 1); % time
u_mss = msspoly('u',model.num_inputs);

%% Scale r_diag


%% Create polynomials V(t,x) and W(x)
V_vars = [t;x];
[prog,V] = prog.newFreePoly(monomials(V_vars,0:degree));

W_vars = x;
[prog,W] = prog.newFreePoly(monomials(W_vars,0:degree));

%% Dynamics
[f,g] = model.controlAffineDynamics([],scale_inv.*x);

nU = size(g,2);
ndgrid_arg = mat2cell(repmat([-1;1],1,nU),2,ones(1,nU)');
[ugrid{1:nU}] = ndgrid(ndgrid_arg{:});
umat = zeros(2^nU,nU);
for i=1:nU,
  umat(:,i) = ugrid{i}(:);
end

f = scale.*f;
g = repmat(scale,1,nU).*g.*g;

T_init = T;
f = f*T;
g = g*T;
T = 1;

%% Goal region

V0p = subs(target,x,scale_inv.*x);

% State constraint
R_diag = scale'.*R_diag;
A = diag(1./(R_diag.^2));
h_X = 1 - x'*A*x;


%% SOS constraints

% (1) w(x) >= 1 + V(0,x) for x in X
[prog, wv_sos] = spotless_add_sprocedure(prog, W-subs(V,t,0) - 1, h_X,W_vars,degree-2);
prog = prog.withSOS(wv_sos);

% (2) w(x) >=0 for x in X
[prog, w_sos] = spotless_add_sprocedure(prog, W, h_X,W_vars,degree-2);
prog = prog.withSOS(w_sos);

% (3) V(x) >=0 for x on boundary of state space
[prog, v_sos] = spotless_add_eq_sprocedure(prog, V, h_X,V_vars,degree-2);
[prog, v_sos] = spotless_add_sprocedure(prog, v_sos, t*(T-t),V_vars,degree-2);
prog = prog.withSOS(v_sos);




% (4) Vdot(x) <= 0 for x in X
for j=1:2^nU
  Vdot = diff(V,x)*(f + g*umat(j,:)') + diff(V,t);
  Vdot_degree = even_degree(Vdot,x);
  [prog, Vdot_sos] = spotless_add_sprocedure(prog,beta*V-Vdot, h_X,V_vars,Vdot_degree-2);
  [prog, Vdot_sos] = spotless_add_sprocedure(prog, Vdot_sos, t*(T-t),V_vars,Vdot_degree-2);
  for k=1:nU,
    B_degree = even_degree(B(k),x);
    [prog, Vdot_sos] = spotless_add_sprocedure(prog, Vdot_sos, umat(j,k)*B(k),x,Vdot_degree - B_degree);
  end
  prog = prog.withSOS(Vdot_sos);
end



% (5) V(T,x) >= 0 for x not in goal region
[prog, Vgoal_sos] = spotless_add_sprocedure(prog,subs(V,t,T), h_X,W_vars,degree-2);
[prog, Vgoal_sos] = spotless_add_sprocedure(prog,Vgoal_sos, -V0p,W_vars,degree-2);
prog = prog.withSOS(Vgoal_sos);

%% Set up cost function -- integration over a sphere

cost = spotlessIntegral(prog,W,x,R_diag,[],[]);

spot_options = spotprog.defaultOptions;
spot_options.verbose = true;
spot_options.do_fr = true;
solver = @spot_mosek;
% solver = @spot_sedumi;
sol = prog.minimize(cost,solver,spot_options);



%% Plotting
Vsol = subs(sol.eval(V),x,scale.*x);
Wsol = subs(sol.eval(W),x,scale.*x);
model.plotfun(0, sum(Vsol), Wsol, subs(h_X,x,scale.*x), R_diag, msspoly('t',1), x, []);
sum(Vsol)
end