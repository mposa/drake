function u_bnd = boundController(model,uCoeff,uMomentBasis,T,R_diag,options)
% Implementation of Korda et. al approach. Finds the controller u_bnd(t,x)
% that approximates u(t,x) with the added restriciton that 1 <= u_bnd <= 1
% for all t in [0,T] and x in X, where X is defined using R_diag
% As in similar scripts, R_diag generates an axis-aligned ellipsoid
% centered at the origin with elliptic radii corresponding to R_diag
%
% see nStepCapturabilitySOS for how uCoeff and uMomentBasis are generated
%   the correspond to the coefficients of u(t,x) and the monomial basis for
%   the moment matrix, aso where 
%   u_i(t,x) = sum_k uCoeff{i}(k)*uMomentBasis{i}(k,1)

if nargin < 6
  options = struct();
end

if ~isfield(options,'zero_origin')
  options.zero_origin = false;
end

%% Compute moment matrix for lebesgue measure
  % check for an error
  for i=2:length(uMomentBasis)
    double(uMomentBasis{i}-uMomentBasis{i-1});
  end
  prog = spotsosprog;
  [prog,z] = prog.newFree(1);
  [prog,x] = prog.newIndeterminate('x',model.num_states);
  [prog,t] = prog.newIndeterminate('t',1);
  
  [vars,alphas,coeff] = decomp(uMomentBasis{1}(:));
  
  %confirm all coeffs equal 1
  assert(~any(sum(abs(coeff),2) ~= 1)); 
  assert(~any(coeff(:) ~= 0 & coeff(:) ~= 1))
  
  x_inds = [];
  x_var_inds = [];
  x_var_inds_c = [];
  % extract dependence on t and x
  for i=1:model.num_states,
    is_var = false;
    for j=1:length(vars),
      if isequal(x(i),vars(j)),
        x_inds(end+1) = j;
        x_var_inds(end+1) = i;
        is_var = true;
        break;
      end
    end
    if ~is_var
      x_var_inds_c(end+1) = i;
    end
  end
  
  t_ind = [];
  for j=1:length(vars),
    if isequal(t,vars(j))
      t_ind = j;
    end
  end
  
  alphas_all = zeros(size(alphas,1),model.num_states+1);
  alphas_all(:,x_var_inds) = alphas(:,x_inds);
  if ~isempty(t_ind)
    alphas_all(:,end) = alphas(:,t_ind);
  end
  
  % ensure all vars are represented
  assert(length([x_inds t_ind]) == length(vars))
  
  M = zeros(size(uMomentBasis{1}));
  for i=1:numel(M),
    I = find(coeff(i,:));    
    M(i) = monomialIntegral(alphas_all(I,:),1:model.num_states,R_diag,model.num_states+1,[0 T]);
  end
  
  %%
  nU = length(uCoeff);
  u_bnd = msspoly();
  
  % State constraint
  A = diag(1./(R_diag.^2));
  h_X = 1 - x'*A*x;
  
  for i=1:nU
    %%
    prog = spotsosprog;
    [prog,x] = prog.newIndeterminate('x',model.num_states);
    

    [prog,t] = prog.newIndeterminate('t',1);
    % check for t dependence
    time_dependent = ~isequal(msspoly(0),sum(diff(uMomentBasis{1}(:,1),t)));
    if time_dependent;
      vars = x;
    else
      vars = [t;x];
    end
    
    [prog,u_bnd_coeff{i}] = prog.newFree(length(uCoeff{i}));
    u_bnd_i = u_bnd_coeff{i}'*uMomentBasis{i}(:,1);
    
    if options.zero_origin
      prog = prog.withEqs(subs(u_bnd_i,vars,vars*0) - .5);
    end
    
    u_deg = even_degree(u_bnd_i,x);
    
    % u>= 0
    [prog, u_pos] = spotless_add_sprocedure(prog, u_bnd_i, h_X,vars,u_deg-2);
    if time_dependent
      [prog, u_pos] = spotless_add_sprocedure(prog, u_pos, t*(T-t),vars,u_deg-2);
    end
    prog = prog.withSOS(u_pos);
    % u<= 1
    [prog, u_max] = spotless_add_sprocedure(prog, 1-u_bnd_i, h_X,vars,u_deg-2);
    if time_dependent
      [prog, u_max] = spotless_add_sprocedure(prog, u_max, t*(T-t),vars,u_deg-2);
    end
    prog = prog.withSOS(u_max);
    
    % objective slack via Schur complement
    L = chol(M);
    [prog,gamma] = prog.newFree(1);
    
    prog = prog.withPSD([eye(size(M)) L*u_bnd_coeff{i};u_bnd_coeff{i}'*L' gamma]);
    cost = gamma - 2*uCoeff{i}'*M*u_bnd_coeff{i};
    % Solve
    spot_options = spotprog.defaultOptions;
    spot_options.verbose = true;
    spot_options.do_fr = true;
    solver = @spot_mosek;
    % solver = @spot_sedumi;
    sol = prog.minimize(cost,solver,spot_options);
    
    u_bnd = [u_bnd;sol.eval(u_bnd_i)];
  end
end