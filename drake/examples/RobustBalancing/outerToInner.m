function outerToInner(model,T,R_diag,V,options)

rho = .7;

degree = options.degree;

% Program:
%  max  gamma
%  s.t  V > a ==> max_u Vdot > -a/T +    gamma/T

nU = model.num_inputs;

prog = spotsosprog;
[prog,x]=prog.newIndeterminate('x', model.num_states); % state
[prog,t]=prog.newIndeterminate('t', 1); % time


[f,g] = model.controlAffineDynamics(t, x);


q_monoms = monomials([t;x],0:degree);
[prog,q,q_coeff] = prog.newFreePoly(q_monoms,nU);
q_coeff = reshape(q_coeff,nU,[]);

% p(i)^2 <= (dV/dx*g_i(x))^2
% if p = c'*m(x) and (dV/dx*g_i(x))^2 = m'*Q*m
% Q - c*c' >= 0
%
% Schur complement: [Q c; c' 1] >= 0
% for i=1:nU,  
%   [prog, Q{i}] = prog.newPSD(length(q_monoms));
%   prog = addPolyEqs(prog,q_monoms'*Q{i}*q_monoms - (diff(V,x)*g(:,i))^2); %computationally not working that well
%   prog = prog.withPSD([Q{i} q_coeff(i,:)'; q_coeff(i,:) 1]);
%   
%   [var,pow,coeff] = decomp((diff(V,x)*g(:,i)));
%   
% end

% dV/dxg_i > 0 ==> q(i) < dV/dx*g_i
% dV/dxg_i < 0 ==> q(i) < -dV/dx*g_i
for i=1:nU,
  qipos_sos = (diff(V,x)*g(:,i) - q(i))*(1+[t;x]'*[t;x]);
  [prog, qipos_sos] = spotless_add_sprocedure(prog, qipos_sos, diff(V,x)*g(:,i),[t;x],degree);
  [prog, qipos_sos] = spotless_add_sprocedure(prog, qipos_sos, t*(T-t),[t;x],degree);
  [prog, qipos_sos] = spotless_add_sprocedure(prog, qipos_sos, V - rho,[t;x],degree);
  prog = prog.withSOS(qipos_sos);

  qineg_sos = (-diff(V,x)*g(:,i) - q(i))*(1+[t;x]'*[t;x]);
  [prog, qineg_sos] = spotless_add_sprocedure(prog, qineg_sos, -diff(V,x)*g(:,i),[t;x],degree);
  [prog, qineg_sos] = spotless_add_sprocedure(prog, qineg_sos, t*(T-t),[t;x],degree);
  [prog, qineg_sos] = spotless_add_sprocedure(prog, qineg_sos, V - rho,[t;x],degree);
  prog = prog.withSOS(qineg_sos);
end


% Program:
%  max  gamma
%  s.t  V > rho ==> max_u Vdot > rho/T -    gamma/T
[prog,gamma] = prog.newFree(1);

Vdot = diff(V,x)*f + sum(q);
Vdot_sos = Vdot + 1/T*(-gamma + rho);
Vdot_sos = Vdot_sos*(1+[t;x]'*[t;x]);
[prog, Vdot_sos] = spotless_add_sprocedure(prog, Vdot_sos, V - rho,[t;x],degree);
[prog, Vdot_sos] = spotless_add_sprocedure(prog, Vdot_sos, t*(T-t),[t;x],degree);

prog = prog.withSOS(Vdot_sos);

%% Solve
spot_options = spotprog.defaultOptions;
spot_options.verbose = true;
spot_options.do_fr = true;
solver = @spot_mosek;
% solver = @spot_sedumi;
sol = prog.minimize(-gamma,solver,spot_options);
% keyboard
end

function pr = addPolyEqs(pr,expr)
decvar = pr.variables;

A = diff(expr,decvar);
b = subs(expr,decvar,0*decvar);
[var,pow,Coeff] = decomp([b A].');

[pr,y] = pr.withEqs(Coeff'*[1;decvar]);
end