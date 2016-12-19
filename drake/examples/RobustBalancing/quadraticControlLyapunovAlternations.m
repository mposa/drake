function [ V,u ] = quadraticControlLyapunovAlternations(x_mss,u_mss,f,V0, R)
N = 1;
nX = length(x_mss);
V = V0;

if nargin < 5
  R = [];
end
Q_init = double(diff(diff(V,x_mss)',x_mss))/2;

for i=1:1,
  %% iter 1
  prog = spotsosprog;
  prog = prog.withIndeterminate(x_mss);
  [prog,gamma] = prog.newPos(1);
  prog = prog.withPos(1 - gamma);
  [prog,u] = prog.newFreePoly(monomials(x_mss,1:3),length(u_mss));
  
  Vdot = diff(V,x_mss)*subs(f,u_mss,u);
  
  [prog, Vdot_sos,mult,coeff] = spotless_add_sprocedure(prog, -Vdot, 1-V,x_mss,4,2);

  prog = prog.withSOS(Vdot_sos - gamma*(x_mss'*x_mss)^3);
  for i=1:length(u),
    [prog, u_sos,upmult{i}] = spotless_add_sprocedure(prog, 1-u(i), 1-V,x_mss,2);
    prog = prog.withSOS(u_sos);
    [prog, u_sos,ummult{i}] = spotless_add_sprocedure(prog, u(i)+1, 1-V,x_mss,2);
    prog = prog.withSOS(u_sos);
  end
  
  
  spot_options = spotprog.defaultOptions;
  spot_options.verbose = true;
%   spot_options.do_fr = true;
  solver = @spot_mosek;
  sol = prog.minimize(-gamma,solver,spot_options);
  u = sol.eval(u)
  
  for i=1:length(u),
    upmult{i} = sol.eval(upmult{i});
    ummult{i} = sol.eval(ummult{i});
  end
%   keyboard
  
  mult = sol.eval(mult);
  
  %% iter 2
  prog = spotsosprog;
  prog = prog.withIndeterminate(x_mss);
  [prog,rho] = prog.newPos(1);
  
  [prog,u] = prog.newFreePoly(monomials(x_mss,1:3),length(u_mss));
  
  Vdot = diff(V,x_mss)*subs(f,u_mss,u);
  
  Vdot_sos = -Vdot - mult*(rho-V);

  prog = prog.withSOS(Vdot_sos + 1e-6*(x_mss'*x_mss)^3);
  for i=1:length(u),
    prog = prog.withSOS(1-u(i) - upmult{i}*(rho-V));
    prog = prog.withSOS(u(i)+1 - ummult{i}*(rho-V));
  end
  
  
  spot_options = spotprog.defaultOptions;
  spot_options.verbose = true;
%   spot_options.do_fr = true;
  solver = @spot_mosek;
  sol = prog.minimize(-rho,solver,spot_options);
  u = sol.eval(u)
  
  %% iter 3
  prog = spotsosprog;
  prog = prog.withIndeterminate(x_mss);
%   [prog,gamma] = prog.newPos(1);
  [prog,Q] = prog.newPSD(nX);
  V = x_mss'*Q*x_mss;
  Vdot = diff(V,x_mss)*subs(f,u_mss,u);
  
  Vdot_sos = -Vdot + mult*(V-1) + (x_mss'*x_mss)^3*1e-6;
  prog = prog.withSOS(Vdot_sos);
  
  for i=1:length(u),
    prog = prog.withSOS(upmult{i}*(V-1) + 1 - u(i));
    prog = prog.withSOS(ummult{i}*(V-1) + 1 + u(i));
  end
  
  if ~isempty(R)
    if iscell(R)
      for j=1:length(R),
        prog = prog.withPSD(Q-R{j});
      end
    else
      prog = prog.withPSD(Q-R);
    end
  end
  
  spot_options = spotprog.defaultOptions;
  spot_options.verbose = true;
%   spot_options.do_fr = true;
  solver = @spot_mosek;
%   cost = trace(Q) + 50*Q(1,1);
  % linearization of determinant
  scale_mat = eye(length(x_mss));
  cost_coeffs = det(scale_mat*Q_init*scale_mat')*inv(scale_mat*Q_init*scale_mat');
  cost_coeffs = cost_coeffs/norm(cost_coeffs(:),inf);
  cost = sum(sum(scale_mat*(Q-Q_init)*scale_mat'.*cost_coeffs));

  sol = prog.minimize(cost,solver,spot_options);
%   
%   keyboard
  
%   solver = @spot_mosek;
%   prog = prog.withPos(1.03*sol.eval(trace(Q)) - trace(Q));
%   sol = prog.minimize(gamma,solver,spot_options);
  
  V = sol.eval(V);
  display(sprintf('Determinant from %f to %f, percent change %f',det(Q_init),det(sol.eval(Q)),100 - 100*det(sol.eval(Q))/det(Q_init)));

  %   keyboard
  
end
end

