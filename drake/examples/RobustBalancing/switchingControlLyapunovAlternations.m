function [ V,B ] = switchingControlLyapunovAlternations(x_mss,f,g,V0,B0,R)
if nargin < 6
  R = [];
end

N = 1;
nX = length(x_mss);
nU = size(g,2);
V = V0;
B = B0;
ndgrid_arg = mat2cell(repmat([-1;1],1,nU),2,ones(1,nU)');
[ugrid{1:nU}] = ndgrid(ndgrid_arg{:});
umat = zeros(2^nU,nU);
for i=1:nU,
  umat(:,i) = ugrid{i}(:);
end


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
  
  for j=1:2^nU
    Vdot = diff(V,x_mss)*(f + g*umat(j,:)')*(1+x_mss'*x_mss)^0;
    [prog, Vdot_sos,mult{j},coeff] = spotless_add_sprocedure(prog, -Vdot, 1-V,x_mss,4);
    for k=1:nU,
      [prog, Vdot_sos,bmult{j}{k},coeff] = spotless_add_sprocedure(prog, Vdot_sos, umat(j,k)*B(k),x_mss,4);
    end
    prog = prog.withSOS(Vdot_sos - gamma*(x_mss'*x_mss)^2);
  end
    
  spot_options = spotprog.defaultOptions;
  spot_options.verbose = true;
%   spot_options.do_fr = true;
  solver = @spot_mosek;
  sol = prog.minimize(-gamma,solver,spot_options);
  
  for j=1:2^nU
    mult{j} = sol.eval(mult{j});
    for k=1:nU,
      bmult{j}{k} = sol.eval(bmult{j}{k});
    end
  end
  
  %% iter 2
  prog = spotsosprog;
  prog = prog.withIndeterminate(x_mss);
  rho = 1;
  [prog,Q] = prog.newPSD(nX);
  V = x_mss'*Q*x_mss; 
%   V = V0;
  
  [prog,B] = prog.newFreePoly(monomials(x_mss,1:2),nU);
%   B = B0;

  if ~isempty(R)
    if iscell(R)
      for j=1:length(R),
        prog = prog.withPSD(Q-R{j});
      end
    else
      prog = prog.withPSD(Q-R);
    end
  end
  
  for j=1:2^nU
    Vdot = diff(V,x_mss)*(f + g*umat(j,:)')*(1+x_mss'*x_mss)^0;
    Vdot_sos = -Vdot - mult{j}*(rho-V);
    for k=1:nU,
      Vdot_sos = Vdot_sos - bmult{j}{k}*umat(j,k)*B(k);
    end
    prog = prog.withSOS(Vdot_sos + 1e-6*(x_mss'*x_mss));
  end    
  
  spot_options = spotprog.defaultOptions;
  spot_options.verbose = true;
%   spot_options.do_fr = true;
  solver = @spot_mosek;
  
  scale_mat = eye(length(x_mss));
  cost_coeffs = det(scale_mat*Q_init*scale_mat')*inv(scale_mat*Q_init*scale_mat');
  cost_coeffs = cost_coeffs/norm(cost_coeffs(:),inf);
  cost = sum(sum(scale_mat*(Q-Q_init)*scale_mat'.*cost_coeffs));
  
  cost = cost + 0.01*Q(1,1);
  
  sol = prog.minimize(cost,solver,spot_options);
  
  display(sprintf('Determinant from %f to %f, percent change %f',det(Q_init),det(sol.eval(Q)),100 - 100*det(sol.eval(Q))/det(Q_init)));
  V = sol.eval(V);
  B = sol.eval(B);
%   u = sol.eval(u)
  
end
end

