function [ V ] = alternateQuadraticLyapunovAlternations(x,f,V0)
% initialize

GO = V0;

N = 1;
nX = length(x);
for i=1:1,
  prog = spotsosprog;
  prog = prog.withIndeterminate(x);
  
  [prog,Q] = prog.newPSD(nX);
  V = x'*Q*x;  
  Vdot = diff(V,x)*f;
  
  [prog, Vdot_sos,Vdot_mult,coeff] = spotless_add_sprocedure(prog, -Vdot, 1-GO,x,4);
  prog = prog.withSOS(Vdot_sos + 0e-6*(1+x'*x));
  
  [prog, V_sos,V_mult,coeff] = spotless_add_eq_sprocedure(prog, V-1, 1-GO,x,4);
  prog = prog.withSOS(V_sos + 0e-6*(1+x'*x));
  
  spot_options = spotprog.defaultOptions;
  spot_options.verbose = true;
  spot_options.do_fr = true;
  solver = @spot_mosek;
  sol = prog.minimize(trace(Q),solver,spot_options);
  V = sol.eval(V);

%   keyboard
  
  V_mult = sol.eval(V_mult);
  Vdot_mult = sol.eval(Vdot_mult);
  
  prog = spotsosprog;
  [prog,gamma] = prog.newFree(1);
  prog = prog.withIndeterminate(x);
%   [prog,Q] = prog.newPSD(nX);
%   V = x'*Q*x;
  Vdot = diff(V,x)*f;
  
  [prog,G] = prog.newPSD(nX);
  GO = x'*G*x;
  
  prog = prog.withSOS(-Vdot - Vdot_mult*(1-GO) + 0e-6*(1+x'*x));
  prog = prog.withSOS(V - 1 - V_mult*(1-GO) + 0e-6*(1+x'*x));
  
  spot_options = spotprog.defaultOptions;
  spot_options.verbose = true;
  spot_options.do_fr = true;
  solver = @spot_mosek;
  sol = prog.minimize(trace(G),solver,spot_options);
  
%   keyboard
  
  solver = @spot_mosek;
%   prog = prog.withPos(1.03*sol.eval(trace(Q)) - trace(Q));
%   sol = prog.minimize(gamma,solver,spot_options);
  GO = sol.eval(GO);
  V = sol.eval(V);

  %   keyboard
  
end
end

