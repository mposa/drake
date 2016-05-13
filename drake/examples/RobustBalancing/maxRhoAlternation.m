function rho = maxRhoAlternation(x_mss,V,a,b,d)
  prog = spotsosprog;
  prog = prog.withIndeterminate(x_mss);
  [prog,rho] = prog.newFree(1);
  
  [prog,ab_sos] = spotless_add_sprocedure(prog, V-rho, -(2*a-b),x_mss,2);
  [prog,ab_sos] = spotless_add_sprocedure(prog, ab_sos, b,x_mss,2);
  prog = prog.withSOS(ab_sos);
  
  [prog,rad_sos] = spotless_add_sprocedure(prog, V-rho, -(b^2-4*a*d),x_mss,2);
  [prog,rad_sos] = spotless_add_sprocedure(prog, rad_sos, b,x_mss,2);
  prog = prog.withSOS(rad_sos);
  
  spot_options = spotprog.defaultOptions;
  spot_options.verbose = true;
  spot_options.do_fr = true;
  solver = @spot_mosek;
  sol = prog.minimize(-rho,solver,spot_options);  
  
  rho = sol.eval(rho);
end