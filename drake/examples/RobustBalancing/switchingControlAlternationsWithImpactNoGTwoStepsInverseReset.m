function [ V,B,rho] = switchingControlAlternationsWithImpactNoGTwoStepsInverseReset(x_mss,s_mss,f,g,V0,rho0,B0,T,r_inv,reset_constraint,Vprev)
% variation on altbounds, but using the inverse reset map
do_backoff = false;
backoff_ratio = .08;

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

V = V0;

rho = rho0;
t = msspoly('t',1);
Q_init = double(diff(diff(subs(V,t,0),x_mss)',x_mss))/2;
Q_init_T = double(diff(diff(subs(V,t,T),x_mss)',x_mss))/2;
for i=1:1,
  
    
%   %% step 0
%   prog = spotsosprog;
%   prog = prog.withIndeterminate(x_mss);
%   prog = prog.withIndeterminate(s_mss);
%   [prog,Q] = prog.newPSD(nX);
%   [prog,c] = prog.newFree(nX);
%   V = x_mss'*Q*x_mss + c'*x_mss;
%  
%   ds = diff(reset_constraint,s_mss);
%   ns = (reset_constraint - s_mss*ds);
%   rsub = diff(r,s_mss)*ns + (r - s_mss*diff(r,s_mss))*ds;
%     
% %   fr = subs(Vprev,x_mss,rsub) - ds^2;
% %   figure(1)
% %   contourSpotless(fr,x_mss(1),x_mss(3),[-1 1],[-3 3],x_mss([2;4;5]),zeros(3,1),0,{'g'});
%   
%   reset_sos = (V-1)*([s_mss;x_mss]'*[s_mss;x_mss])^3;
%   Vprev_reset = subs(Vprev,x_mss,r);
% %   [prog, reset_sos] = spotless_add_sprocedure(prog, reset_sos,Vprev - 1,[s_mss;x_mss],4);
%   [prog, reset_sos] = spotless_add_sprocedure(prog, reset_sos,Vprev_reset - 1,[s_mss;x_mss],4);
%   [prog, reset_sos] = spotless_add_eq_sprocedure(prog, reset_sos,reset_constraint,[s_mss;x_mss],5);
%   
% %   reset_sos = (V-1)*(1 + [x_mss]'*[x_mss])^3;
% %   fr = x_mss'*x_mss-1;
% %   [prog, reseet_sos] = spotless_add_sprocedure(prog, reset_sos,fr,[x_mss],4);
%   
% %   reset_sos = subs(reset_sos,x_mss([2;4;5]),zeros(3,1));
% %   reset_sos = subs(reset_sos,s_mss,1);
%   prog = prog.withSOS(reset_sos);
% 
%   
%   spot_options = spotprog.defaultOptions;
%   spot_options.verbose = true;
%   spot_options.do_fr = false;
%   spot_options.sos_slack = 0e-6;
%   solver = @spot_mosek;
%   
% %   cost = Q(1,1);
%   cost = trace(Q);
%   sol = prog.minimize(cost,solver,spot_options);
%   
%   figure(1)
%   contourSpotless(sol.eval(V),x_mss(1),x_mss(3),[-1 1],[-3 3],x_mss([2;4;5]),zeros(3,1),1,{'g'});
  
  
  %% step 1
  prog = spotsosprog;
  prog = prog.withIndeterminate(x_mss);
  prog = prog.withIndeterminate(t);
  sproc_vars = [t;x_mss];

  [prog,gamma] = prog.newPos(1);
  prog = prog.withPos(1 - gamma);
  
  rhodot = diff(rho,t);  
  
  for j=1:2^nU
    Vdot = diff(V,x_mss)*(f + g*umat(j,:)') + diff(V,t);
    
    [prog, Vdot_sos,mult{j},coeff] = spotless_add_eq_sprocedure(prog, rhodot-Vdot, rho-V,x_mss,2);
    [prog, Vdot_sos] = spotless_add_sprocedure(prog, Vdot_sos, t*(T-t),sproc_vars,4);
    for k=1:nU,
      [prog, Vdot_sos,bmult{j}{k},coeff] = spotless_add_sprocedure(prog, Vdot_sos, umat(j,k)*B(k),x_mss,4);
    end
    Vdot_degree = even_degree(Vdot_sos,sproc_vars)-2;
    prog = prog.withSOS(Vdot_sos - gamma*(x_mss'*x_mss)^(Vdot_degree/2));
  end  
  
%   rad_sos = subs(V-rho,t,T)*(1+x_mss'*x_mss);
%   [prog, rad_sos] = spotless_add_sprocedure(prog, rad_sos, 4*a*d-b^2,x_mss,2);
%   [prog, rad_sos] = spotless_add_sprocedure(prog, rad_sos, b,x_mss,2);
%   prog = prog.withSOS(rad_sos);
%   
%   ab_sos = subs(V-rho,t,T)*(1+x_mss'*x_mss);
%   [prog, ab_sos] = spotless_add_sprocedure(prog, ab_sos, 2*a-b,x_mss,2);
%   [prog, ab_sos] = spotless_add_sprocedure(prog, ab_sos, -a+b-d,x_mss,2);
%   [prog, ab_sos] = spotless_add_sprocedure(prog, ab_sos, -b,x_mss,2);
%   prog = prog.withSOS(ab_sos);  

  spot_options = spotprog.defaultOptions;
  spot_options.verbose = true;
  spot_options.do_fr = false;
  spot_options.clean_primal = true;
  solver = @spot_mosek;
  sol = prog.minimize(-gamma,solver,spot_options);
  
  for j=1:2^nU
    mult{j} = sol.eval(mult{j});
    for k=1:nU,
      bmult{j}{k} = sol.eval(bmult{j}{k});
    end
  end

  

  
  %% step 2 of 2
  prog = spotsosprog;
  prog = prog.withIndeterminate(x_mss);
  prog = prog.withIndeterminate(t);
  prog = prog.withIndeterminate(s_mss);
  
%   [prog,gamma] = prog.newPos(1);

  rho_T_nom = double(subs(rho,t,T));

  [prog,rho] = prog.newFreePoly(monomials(t,1:2));
  rho = rho  + 1; % set rho(0) = 1

  [prog,V] = prog.newFreePoly(reshape(monomials(x_mss,0:2)*monomials(t,(0:1))',[],1));
  [prog,B] = prog.newFreePoly(monomials(x_mss,1:2),nU);
% V = V0
  S0 = diff(diff(subs(V,t,0),x_mss)',x_mss)/2;
  ST = diff(diff(subs(V,t,T),x_mss)',x_mss)/2;
%   [prog,Q] = prog.newPSD(nX);
%   V = x_mss'*Q*x_mss;
%   S0 = Q;
%   ST = Q;

  % V >= 0
%   [prog,V_pos] = spotless_add_sprocedure(prog, V,t*(T-t),[t;x_mss],4);
%   prog = prog.withSOS(V_pos);

  prog = prog.withSOS(subs(V,t,0));
  prog = prog.withSOS(subs(V,t,T));

%   prog = prog.withPSD(S0);
%   prog = prog.withPSD(ST);
  
  rhodot = diff(rho,t);
    
  for j=1:2^nU
    Vdot = diff(V,x_mss)*(f + g*umat(j,:)') + diff(V,t);    Vdot_sos = rhodot-Vdot - mult{j}*(rho-V);
    for k=1:nU,
      Vdot_sos = Vdot_sos - bmult{j}{k}*umat(j,k)*B(k);
    end
    [prog, Vdot_sos] = spotless_add_sprocedure(prog, Vdot_sos, t*(T-t),sproc_vars,4);
    Vdot_degree = even_degree(Vdot_sos,sproc_vars);
    prog = prog.withSOS(Vdot_sos + 1e-6*(x_mss'*x_mss)^(Vdot_degree/2));
  end            
 
  
  V_inv_reset = subs(V,x_mss,r_inv);
  reset_sos = subs(V_inv_reset - rho,t,T)*(1 + [s_mss;x_mss]'*[s_mss;x_mss])^2;
%   [prog, reset_sos] = spotless_add_sprocedure(prog, reset_sos,Vprev - 1,[s_mss;x_mss],4);
  [prog, reset_sos] = spotless_add_sprocedure(prog, reset_sos,Vprev - 1,[s_mss;x_mss],4);
  [prog, reset_sos] = spotless_add_eq_sprocedure(prog, reset_sos,reset_constraint,[s_mss;x_mss],3);
  prog = prog.withSOS(reset_sos);
  
%   spot_options = spotprog.defaultOptions;
%   spot_options.verbose = true;
%   spot_options.do_fr = false;
  spot_options.sos_slack = 0e-6;
  solver = @spot_mosek;
  
%   cost = 0.01*trace(S0)+5*S0(2,2);

  scale_mat = eye(length(x_mss));
  scale_mat(1) = 1;
%   scale_mat(5) = 0;
  
  det_init = det(scale_mat*Q_init*scale_mat');
  det_init_T = det(scale_mat*Q_init_T*scale_mat');

  % linearization of determinant
  cost_coeffs = det_init*inv(scale_mat*Q_init*scale_mat');
  cost_coeffs_T = det_init_T*inv(scale_mat*Q_init_T*scale_mat');
  
  cost_rho = -length(x_mss)/rho_T_nom*det_init_T*subs(rho,t,T); 
   
  cost = 1*sum(sum(scale_mat*(S0-Q_init)*scale_mat'.*cost_coeffs));
  cost = cost + 1*rho_T_nom^(-length(x_mss))*(sum(sum(scale_mat*(ST-Q_init_T)*scale_mat'.*cost_coeffs_T)) + cost_rho);
  
  cost = cost/norm(cost_coeffs(:),inf);
  
%   xstar = [sqrt(10);0;0;-1;zeros(2,1)];
%   cost = cost + .01*xstar'*S0*xstar;
  
%   cost = cost + 1*S0(1,1);
  
  sol = prog.minimize(cost,solver,spot_options);
  
  if do_backoff
    display('Backing off in iter 2 and resolving');
    prog = prog.withPos(cost + backoff_ratio*abs(double(sol.eval(cost))));
    sol = prog.minimize(0,solver,spot_options);
  end
  
 
  det_new = det(scale_mat*double(sol.eval(S0))*scale_mat');
  det_new_T = det(scale_mat*double(sol.eval(ST))*scale_mat');
  display(sprintf('Determinant from %f to %f, percent change %f',det_init,det_new,100-100*det_new/det_init));
  display(sprintf('T-determinant from %f to %f, percent change %f',det_init_T,det_new_T,100-100*det_new_T/det_init_T));
  
  V = sol.eval(V);
  rho = sol.eval(rho);
  B = sol.eval(B);  
end

