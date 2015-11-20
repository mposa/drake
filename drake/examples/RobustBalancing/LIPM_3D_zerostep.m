prog = spotsosprog;
degree = 6;

[prog,t]=prog.newIndeterminate('t',1);
[prog,q]=prog.newIndeterminate('q',2);
[prog,v]=prog.newIndeterminate('v',2);


do_backoff = false;
time_varying = true;


% problem data
T = 2;  %step time
g = 10; %gravity
z_nom = 1; %nominal height
goal_radius = .01;
R_diag = [2 2 2 2];  % state space ball

x = [q;v];

if time_varying
  V_vars = [t;x];
else
  V_vars = x;
end
W_vars = x;

[prog,V] = prog.newFreePoly(monomials(V_vars,0:degree));
[prog,W] = prog.newFreePoly(monomials(W_vars,0:degree));

f = [v;q*g/z_nom];

% tau = t/T
% dx/dtau = dx/dt * dt/dtau = dx/dt*T

% time rescaling
f = f*T;
T = 1;


Vdot = diff(V,x)*f + diff(V,t);

sos = [V;-Vdot; W;W - subs(V,t,0) - 1];

A = diag(1./(R_diag.^2));

h_goal = goal_radius^2 - x'*x;
h_X = 1 - x'*A*x;

[prog, sos(1)] = spotless_add_sprocedure(prog, sos(1), h_goal,V_vars,degree-2);
[prog, sos(2)] = spotless_add_sprocedure(prog, sos(2), h_X,V_vars,degree-2);
[prog, sos(3)] = spotless_add_sprocedure(prog, sos(3), h_X,W_vars,degree-2);
[prog, sos(4)] = spotless_add_sprocedure(prog, sos(4), h_X,W_vars,degree-2);

if time_varying
  [prog, sos(1)] = spotless_add_sprocedure(prog, sos(1), T^2-t^2,V_vars,degree-2);
  [prog, sos(2)] = spotless_add_sprocedure(prog, sos(2), T^2-t^2,V_vars,degree-2);
end

%% Setup cost function
[vars,alphas,coeff] = decomp(W, prog.freeVar);
inds = [1;3;2;4];
alphas = alphas(:,inds);
assert(isequal(vars(inds),x))

nX = 2;

% [~,alphas] = monomials(x_vars,0:degree);
betas = 0.5*(alphas + 1);
Ra = (1.^(sum(alphas,2) + nX))./(sum(alphas,2) + nX);
IS = 2*prod(gamma(betas),2)./(gamma(sum(betas,2)));
l = Ra.*IS;
alphaszero = (mod(alphas,2) ~= 0);
alphaszero = any(alphaszero,2);
l(alphaszero) = 0;
l = l.*prod(repmat(R_diag,size(alphas,1),1).^(alphas+1),2);

for i=1:length(sos)
  prog = prog.withSOS(sos(i));
end

options = spotprog.defaultOptions;
options.verbose = true;
options.do_fr = true;
sol = prog.minimize(coeff*l,@spot_mosek,options);

if do_backoff
  
  prog = prog.withPos(sol.eval(coeff*l)*1.01 - coeff*l);
  sol = prog.minimize(0,@spot_mosek,options);
end

%% plotting
Vsol = sol.eval(V);
Wsol = sol.eval(W);

vars_elim = [q(2);v(2)];
vars_elim_val = [0;0];

Vplot = subs(Vsol,vars_elim,vars_elim_val);
Wplot = subs(Wsol,vars_elim,vars_elim_val);
h_goalplot = subs(h_goal,vars_elim,vars_elim_val);
h_Xplot = subs(h_X,vars_elim,vars_elim_val);


vars_plot = [q(1);v(1)];

[AXIS1,AXIS2] = meshgrid(linspace(-R_diag(1),R_diag(1),400),linspace(-R_diag(2),R_diag(2),400));
h_goal_val = reshape(msubs(h_goalplot,vars_plot,[AXIS1(:)';AXIS2(:)']),size(AXIS1,1),[]);
h_X_val = reshape(msubs(h_Xplot,vars_plot,[AXIS1(:)';AXIS2(:)']),size(AXIS1,1),[]);

Wval = reshape(msubs(Wplot,vars_plot,[AXIS1(:)';AXIS2(:)']),size(AXIS1,1),[]);
figure(1)
hold off
[cl,h]=contour(AXIS1,AXIS2,Wval,[1 1],'b');
clabel(cl,h);
hold on
[cl,h]=contour(AXIS1,AXIS2,h_goal_val,[0 0],'k');
[cl,h]=contour(AXIS1,AXIS2,h_X_val,[0 0],'r');

Vval = reshape(msubs(subs(Vplot,t,0),vars_plot,[AXIS1(:)';AXIS2(:)']),size(AXIS1,1),[]);
figure(2)
hold off
[cl,h]=contour(AXIS1,AXIS2,Vval,[0 0],'b');
clabel(cl,h);
hold on
[cl,h]=contour(AXIS1,AXIS2,h_goal_val,[0 0],'k');
[cl,h]=contour(AXIS1,AXIS2,h_X_val,[0 0],'r');


%% sim test
% x0 = [-.701;.7];
% dt = 1e-4;
% t = 0:dt:5;
% x_vec = zeros(2,length(t));
% x_vec(:,1) = x0;
% for i=2:length(t),
%   x_vec(:,i) = [1 dt;dt 1]*x_vec(:,i-1);
% end
% figure(3)
% plot(t,x_vec)

%%
V0 = Vsol;
save V0_LIPM V0
