g = 10;
z_nom = 1;
R_diag = [1, 1, 1, 2, 2, 2];
% R_diag = [2, 1, 1, 4, 2, 2];
% f_min = .5;
% f_max = 1.5;
uz_bnd = .3;
ux_bnd = .1;
foot_radius = .05;
% inertia_ratio = .3^2/2; 
inertia_ratio = .6^2/2; 
step_max = .7;
step_time = .3;
z_inv_degree = 1;

model = LIPMVariation2D(g, z_nom, step_max, step_time, inertia_ratio, ux_bnd, uz_bnd, foot_radius, z_inv_degree);

%% Get an initial quadratic Lyapunov candidate
x = msspoly('x',model.num_states);
t = msspoly('t',1);
u = msspoly('u',model.num_inputs);
f = model.dynamics(t,x,u);

A = double(subs(diff(f,x),[t;x;u],zeros(1+model.num_states+model.num_inputs,1)));
B = double(subs(diff(f,u),[t;x;u],zeros(1+model.num_states+model.num_inputs,1)));

Q = 10*eye(model.num_states);
R = eye(model.num_inputs);
[K,Q] = lqr(A,B,Q,R);

%%
% A_state = diag([0;0;.5;pi/2;0;0]);
A_state{1} = diag([0;1/.5^2;0;0;0;0]);
A_state{2} = diag([0;0;1/(pi/2)^2;0;0;0]);
V0 = x'*Q*x;
[V,u_fn] = quadraticControlLyapunovAlternations(x,u,f,V0*100,A_state)
figure(1)
contourSpotless(V,x(1),x(4),[-1 1],[-2 2],[t;x([2;3;5;6])],zeros(model.num_states-1,1),1,{'r'});
%%
for i=1:20,
  [V,u_fn] = quadraticControlLyapunovAlternations(x,u,f,V,A_state);
%   figure(1)
  hold on
  if mod(i,2) == 0,
    contourSpotless(V,x(1),x(4),[-1 1],[-2 2],[t;x([2;3;5;6])],zeros(model.num_states-1,1),1,{'r'});
  else
    contourSpotless(V,x(1),x(4),[-1 1],[-2 2],[t;x([2;3;5;6])],zeros(model.num_states-1,1),1,{'b'});
  end
end;
% keyboard

%%
V_0step = V;
s = msspoly('s',model.num_reset_inputs);
r = model.reset(t,x,s);

% get inner approximation of pre-reset
% s = (c'SA)^-1 * c'Sx
c=double(diff(r,s));
Q_V0=double(diff(diff(V,x)',x))/2;

a = c'*Q_V0*c;
b = 2*x'*Q_V0*c;
d = x'*Q_V0*x - 1;
%%
rho2 = msspoly(1) + 10*t;
% V2 = x'*Q*x*100;
V2 = V_0step*10;
u2 = u_fn;
% [V2,u2,rho2] =  quadraticControlAlternationsWithResetNoGThreeSteps(x,u,f,V2,step_time,rho2,a,b,d);
figure(2)
%% construct alt-model
f_alt = subs(f,u(1),1);
u_alt = u(2:3);
constraint_alt = (x(4) + sqrt(model.gravity/model.z_nom)*x(1));

%%
for i=1:21,
  V2_i{i} = V2;
  u2_i{i} = u2;
  rho2_i{i} = rho2;
  [V2,u2,rho2] =  quadraticControlAlternationsWithResetNoGThreeSteps(x,u,f,V2,step_time,rho2,a,b,d,[],false);
% [V2,u2,rho2] =  quadraticControlAlternationsWithResetNoGThreeSteps(x,u_alt,f_alt,V2,step_time,rho2,a,b,d,constraint_alt);
%   figure(2)

  hold off
  contourSpotless(V,x(1),x(4),[-1 1],[-3 3],[t;x([2;3;5;6])],zeros(model.num_states-1,1),1,{'g'});
  hold on
  contourSpotless(V2,x(1),x(4),[-1 1],[-3 3],[t;x([2;3;5;6])],zeros(model.num_states-1,1),dmsubs(rho2,t,0),{'k'});
  contourSpotless(V2,x(1),x(4),[-1 1],[-3 3],[t;x([2;3;5;6])],[step_time;zeros(model.num_states-2,1)],dmsubs(rho2,t,step_time),{'r'});
  contourSpotless(b^2-4*a*d,x(1),x(4),[-1 1],[-3 3],[t;x([2;3;5;6])],[step_time;zeros(model.num_states-2,1)],0,{'y'});
  contourSpotless(b,x(1),x(4),[-1 1],[-3 3],[t;x([2;3;5;6])],[step_time;zeros(model.num_states-2,1)],0,{'b'});
  contourSpotless(2*a-b,x(1),x(4),[-1 1],[-3 3],[t;x([2;3;5;6])],[step_time;zeros(model.num_states-2,1)],0,{'b'});

  hold off
  
  sqrt(1./diag(.5*double(diff(diff(subs(V2,t,0),x)',x))))
end
%%
% s = msspoly('s',model.num_reset_inputs);
% r = model.reset(t,x,s);
% 
% % get inner approximation of pre-reset
% % s = (c'SA)^-1 * c'Sx
% c=double(diff(r,s));
% 
% % r = x + c*s! needs modification otherwise
% assert(isequal(r - c*s - x,msspoly(zeros(model.num_states,1)))); 
% 
% Q=double(diff(diff(V,x)',x))/2;
% T=eye(model.num_states) - c*inv(c'*Q*c)*c'*Q;
% % T=T'*Q*T;
% 
% % S < T
% % s = (c'Qc)^-1 * c'Qx > 1 AND x'Sx < 1 ==> (x+cc'Qc)'Q(x+cc'Qc) < 1
% 
% s_opt = -inv(c'*Q*c)*c'*Q*x;
% xp1 = x - c;
% 
% S = Q;
% 
% for i=1:0,
%   prog = spotsosprog();
%   prog = prog.withIndeterminate(x);
%   [prog,gamma] = prog.newPos(1);
% 
%   [prog, ssos] = spotless_add_sprocedure(prog, 1 - xp1'*Q*xp1, s_opt - 1,x_mss,2);
%   [prog, ssos,smult] = spotless_add_sprocedure(prog, ssos, 1 - x'*S*x,x_mss,2);
%   prog = prog.withSOS(ssos - gamma*(1+x'*x));
%   
%   spot_options = spotprog.defaultOptions;
%   spot_options.verbose = false;
%   spot_options.do_fr = true;
%   solver = @spot_mosek;
%   sol = prog.minimize(gamma,solver,spot_options);
%   
%   smult = sol.eval(smult);
%   
%   prog = spotsosprog();
%   prog = prog.withIndeterminate(x);
%   
%   [prog,S] = prog.newPSD(model.num_states);
%   prog = prog.withPSD(S - T'*Q*T);
%   
%   [prog, ssos,smult2] = spotless_add_sprocedure(prog, 1 - xp1'*Q*xp1 - smult*(1-x'*S*x), s_opt - 1,x_mss,2);
%   prog = prog.withSOS(ssos);
%   
%   spot_options = spotprog.defaultOptions;
%   spot_options.verbose = true;
%   spot_options.do_fr = true;
%   solver = @spot_mosek;
%   sol = prog.minimize(trace(S),solver,spot_options);
%   S = sol.eval(S);
% end
% % Vm(x) <= 1 ==> x'Tx <= 1

%%
% Qlqr = 100*eye(model.num_states);
% Rlqr = eye(model.num_inputs);
% [K,Slqr] = lqr(A,B,Qlqr,Rlqr);
% % S2 = Slqr*100;
% S2 = Q*10;
% rho2 = msspoly(1) + 20*t + 20*t^2;
% V2 = x'*S2*x;
% 
% a = c'*Q*c;
% b = 2*x'*Q*c;
% d = x'*Q*x - 1;
% % [gmult_v,gmult_vdot,gmult_rad,gmult_ab] = initQuadraticControlAlternationWithReset(x,u,f,V2,step_time,rho2,a,b,d);

%%
% for i=1:1,
%   [V2,u2,rho2,gmult_v,gmult_vdot,gmult_rad,gmult_ab] =  quadraticControlAlternationsWithReset(x,u,f,V2,step_time,rho2,a,b,d,gmult_v,gmult_vdot,gmult_rad,gmult_ab);
%   figure(2)
%   hold off
%   contourSpotless(V_inner,x(1),x(4),[-2 2],[-2 2],[t;x([2;3;5;6])],zeros(model.num_states-1,1),1,{'g'});
%   hold on
%   contourSpotless(V2,x(1),x(4),[-2 2],[-2 2],[t;x([2;3;5;6])],zeros(model.num_states-1,1),dmsubs(rho2,t,0),{'k'});
%   contourSpotless(V2,x(1),x(4),[-2 2],[-2 2],[t;x([2;3;5;6])],[step_time;zeros(model.num_states-2,1)],dmsubs(rho2,t,step_time),{'r'});
%   hold off
% end