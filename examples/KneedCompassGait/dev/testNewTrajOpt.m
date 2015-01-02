function [p,xtraj,utraj,ltraj,ljltraj,z,F,info,traj_opt] = testNewTrajOpt(xtraj,utraj,ltraj,ljltraj,scale)
warning('off','Drake:RigidBodyManipulator:UnsupportedContactPoints');
warning('off','Drake:RigidBodyManipulator:WeldedLinkInd');
options.terrain = RigidBodyFlatTerrain();
options.floating = true;
options.ignore_self_collisions = true;
p = PlanarRigidBodyManipulator('../KneedCompassGait.urdf',options);
% trajopt = ContactImplicitTrajectoryOptimization(p,[],[],[],10,[1 1]);

%todo: add joint limits, periodicity constraint

N = 15;
T = 1.5;
T0 = 1;

% periodic constraint
R_periodic = zeros(p.getNumStates,2*p.getNumStates);
R_periodic(2,2) = 1; %z
R_periodic(3,3) = 1; %pitch-hip w/symmetry
R_periodic(3,5) = 1; %pitch-hip w/symmetry
R_periodic(4,6) = 1; %knee w/symmetry
R_periodic(6,4) = 1; %knee w/symmetry
R_periodic(5,5) = -1; %hip w/symmetry

R_periodic(7,7) = 1; %x-vel
R_periodic(8,8) = 1; %z-vel
R_periodic(9,9) = 1; %pitch-hip w/symmetry
R_periodic(9,11) = 1; %pitch-hip w/symmetry
R_periodic(10,12) = 1; %knee w/symmetry
R_periodic(12,10) = 1; %knee w/symmetry
R_periodic(11,11) = -1; %hip w/symmetry

R_periodic(2:end,p.getNumStates+2:end) = -eye(p.getNumStates-1);

periodic_constraint = LinearConstraint(zeros(p.getNumStates,1),zeros(p.getNumStates,1),R_periodic);

% x0 = [0;0;1;zeros(15,1)];
% xf = [0;0;1;zeros(15,1)];
x0 = [0;1;0;0;0;1;.4;zeros(5,1)];
xf = [.2;1;0;1;0;0;.4;zeros(5,1)];

N2 = floor(N/2);

if nargin < 2
  %Try to come up with a reasonable trajectory
  x1 = [.1;1;pi/8-pi/16;pi/8;-pi/8;pi/8;zeros(6,1)];
  t_init = linspace(0,T0,N);
%   traj_init.x = PPTrajectory(foh(t_init,linspacevec(x0,xf,N)));
  traj_init.x = PPTrajectory(foh(t_init,[linspacevec(x0,x1,N2), linspacevec(x1,xf,N-N2)]));
  traj_init.u = PPTrajectory(foh(t_init,zeros(3,N)));
  traj_init.l = PPTrajectory(foh(t_init,[repmat([1;zeros(7,1)],1,N2) repmat([zeros(4,1);1;zeros(3,1)],1,N-N2)]));
  traj_init.ljl = PPTrajectory(foh(t_init,zeros(p.getNumJointLimitConstraints,N)));
  scale = 1;
else
  t_init = xtraj.pp.breaks;
  traj_init.x = xtraj;
  traj_init.u = utraj;
  traj_init.l = ltraj;
  traj_init.ljl = ljltraj;
end
T_span = [.2 T];


x0_min = [x0(1:5);.5; -inf; 0; -inf(4,1)];
x0_max = [x0(1:5);inf;  inf; 0; inf(4,1)];
xf_min = [.4;-inf(11,1)];
xf_max = inf(12,1);

to_options.compl_slack = scale*.1;
to_options.lincompl_slack = scale*.001;
to_options.jlcompl_slack = scale*.01;

to_options.nlcc_mode = 2;
to_options.lincc_mode = 1;
to_options.lambda_mult = p.getMass*9.81*T0/N;
to_options.lambda_jl_mult = T0/N;

to_options.integration_method = ContactImplicitTrajectoryOptimization.MIXED;

traj_opt = ContactImplicitTrajectoryOptimization(p,N,T_span,to_options);
traj_opt = traj_opt.addRunningCost(@running_cost_fun);
traj_opt = traj_opt.addRunningCost(@foot_height_fun);
traj_opt = traj_opt.addFinalCost(@final_cost_fun);
traj_opt = traj_opt.addStateConstraint(BoundingBoxConstraint(x0_min,x0_max),1);
traj_opt = traj_opt.addStateConstraint(BoundingBoxConstraint(xf_min,xf_max),N);
traj_opt = traj_opt.addStateConstraint(periodic_constraint,{[1 N]});

% min force
min_force = LinearConstraint(.25*ones(N,1),inf(N,1),kron(eye(N),[1 1]));
traj_opt = traj_opt.addConstraint(min_force,reshape(traj_opt.l_inds(1:traj_opt.nD+2:end,:),[],1));

xr = randn(12,1);
ur = randn(3,1);
[f,df] = foot_height_fun(0,xr,ur);
[f2,df2] = geval(@foot_height_fun,0,xr,ur,struct('grad_method','numerical'));

l1 = traj_opt.l_inds(5:end,1:5);
l2 = traj_opt.l_inds(1:4,end-4:end);
% traj_opt = traj_opt.addConstraint(ConstantConstraint(0*l1(:)),l1);
% traj_opt = traj_opt.addConstraint(ConstantConstraint(0*l2(:)),l2);

% traj_opt = traj_opt.setCheckGrad(true);
snprint('snopt.out');
traj_opt = traj_opt.setSolverOptions('snopt','MajorOptimalityTolerance',1e-5);
traj_opt = traj_opt.setSolverOptions('snopt','MajorIterationsLimit',200);
traj_opt = traj_opt.setSolverOptions('snopt','MinorIterationsLimit',200000);
traj_opt = traj_opt.setSolverOptions('snopt','IterationsLimit',200000);

% traj_opt = traj_opt.setSolverOptions('snopt','New Basis File',1);
% traj_opt = traj_opt.setSolverOptions('snopt','Old Basis File',0);


% prevent knee from bending backwards
knee_constraint = LinearConstraint([-inf(3,1);0;-inf;0;-inf(6,1)],inf(12,1),eye(12));
traj_opt = traj_opt.addStateConstraint(knee_constraint,mat2cell(1:N,1,ones(N,1)));

% contact delta fun
lz_inds = traj_opt.l_inds(1:4:end,:);
contactDeltaCost = FunctionHandleConstraint(-inf,inf,size(lz_inds,2),@contact_delta_cost_fun,1);
traj_opt = traj_opt.addCost(contactDeltaCost,lz_inds(1,:)');
traj_opt = traj_opt.addCost(contactDeltaCost,lz_inds(2,:)');

[xtraj,utraj,ltraj,ljltraj,z,F,info] = traj_opt.solveTraj(t_init,traj_init);

function [f,df] = running_cost_fun(h,x,u)
  f = u'*u;
  df = [0 zeros(1,12) 2*u'];
end

function [f,df] = foot_height_fun(h,x,u)
  q = x(1:6);
  K = 100;
  [phi,~,~,~,~,~,~,~,n] = p.contactConstraints(q,false,struct('terrain_only',false));
  phi0 = [.1;.1];
  f = K*(phi - phi0)'*(phi - phi0);
  % phi: 2x1
  % n: 2xnq
  df = [0 2*K*(phi-phi0)'*n zeros(1,9)];
end

  function [f,df] = contact_delta_cost_fun(l)
    K = 10;
    f = K*sum(diff(l).^2);
    n = length(l);
    R = sparse([1:n-1 2:n], [1:n-1 1:n-1], [ones(n-1,1);-ones(n-1,1)]);
    df = -2*K*diff(l)'*R';
  end

function [f,df] = final_cost_fun(T,x)
  K = 100;
  f = K*T;
  df = [K zeros(1,12)];
end

end
