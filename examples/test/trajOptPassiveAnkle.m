function [p,xtraj,utraj,ltraj,ljltraj,z,F,info,traj_opt] = trajOptPassiveAnkle(xtraj,utraj,ltraj,ljltraj,scale)
warning('off','Drake:RigidBodyManipulator:UnsupportedContactPoints');
warning('off','Drake:RigidBodyManipulator:WeldedLinkInd');
warning('off','Drake:RigidBodyManipulator:UnsupportedJointLimits');
options.terrain = RigidBodyFlatTerrain();
options.floating = true;
options.ignore_self_collisions = true;
p = PlanarRigidBodyManipulator('OneLegHopper_passiveankle.urdf',options);
% trajopt = ContactImplicitTrajectoryOptimization(p,[],[],[],10,[1 1]);

%todo: add joint limits, periodicity constraint

N = 30;

q0 = [0;.7;.2;-.4;.2+pi/2];
x0 = [q0;zeros(5,1)];

qf = [0;0;.6;-1.2;.6+pi/2];
phi_f = p.contactConstraints(qf);
qf(2) = -phi_f(1);

N1 = floor(N/2);
N2 = N-N1;
d = floor(N/4);
tf0 = .5;
if nargin < 2
  t_init = linspace(0,tf0,N);
  traj_init.x = PPTrajectory(foh(t_init,linspacevec([q0;zeros(5,1)],[qf;zeros(5,1)],N)));
  traj_init.x = traj_init.x.setOutputFrame(p.getStateFrame);
 
  traj_init.u = PPTrajectory(foh(t_init,randn(1,N)));
  
  lp = [1;0;0;0];
  ln = zeros(4,1);
  traj_init.l = PPTrajectory(foh(t_init,[repmat([ln;ln],1,N1-d) repmat([lp;lp],1,N2+d)]));
  traj_init.ljl = [];
  
  scale = 1;
else
  t_init = xtraj.pp.breaks;
  if length(t_init) ~= N
    t_init = linspace(0,t_init(end),N);
  end
  traj_init.x = xtraj;
  traj_init.u = utraj;
  traj_init.l = ltraj;
  traj_init.ljl = ljltraj;
end


T_span = [1 1];


x0_min = [q0;zeros(5,1)];
x0_max = [q0;zeros(5,1)];

xf_min = [qf;zeros(5,1)] - [.1;.01;zeros(8,1)];
xf_max = [qf;zeros(5,1)] + [.1;.01;zeros(8,1)];

to_options.compl_slack = scale*.01;
to_options.lincompl_slack = scale*.001;
to_options.jlcompl_slack = scale*.01;

to_options.nlcc_mode = 2;
to_options.lincc_mode = 1;
to_options.lambda_mult = p.getMass*9.81*tf0/N/2;
to_options.lambda_jl_mult = tf0/N;

% to_options.integration_method = ContactImplicitTrajectoryOptimization.MIDPOINT;
to_options.integration_method = ContactImplicitTrajectoryOptimization.BACKWARD_EULER;

traj_opt = ContactImplicitTrajectoryOptimization(p,N,T_span,to_options);
traj_opt = traj_opt.addRunningCost(@running_cost_fun);
% traj_opt = traj_opt.addRunningCost(@foot_height_fun);
% traj_opt = traj_opt.addFinalCost(@final_cost_fun);
traj_opt = traj_opt.addStateConstraint(BoundingBoxConstraint(x0_min,x0_max),1);
traj_opt = traj_opt.addStateConstraint(BoundingBoxConstraint(xf_min,xf_max),N);

% traj_opt = traj_opt.setCheckGrad(true);
snprint('snopt.out');
traj_opt = traj_opt.setSolverOptions('snopt','MajorIterationsLimit',100);
traj_opt = traj_opt.setSolverOptions('snopt','MinorIterationsLimit',200000);
traj_opt = traj_opt.setSolverOptions('snopt','IterationsLimit',200000);
traj_opt = traj_opt.setSolverOptions('snopt','SuperbasicsLimit',5000);
traj_opt = traj_opt.setSolverOptions('snopt','MajorOptimalityTolerance',1e-5);
[xtraj,utraj,ltraj,ljltraj,z,F,info] = traj_opt.solveTraj(t_init,traj_init);

  function [f,df] = contact_delta_cost_fun(l)
    K = 500;
    f = K*sum(diff(l).^2);
    n = length(l);
    R = sparse([1:n-1 2:n], [1:n-1 1:n-1], [ones(n-1,1);-ones(n-1,1)]);
    df = -2*K*diff(l)'*R';
  end

  function [f,df] = running_cost_fun(h,x,u)
    K = 10;
    R = eye(1);
    f = K*u'*R*u;
    df = [0 zeros(1,10) 2*K*u'*R];
  end

  function [f,df] = foot_height_fun(h,x,u)
    q = x(1:9);
    K = 50;
    [phi,~,~,~,~,~,~,~,n] = p.contactConstraints(q,false,struct('terrain_only',false));
    phi0 = [.1;.2;.1;.2];
    f = K*(phi - phi0)'*(phi - phi0);
    % phi: 2x1
    % n: 2xnq
    df = [0 2*K*(phi-phi0)'*n zeros(1,13)];
  end

  function [f,df] = final_cost_fun(T,x)
    K = 1;
    f = K*T;
    df = [K zeros(1,18)];
    
  end

end
