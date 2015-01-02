function [p,xtraj,utraj,ltraj,ljltraj,z,F,info,traj_opt] = testRod(z0)
warning('off','Drake:RigidBodyManipulator:UnsupportedContactPoints');
warning('off','Drake:RigidBodyManipulator:WeldedLinkInd');
options.terrain = RigidBodyFlatTerrain();
options.floating = true;
options.ignore_self_collisions = true;
p = PlanarRigidBodyManipulator('Rod.urdf',options);
% trajopt = ContactImplicitTrajectoryOptimization(p,[],[],[],10,[1 1]);

%todo: add joint limits, periodicity constraint

N = 20;
T = 1;
T0 = 1;

x0 = [.5;1;1;-2;0;1];
xf = zeros(6,1);

N2 = floor(N/2);

  t_init = linspace(0,T0,N);
%   traj_init.x = PPTrajectory(foh(t_init,linspacevec(x0,xf,N)));
  traj_init.x = PPTrajectory(foh(t_init,randn(6,N)));
  traj_init.u = [];%PPTrajectory(foh(t_init,randn(0,N)));
  
  traj_init.lm = PPTrajectory(foh(t_init,rand(6,N)));
  traj_init.lp = PPTrajectory(foh(t_init,rand(6,N)));
  traj_init.L = PPTrajectory(foh(t_init,repmat([rand(3,1)],2,N)));
  
T_span = [T T];


x0_min = [x0(1:3);-inf(3,1)];
x0_max = [inf(2,1);x0(3);inf(3,1)];
% xf_min = xf;
% xf_max = xf;

x0_min = x0 - .0*ones(6,1);
x0_max = x0+.0*ones(6,1);
xf_min = -inf(6,1);
xf_max = inf(6,1);

to_options.nlcc_mode = 2;
to_options.lincc_mode = 1;
to_options.compl_slack = .000;
to_options.lincompl_slack = .0;
to_options.jlcompl_slack = .01;
to_options.lambda_mult = p.getMass*9.81;
% to_options.Lambda_mult = p.getMass;
to_options.lambda_jl_mult = T0/N;
% to_options.time_option = 1;

to_options.time_option = 2;



traj_opt = ColocatedContactImplicitTrajectoryOptimization(p,N,T_span,to_options);
traj_opt = traj_opt.addRunningCost(@running_cost_fun);
traj_opt = traj_opt.addStateConstraint(BoundingBoxConstraint(x0_min,x0_max),1);
traj_opt = traj_opt.addStateConstraint(BoundingBoxConstraint(xf_min,xf_max),N);


h_cons = LinearConstraint(.01*ones(N-1,1),2*T/(N-1)*ones(N-1,1),eye(N-1));
traj_opt = traj_opt.addConstraint(h_cons,traj_opt.h_inds);

% traj_opt = traj_opt.setCheckGrad(true);
snprint('snopt.out');
traj_opt = traj_opt.setSolverOptions('snopt','MajorIterationsLimit',300);
traj_opt = traj_opt.setSolverOptions('snopt','MinorIterationsLimit',20000);
traj_opt = traj_opt.setSolverOptions('snopt','IterationsLimit',100000);
traj_opt = traj_opt.setSolverOptions('snopt','Major Optimality Tolerance',1e-4);

if nargin > 0
  [xtraj,utraj,ltraj,ljltraj,z,F,info] = traj_opt.solveTrajFromZ0(z0);
else
  z0 = traj_opt.getInitialVars(t_init,traj_init);
  [f,df] = traj_opt.objectiveAndNonlinearConstraints(z0);
  [xtraj,utraj,ltraj,ljltraj,z,F,info] = traj_opt.solveTraj(t_init,traj_init);
end
[f,df] = traj_opt.objectiveAndNonlinearConstraints(z);

function [f,df] = running_cost_fun(h,x,u)
  f = h*x'*x + 1000*h^2;
  df = [2000*h+x'*x 2*h*x'];
% f = h^2;
% df = [2*h 0*x']; 
end

end
