function [p,xtraj,utraj,ltraj,ljltraj,z,F,info,traj_opt] = testColocatedTrajOpt;
warning('off','Drake:RigidBodyManipulator:UnsupportedContactPoints');
warning('off','Drake:RigidBodyManipulator:WeldedLinkInd');
% options.terrain = RigidBodyFlatTerrain();
options.floating = false;
% options.ignore_self_collisions = true;
p = PlanarRigidBodyManipulator('../Acrobot.urdf',options);
% trajopt = ContactImplicitTrajectoryOptimization(p,[],[],[],10,[1 1]);

%todo: add joint limits, periodicity constraint

N = 20;
T = 3;
T0 = 3;

% x0 = [0;0;1;zeros(15,1)];
% xf = [0;0;1;zeros(15,1)];
x0 = [0;0;0;0];
xf = [pi;0;0;0];

t_init = linspace(0,T0,N);
traj_init.x = PPTrajectory(foh(t_init,linspacevec(x0,xf,N)));
traj_init.u = PPTrajectory(foh(t_init,randn(1,N)));

T_span = [1 T];


x0_min = x0;
x0_max = x0;
xf_min = xf;
xf_max = xf;

to_options.nlcc_mode = 2;
to_options.lincc_mode = 1;
to_options.compl_slack = .01;
to_options.lincompl_slack = .1;
to_options.jlcompl_slack = .01;
to_options.lambda_mult = p.getMass*9.81;
% to_options.Lambda_mult = p.getMass;
to_options.lambda_jl_mult = T0/N;
to_options.time_option = 1;


traj_opt = ColocatedContactImplicitTrajectoryOptimization(p,N,T_span,to_options);
traj_opt = traj_opt.addRunningCost(@running_cost_fun);
traj_opt = traj_opt.addStateConstraint(BoundingBoxConstraint(x0_min,x0_max),1);
traj_opt = traj_opt.addStateConstraint(BoundingBoxConstraint(xf_min,xf_max),N);
% traj_opt = traj_opt.setCheckGrad(true);
snprint('snopt.out');
traj_opt = traj_opt.setSolverOptions('snopt','MajorIterationsLimit',100);
traj_opt = traj_opt.setSolverOptions('snopt','MinorIterationsLimit',20000);
traj_opt = traj_opt.setSolverOptions('snopt','IterationsLimit',100000);
z0 = traj_opt.getInitialVars(t_init,traj_init);
% [f,df] = traj_opt.objectiveAndNonlinearConstraints(z0);
[xtraj,utraj,ltraj,ljltraj,z,F,info] = traj_opt.solveTraj(t_init,traj_init);

[f,df] = traj_opt.objectiveAndNonlinearConstraints(z);

function [f,df] = running_cost_fun(h,x,u)
  f = h*u'*u;
  df = [u'*u zeros(1,4) 2*h*u'];
end

end
