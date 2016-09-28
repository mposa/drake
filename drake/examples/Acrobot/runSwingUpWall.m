function [p,xtraj,utraj,ltraj,traj_opt,z] = runSwingUpWall(xtraj,utraj,ltraj)
options.terrain = RigidBodyFlatTerrain();
p = PlanarRigidBodyManipulator('AcrobotContact.urdf',options);
p = p.setGravity([-9.81;0;0]);
p = p.compile();

N = 30;
T = 15;
T0 = 5;

x0 = [pi/2;0;0;0];
xf = [3*pi/2;0;0;0];

N2 = floor(N/2);

if nargin < 2
  %Try to come up with a reasonable trajectory

  t_init = linspace(0,T0,N);
  traj_init.x = PPTrajectory(foh(t_init,linspacevec(x0,xf,N)));
  traj_init.u = PPTrajectory(foh(t_init,randn(1,N)));
  traj_init.l = PPTrajectory(foh(t_init,zeros(4,N)));
else
  t_init = xtraj.pp.breaks;
  traj_init.x = xtraj;
  traj_init.u = utraj;
  traj_init.l = ltraj;
end
T_span = [1 T];

to_options.nlcc_mode = 2;
to_options.lincc_mode = 1;
to_options.compl_slack = .0;
to_options.lincompl_slack = .0;
to_options.jlcompl_slack = .0;
to_options.lambda_mult = p.getMass*9.81*T0/N;
to_options.lambda_jl_mult = T0/N*100;

traj_opt = ContactImplicitTrajectoryOptimization(p,N,T_span,to_options);
traj_opt = traj_opt.addRunningCost(@running_cost_fun);
traj_opt = traj_opt.addStateConstraint(ConstantConstraint(x0),1);
traj_opt = traj_opt.addStateConstraint(ConstantConstraint(xf),N);
% traj_opt = traj_opt.setCheckGrad(true);
traj_opt = traj_opt.setSolverOptions('snopt','print','snopt.out');
traj_opt = traj_opt.setSolverOptions('snopt','MajorIterationsLimit',200);
traj_opt = traj_opt.setSolverOptions('snopt','MinorIterationsLimit',200000);
traj_opt = traj_opt.setSolverOptions('snopt','IterationsLimit',200000);
[xtraj,utraj,ltraj,ljltraj,z,F,info] = traj_opt.solveTraj(t_init,traj_init);

function [f,df] = running_cost_fun(h,x,u)
  f = h*u'*u;
  df = [u'*u zeros(1,4) 2*h*u'];
end

end
