
N = 10;
T = .3;
duration = [T T];

% Trajectory optimization
g = 10;
z_nom = 1;
step_max = .7;
step_time = 0.3;
cop_max = .1; % set to 0 to get point foot model with no continuous inputs
fz_range = .3;
fx_range = .3;
inertia_ratio = .6^2/2;

model = TransformedFull2DModel(g, inertia_ratio, z_nom, step_max, step_time, fz_range, fx_range, cop_max);

plant = NStepCapturabilityPlant(model,false);
plant = plant.setInputLimits(-ones(3,1),ones(3,1));

[Y1,Y2] = meshgrid(linspace(-3,3,50),linspace(-2,2,50));
Y1 = Y1(:);
Y2 = Y2(:);
info_vec = zeros(size(Y1));
cost_vec = zeros(size(Y1));
for i=1:numel(Y1),
  
  x0 = [Y1(i);Y2(i);0;0;0;0];
  %%
  traj_opt = DircolTrajectoryOptimization(plant,N,duration);  
  traj_opt = traj_opt.addDecisionVariable(2);
  z_inds = traj_opt.num_vars + [-1;0];
  % traj_opt = traj_opt.addRunningCost(
  
  init_constraint = ConstantConstraint(x0);
  traj_opt = traj_opt.addStateConstraint(init_constraint,1);
  
  final_state = QuadraticConstraint(0,1,2*Q,zeros(6,1));
  traj_opt = traj_opt.addConstraint(final_state,[z_inds;traj_opt.x_inds(3:end,end)]);
  
  final_cost = QuadraticConstraint(-inf,inf,100*[1 -1;-1 1],zeros(2,1));
  traj_opt = traj_opt.addCost(final_cost,[z_inds(1);traj_opt.x_inds(1,end)]);
  traj_opt = traj_opt.addCost(final_cost,[z_inds(2);traj_opt.x_inds(2,end)]);
  
  tic
  [xtraj,utraj,z,F,info] = traj_opt.solveTraj(linspace(0,T,N));
  toc
  [info F(1)]
  
  info_vec(i) = info;
  cost_vec(i) = F(1);
  i
end
%%
figure(3)
Y1 = reshape(Y1,50,[]);
Y2 = reshape(Y2,50,[]);
cost_vec = reshape(cost_vec,50,[]);
hold off
[c,h]=contour(Y1,Y2,cost_vec,[step_max^2*50 step_max^2*50])
clabel(c,h)
hold on
[c,h]=contour(Y1,Y2,cost_vec,[1e-3 1e-3])
clabel(c,h)
%%
v = NStepCapturabilityVisualizer(plant);
v = v.setInputFrame(plant.getStateFrame);
v.playback_speed = T;
figure(25)
v.playback(xtraj);