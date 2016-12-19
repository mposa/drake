
N = 10;
T = 5;
duration = [T T];

% Trajectory optimization
g = 1;
z_nom = 1;
R_diag = [2, 1, 1, 2, 2, 2];
fx_range = 1;
fz_range = .5;
inertia_ratio = .6^2/2; 
step_max = .7;
step_time = .3;
foot_size = .4;

model = ScaledFull2D(g, z_nom, step_max, step_time,  inertia_ratio, fx_range, fz_range, foot_size);

A = diag(1./(R_diag.^2));

plant = NStepCapturabilityPlant(model,false);
plant = plant.setInputLimits(-ones(3,1),ones(3,1));

[X,XD] = meshgrid(linspace(-2,2,50),linspace(-2,2,50));
X = X(:);
XD = XD(:);
info_vec = zeros(size(X));
cost_vec = zeros(size(X));
for i=1:length(X),
  
  x0 = [X(i);0;0;XD(i);0;0];
  
  if x0'*A*x0 < 1
  
  traj_opt = DircolTrajectoryOptimization(plant,N,duration);
  
  % traj_opt = traj_opt.addRunningCost(
  
  init_constraint = ConstantConstraint(x0);
  traj_opt = traj_opt.addStateConstraint(init_constraint,1);
  
  final_cost = QuadraticConstraint(-inf,inf,100*eye(6),zeros(6,1));
  traj_opt = traj_opt.addCost(final_cost,traj_opt.x_inds(:,end));
  
  state_ball = QuadraticConstraint(0,1,A,zeros(6,1));
  traj_opt = traj_opt.addStateConstraint(state_ball,2:N);
  
  tic
  [xtraj,utraj,z,F,info] = traj_opt.solveTraj(linspace(0,T,N));
  toc
  [info F(1)]
  
  info_vec(i) = info;
  cost_vec(i) = F(1);
  else
    info_vec(i) = inf;
    cost_vec(i) = inf;
  end
  i
end
%%
v = NStepCapturabilityVisualizer(plant);
v = v.setInputFrame(plant.getStateFrame);
v.playback_speed = T;
figure(25)
v.playback(xtraj);