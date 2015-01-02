function [p,xtraj,utraj,ltraj,ljltraj,z,F,info,traj_opt] = trajRefinement(xtraj,utraj,ltraj,ljltraj,scale,t0,tf,is_iter)
warning('off','Drake:RigidBodyManipulator:UnsupportedContactPoints');
warning('off','Drake:RigidBodyManipulator:WeldedLinkInd');
options.terrain = RigidBodyFlatTerrain();
options.floating = true;
options.ignore_self_collisions = true;
p = PlanarRigidBodyManipulator('../KneedCompassGait.urdf',options);
% trajopt = ContactImplicitTrajectoryOptimization(p,[],[],[],10,[1 1]);

if nargin < 8
  is_iter = false;
end

%todo: add joint limits, periodicity constraint

N = 3;
T = tf-t0;

% x0 = [0;0;1;zeros(15,1)];
% xf = [0;0;1;zeros(15,1)];


N2 = floor(N/2);

t_init = linspace(0,T,N);
t_sample = linspace(t0,tf,N);

if is_iter
  traj_init.x = xtraj;
  traj_init.u = utraj;
  traj_init.l = ltraj;
  traj_init.ljl = ljltraj;
  
  x0 = xtraj.eval(0);
  xf = xtraj.eval(T);
else
  x0 = xtraj.eval(t0);
  xf = xtraj.eval(tf);
  traj_init.x = PPTrajectory(foh(t_init,xtraj.eval(t_sample)));
  traj_init.u = PPTrajectory(foh(t_init,utraj.eval(t_sample)));
  traj_init.l = PPTrajectory(foh(t_init,ltraj.eval(t_sample)));
  traj_init.ljl = [];%PPTrajectory(foh(t_init,ljltraj.eval(t_sample)));
end
T_span = [T T];


x0_min = x0;% - .001*ones(12,1);
x0_max = x0;% + .001*ones(12,1);
xf_min = [xf(1:6);-inf(6,1)];%xf - .001*ones(12,1);
xf_max = [xf(1:6);inf(6,1)];%xf + .001*ones(12,1);

to_options.compl_slack = scale*.1;
to_options.lincompl_slack = scale*.001;
to_options.jlcompl_slack = scale*.01;

to_options.nlcc_mode = 2;
to_options.lincc_mode = 1;
to_options.lambda_mult = p.getMass*9.81*T/N;
to_options.lambda_jl_mult = T/N;

traj_opt = ContactImplicitTrajectoryOptimization(p,N,T_span,to_options);
traj_opt = traj_opt.addRunningCost(@running_cost_fun);
% traj_opt = traj_opt.addRunningCost(@foot_height_fun);
traj_opt = traj_opt.addFinalCost(@(T,x) final_cost_fun(T,x,xf));
traj_opt = traj_opt.addStateConstraint(BoundingBoxConstraint(x0_min,x0_max),1);
traj_opt = traj_opt.addStateConstraint(BoundingBoxConstraint(xf_min,xf_max),N);

% min force
min_force = LinearConstraint(.25*ones(N,1),inf(N,1),kron(eye(N),[1 1]));
traj_opt = traj_opt.addConstraint(min_force,reshape(traj_opt.l_inds(1:traj_opt.nD+2:end,:),[],1));

% xr = randn(12,1);
% ur = randn(3,1);
% [f,df] = foot_height_fun(0,xr,ur);
% [f2,df2] = geval(@foot_height_fun,0,xr,ur,struct('grad_method','numerical'));

% l1 = traj_opt.l_inds(5:end,1:5);
% l2 = traj_opt.l_inds(1:4,end-4:end);
% traj_opt = traj_opt.addConstraint(ConstantConstraint(0*l1(:)),l1);
% traj_opt = traj_opt.addConstraint(ConstantConstraint(0*l2(:)),l2);
% traj_opt = traj_opt.setCheckGrad(true);
snprint('snopt.out');
traj_opt = traj_opt.setSolverOptions('snopt','MajorIterationsLimit',200);
traj_opt = traj_opt.setSolverOptions('snopt','MinorIterationsLimit',50000);
traj_opt = traj_opt.setSolverOptions('snopt','IterationsLimit',50000);
v = p.constructVisualizer; 

zerou = ConstantConstraint(zeros(3,1));
traj_opt = traj_opt.addInputConstraint(zerou,N);


% prevent knee from bending backwards
knee_constraint = LinearConstraint([-inf(3,1);0;-inf;0;-inf(6,1)],inf(12,1),eye(12));
traj_opt = traj_opt.addStateConstraint(knee_constraint,mat2cell(1:N,1,ones(N,1)));

% contact delta fun
% lz_inds = traj_opt.l_inds(1:4:end,:);
% contactDeltaCost = FunctionHandleConstraint(-inf,inf,size(lz_inds,2),@contact_delta_cost_fun,1);
% traj_opt = traj_opt.addCost(contactDeltaCost,lz_inds(1,:)');
% traj_opt = traj_opt.addCost(contactDeltaCost,lz_inds(2,:)');

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

  function [f,df] = final_cost_fun(T,x,xf)
    K = 1000;
    v_err = x - xf;
    f = K*(v_err'*v_err);
    df = [0 2*K*v_err'];
  end


end
