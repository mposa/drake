clear all
load lipmSwingAndHeightVariation_switching_1
% B2 = [Bu;Bu(1)];
%%
N = 20;
[X,XD] = meshgrid(linspace(-.5,.5,N),linspace(-1.5,1.5,N)); 

for i=1:numel(X),
  x0 = [X(i);0;XD(i);0;0];
  [xf,traj]=simLIPMSwingAndHeightVariationSwitching(model_swing,B2,false,x0);
  ts = traj.pp.breaks;
  xs = traj.eval(ts);
  constraint_i = msubs(reset_constraint,x,xs(1:5,:))';
  a = double(diff(constraint_i,s));
  b = dmsubs(constraint_i,s,0);
  s_i = -b./a;
  Vp_i = dmsubs(V,x,msubs(r,[x;s],[xs(1:5,:);s_i']));
  Vp{i} = Vp_i;
  traj_list{i} = traj;
  i
end
%%
for i=1:numel(X),
  traj = traj_list{i};
  ts = traj.pp.breaks;
  xs = traj.eval(ts);
  V2_i = dmsubs(V2,[t;x],[ts;xs(1:5,:)]);
  V2_list(:,i) = V2_i;
  V2_end(i) = V2_i(end);
  V2_0(i) = V2_i(1);
  xs_end(:,i) = xs(1:5,end);
end

%%
for i=1:numel(X),
  [Vpmin(i),j] = min(Vp{i});
  ts_safe{i} = find(Vp{i} < 1);
  
  if ~isempty(ts_safe{i})
    ts_min(i) = ts(min(ts_safe{i}));
    ts_max(i) = ts(max(ts_safe{i}));
    ts_opt(i) = ts(j);
  else
    ts_min(i) = inf;
    ts_max(i) = inf;
    ts_opt(i) = inf;
  end
end

%%
figure(4)
colormap default
[cl,h] = contour(X,XD,reshape(ts_opt,N,N));
set(h,'Fill','On');
colorbar

set(gca,'LooseInset',get(gca,'TightInset'))
xlabel('x_c_m','FontSize',24)
ylabel('xdot_c_m','FontSize',24)
title('Stepping Times (Sampled)','FontSize',24)
hold off

%%
figure(3)
  xf_plot = .23333;
  t_plot = step_time;
  hold off
  contourSpotless(V,x(1),x(3),[-1 1],[-3 3],[t;x([2;4;5])],zeros(model_swing.num_states-1,1),1,{'g'});
  hold on
  contourSpotless(V2,x(1),x(3),[-1 1],[-3 3],[t;x([2;4;5])],[zeros(model_swing.num_states-2,1);xf_plot],dmsubs(rho2,t,0),{'k'});
  contourSpotless(V2,x(1),x(3),[-1 1],[-3 3],[t;x([2;4;5])],[step_time;zeros(model_swing.num_states-3,1);xf_plot],dmsubs(rho2,t,step_time),{'r'});