function [xf,traj,v]=simLIPMSwingAndHeightVariationSwitching(model,B,doPlot,x0)

if nargin < 3
  doPlot = true;
end
if nargin < 4
  x0 = [.1;0;0;0;0];
%   x0 = [
%    -0.2388
%    -0.0154
%     0.2043
%     0.5761
%    -0.0963
%     1.2199
% ]
end
% clear all
% load V1_LIPM3D
% V0 = load('V0_TransformedLIPM2D');
% data = {V0};
% model = V0.model;
% num = u;
% den = msspoly(1);
% data{1}.u_sol = usol;
% p = HybridCapturabilityPlant(model,data);

%%
t = msspoly('t',1);
x = msspoly('x',model.num_states);

plant = NStepCapturabilityPlant(model,true);
controller = SwitchingController(plant,t,x,B,true);
p = feedback(plant,controller);
p = p.setSimulinkParam('Solver', 'ode1');
p = p.setSimulinkParam('FixedStep', '0.01');    

% x0 = [1;-.1;.5];
% x0 = [1;-.29;.31];

% x0 = [1;-0;.1];
% x0 = [1;0;0];
% x0 = [1;-.2;0;.75;0];
% x0(2) = -x0(4)/sqrt(model.g);
T = .3;
% T = 2;

traj = p.simulate([0 T],[x0;0;0;0]);
if doPlot
plant = NStepCapturabilityPlant(model);
v = NStepCapturabilityVisualizer(plant);
v = v.setInputFrame(p.getOutputFrame);
v.playback_speed = 1;
figure(25)
v.playback(traj);
end
xf = traj.eval(T);
xf = xf(1:end-3);

%% 6

% t_sim = traj.pp.breaks;
% x_sim = traj.eval(t_sim);
% V0 = msubs(V0.Vsol,[t;x],[t_sim;x_sim(1:model.num_states,:)]);
% figure(5)
% subplot(3,1,1)
% plot(t_sim,x_sim(1:model.num_states,:))
% subplot(3,1,2)
% plot(t_sim,V0)
% grid on
% subplot(3,1,3)
% plot(t_sim,sqrt(sum(x_sim(1:2,:).*x_sim(1:2,:),1)))
% grid on

% min(sqrt(sum(x_sim(1:2,:).*x_sim(1:2,:),1)))
end