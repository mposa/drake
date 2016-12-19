clear all
% load V1_LIPM3D
V0 = load('V0_LIPM3D');
data = {V0};
model = V0.model;

p = HybridCapturabilityPlant(model,data);

%%


x0 = [1;-.0846;0;.4;0];
% x0 = [1;-.2;0;.75;0];
% x0(2) = -x0(4)/sqrt(model.g);
T = 1;

traj = p.simulate([0 T],[x0;0;0;0]);
plant = NStepCapturabilityPlant(model);
v = NStepCapturabilityVisualizer(plant);
v = v.setInputFrame(p.getOutputFrame);
v.playback_speed = .5;
figure(25)
v.playback(traj);

%%
t = msspoly('t',1)
x = msspoly('x',model.num_states);
t_sim = traj.pp.breaks;
x_sim = traj.eval(t_sim);
V0 = msubs(V0.Vsol,[t;x],[t_sim;x_sim(1:model.num_states,:)]);
figure(5)
subplot(2,1,1)
plot(t_sim,x_sim)
subplot(2,1,2)
plot(t_sim,V0)
grid on