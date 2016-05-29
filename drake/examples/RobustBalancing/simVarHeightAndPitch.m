data_0 = load('V0_VariableHeightandPitch2D');
% data_1 = load('V1_VariableHeightandPitch2D');
% data = {data_0;data_1};
data = {data_0};
model = data_0.model;
p = HybridCapturabilityPlant(model,data);
%%

% x0 = [2;.0;0;0;1;0;0];
% x0 = [2;.7;0;0;0;0;0];
x0 = [1;.05;0;0;.0;0;0];

traj = p.simulate([0 .3*(x0(1)-1)+1],[x0;0;0;0]);
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
V0 = msubs(data_0.Vsol,[t;x],[t_sim;x_sim(1:model.num_states,:)]);
% V1 = msubs(data_1.Vsol,[t;x],[t_sim;x_sim(1:model.num_states,:)]);

figure(4)
subplot(3,1,1)
plot(t_sim,x_sim)
subplot(3,1,2)
plot(t_sim,V0);
subplot(3,1,3)
% plot(t_sim,V1);