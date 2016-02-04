data_0 = load('V0_VariableHeightPointMass2D');
% data_1 = load('V1_VariableHeightPointMass2D');
data = {data_0;data_1};
data = {data_0};
% data = {data_1};
model = data_0.model;
p = HybridCapturabilityPlant(model,data);
%%
% x0 = [1;-.17;0;.57;0];
x0 = [1;-.1;0;.3;0];
% x0 = [2;-.26;0;1.8;0];
% x0 = [2;.0;0;1;0];
% x0 = [2;0;0;0;0];
% x0 = [3;-.9;0;1.4;0];
% x0 = [3;0;0;1;0];
% T = .3*(x0(1)-1)+1;
T = 1;
traj = p.simulate([0 T],[x0;0;0;0]);
plant = NStepCapturabilityPlant(model);
v = NStepCapturabilityVisualizer(plant);
v = v.setInputFrame(p.getOutputFrame);
v.playback_speed = .5;

v.playback(traj);


%%
% sample sim
%%
t = msspoly('t',1)
x = msspoly('x',model.num_states);
t_sim = traj.pp.breaks;
x_sim = traj.eval(t_sim);
V0 = msubs(data_0.Vsol,[t;x],[t_sim;x_sim(1:model.num_states,:)]);