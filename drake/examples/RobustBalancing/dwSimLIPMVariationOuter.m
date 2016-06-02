data_0 = load('V0_LIPMVariation2D');
data_1 = load('V1_LIPMVariation2D');
data = {data_1};
model = data_1.model;
N = 15;
[X,XD] = meshgrid(linspace(-.5,.5,N),linspace(-1.5,1.5,N));

xf = zeros(9,N^2);
for i=1:numel(X),
  x0 = [2;X(i);0;0;XD(i);0;0];
  traj = p.simulate([0 .3],[x0;0;0;0]);
  xf(:,i) = traj.eval(.3);
  i
end
xf = xf(1:6,:);
% plant = NStepCapturabilityPlant(model);
% v = NStepCapturabilityVisualizer(plant);
% v = v.setInputFrame(p.getOutputFrame);
% v.playback_speed = .5;

% figure(25)
% v.playback(traj);

%%
data_inner = load('lipmVariationWorkspace_7');
for i=1:numel(X)
  Vfi = subs(data_inner.V,data_inner.x,xf(:,i) + data_inner.x - data_inner.r);
  Vfo = subs(data_0.Vsol,[data_inner.t;data_inner.x],[0;xf(:,i) + data_inner.x - data_inner.r]);
  V0(i) = double(min(dmsubs(Vfi,data_inner.s,linspace(-1,1,100))));
  V0o(i) = double(max(dmsubs(Vfo,data_inner.s,linspace(-1,1,100))));
%   i
end