th_range = .1;
% x0 = [0;.049;0;0;0;0];
T = .1;
x0 = [0;.549;pi/2;0;0;0];
sample_range = [.02;.02;.3;.02;.02;.02];
T = .5;
rbmoptions.floating = true;
rbmoptions.ignore_self_collisions = true;
rbmoptions.terrain = RigidBodyFlatTerrain();
rbmoptions.twoD = true;
r = TimeSteppingRigidBodyManipulator('Rod.urdf',1e-3,rbmoptions);

%%
N = 500;
N = 1;
xf = zeros(6,N);
for i=1:N,
  x_sample = (rand(6,1) - .5)*2.*sample_range + x0;
  phi = p.contactConstraints(x_sample(1:3));
  while any(phi) < 0
    x_sample = (rand(6,1) - .5)*2.*sample_range + x0;
    phi = p.contactConstraints(x_sample(1:3));
  end
  traj = r.simulate([0 T],x_sample);
  xf(:,i) = traj.eval(T);
end

% v = r.constructVisualizer;

%%
  x_sample = (rand(6,1) - .5)*2.*sample_range + x0;
  traj = r.simulate([0 T],x_sample);
  v.playback(traj)
  traj.eval(.1)*1