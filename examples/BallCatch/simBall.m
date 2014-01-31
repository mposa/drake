p = BallPlant;
h = 1e-4;
r = TimeSteppingRigidBodyManipulator(p,h);
[H,C,B] = p.manipulatorDynamics(zeros(9,1),zeros(9,1));
K = [     0     0   500     0     0     0     0     0     0     0     0    30     0     0     0     0     0     0
     0     0     0   500     0     0     0     0     0     0     0     0    30     0     0     0     0     0
     0     0     0     0     0     0     0   10     0     0     0     0     0     0     0     0    3     0
     0     0     0     0   500     0     0     0     0     0     0     0     0    30     0     0     0     0
     0     0     0     0     0   500     0     0     0     0     0     0     0     0    30     0     0     0
     0     0     0     0     0     0     0     0   10     0     0     0     0     0     0     0     0    3];
v = p.constructVisualizer;

xnom = [0;0;.05;0;-.05;0;0;0;0;zeros(9,1)];
x0 = xnom;
x0(2) = .15;
x0(11) = 0;
x = x0;
for i=1:5000,
  u = -K*(x - xnom);
  x = r.update(0,x,u);
  if rem(i,10) == 0
    v.drawWrapper(i*h,x);
  end
end