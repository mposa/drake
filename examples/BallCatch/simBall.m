p = BallPlant;
h = 1e-4;
r = TimeSteppingRigidBodyManipulator(p,h);
[H,C,B] = p.manipulatorDynamics(zeros(9,1),zeros(9,1));
K = [     0     0   100     0     0     0     0     0     0     0     0    30     0     0     0     0     0     0
  0     0     0   100     0     0     0     0     0     0     0     0    30     0     0     0     0     0
  0     0     0     0     0     0     0   10     0     0     0     0     0     0     0     0    3     0
  0     0     0     0   100     0     0     0     0     0     0     0     0    30     0     0     0     0
  0     0     0     0     0   100     0     0     0     0     0     0     0     0    30     0     0     0
  0     0     0     0     0     0     0     0   10     0     0     0     0     0     0     0     0    3];

T = zeros(9);
T(1:2,1:2) = eye(2);
T(3:4,4:5) = eye(2);
T(5:6,7:8) = eye(2);
T(7:9,[3 6 9]) = eye(3);

g = 9.81;

v = p.constructVisualizer;

% xnom = [0;0;.05;0;-.05;0;0;0;0;zeros(9,1)];
xnom =  [0;.2;0;-.12;0;0;.12;0;0;zeros(9,1)];
xnom_frame = kron(eye(2),T')*xnom;
x0 = xnom_frame;
N = 10000;

x = zeros(18,N);
x(:,1) = x0;
for i=2:N,
  u = -K*(x(:,i-1) - xnom_frame) + [0;1.5*g;0;0;1.5*g;0];
  x(:,i) = r.update(0,x(:,i-1),u);
  if rem(i,50) == 0
    v.drawWrapper(i*h,x(:,i));
  end
end