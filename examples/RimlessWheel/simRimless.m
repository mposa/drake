N_trials = 1;

h = 1e-3;
T = 1;
N = floor(T/h);

x0_vec = zeros(6,N_trials);
xf_vec = zeros(6,N_trials);
x_vec = zeros(6,N,N_trials);

p = PlanarRigidBodyManipulator('RimlessWheel.urdf',struct('floating',true));
r = TimeSteppingRigidBodyManipulator(p,h);
v = p.constructVisualizer;
v.playback_speed = 1/2;
v.axis = [-1 1 -.1 1.5];
% A = [0;.4;1;1;1;1;.5;.5];

for trial_ind=1:N_trials,
%   x0 = .1*diag(A)*randn(8,1);
% x0 = [0;.075;.05;.08;zeros(4,1)];
% x0 = [0;.055;.05;.13;zeros(4,1)];
% x0 = [0;.055;.05;.4;zeros(4,1)]*2;
% x0 = [0;.09;.06;.29;zeros(4,1)]*3;
% x0 = [0;.00;.0;.87;zeros(4,1)];
x0 = [0;.075 + cos(pi/8);.44 - pi/8;zeros(3,1)];
% x0 = [0;.066;.23;0;zeros(4,1)]*2;
%   %check feasibility
%   s = sin(x0(3));
%   c = cos(x0(3));
%   z = x0(2);
%   phi(1) =  z - (8321567036706119*c)/9007199254740992 + (215431620425035*s)/562949953421312 + 1040195879588265/1125899906842624;
%   phi(2) =  z - (2^(1/2)*((8321567036706119*c)/9007199254740992 - (215431620425035*s)/562949953421312))/2 - (2^(1/2)*((215431620425035*c)/562949953421312 + (8321567036706119*s)/9007199254740992))/2 + 1040195879588265/1125899906842624;
%   
%   while sum(phi < 0) > 0
%     x0 = .1*diag(A)*randn(8,1);
%     
%     %check feasibility
%     s = sin(x0(3));
%     c = cos(x0(3));
%     z = x0(2);
%     phi(1) =  z - (8321567036706119*c)/9007199254740992 + (215431620425035*s)/562949953421312 + 1040195879588265/1125899906842624;
%     phi(2) =  z - (2^(1/2)*((8321567036706119*c)/9007199254740992 - (215431620425035*s)/562949953421312))/2 - (2^(1/2)*((215431620425035*c)/562949953421312 + (8321567036706119*s)/9007199254740992))/2 + 1040195879588265/1125899906842624;
%   end
%   x0
  v.draw(0,x0);
  
  t0 = 0.0;
  x = zeros(6,N);
  x(:,1) = x0;
  t = (0:N-1)*h;
  for i=2:N,
    [x(:,i)] = r.update(t(i-1),x(:,i-1),[]);
  end
  
%   v = p.constructVisualizer;
  xtraj = PPTrajectory(foh(t,x));
  xtraj = xtraj.setOutputFrame(r.getStateFrame);
  v.playback(xtraj);
  
  x0_vec(:,trial_ind) = x0;
  xf_vec(:,trial_ind) = x(:,N);
  xf = xf_vec(:,trial_ind)
  norm(xf_vec(2:end,trial_ind))
  
  x_vec(:,:,trial_ind) = x;
  
  X = zeros(6,N);
  X(1,:) = x(2,:) - cos(pi/8);
  X(2,:) = -sin(x(3,:) + pi/8);
  X(3,:) = cos(x(3,:) + pi/8);
  X(4,:) = x(4,:);
  X(5,:) = x(5,:);
  X(6,:) = -x(6,:);
  save torso_tmp_data xf_vec x0_vec x_vec
end