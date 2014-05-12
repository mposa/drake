N_trials = 1;

pitchd_val = linspace(-5,5,100);
thetad_val = linspace(-5,5,100);
[PITCHD,THETAD] = meshgrid(pitchd_val,thetad_val);
N_trials = length(PITCHD(:));

h = 5e-3;
T = 5;
N = floor(T/h);

x0_vec = zeros(8,N_trials);
xf_vec = zeros(8,N_trials);
x_vec = zeros(8,N,N_trials);
load(datapath('zscale_skinny_cubic_controller_taylor_iter_68'));
z_scale = .1;
Vsol = subs(Vsol,[z;zd],[z;zd]/z_scale);
q = msspoly('q',4);
qd = msspoly('qd',4);
[~,u_fun,u_coef,u_perm] = dmsubs(controllersol,[q;qd],zeros(8,1));

p = PlanarRigidBodyManipulator('SkinnyTorsoBalance.urdf');
% p = PlanarRigidBodyManipulator('TorsoBalance.urdf');
% p = PlanarRigidBodyManipulator('ReactionRimless.urdf');
r = TimeSteppingRigidBodyManipulator(p,h);
v = p.constructVisualizer;
v.playback_speed = 1;
v.axis = [-1 1 -.1 1.5];
% A = [0;.4;1;1;1;1;.5;.5];



for trial_ind=1:N_trials,
%   x0 = .1*diag(A)*randn(8,1);
% x0 = [0;.075;.05;.08;zeros(4,1)];
% x0 = [0;.055;.05;.13;zeros(4,1)]*1.8;
% x0 = [0;.055;.14;.24;zeros(4,1)]/2;
% x0 = [0;.015;.08;.22;zeros(4,1)];
% x0 = [0;.017;.1;.12;zeros(4,1)]*1.2;
% x0 = [0;.075;-.14;.18;zeros(4,1)];
% x0 = [0;.02;.1;.06;zeros(4,1)];
% x0 = [0;0;0;0;0;0;2;-2];
% x0 = [0;.05;-.09;-.16;zeros(4,1)]*1.2;
% x0 = [0;        0.0027
%     0.0089
%     0.0424
%     0.2719
%    -0.3734
%    -0.6247
%     0.1427];
% x0 = [0;.055;.05;.4;zeros(4,1)];
% x0 = [0;.09;.06;.29;zeros(4,1)]r;
% x0 = [0;.00;.0;.87;zeros(4,1)];
% x0 = [0;.066;.23;0;zeros(4,1)];

x0 = zeros(8,1);
x0(7) = PITCHD(trial_ind);
x0(8) = THETAD(trial_ind);
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
%   v.draw(0,x0);
  
  t0 = 0.0;
  x = zeros(8,N);
  x(:,1) = x0;
  u = zeros(1,N);
  t = (0:N-1)*h;
  doBreak = false;
  for i=2:N,
%     u = 0;
    K = [0 0 0 -10 0 0 0 -1];
%     u = -[25 1]*[sin(x(4,i-1)); x(8,i-1)];
    
%     u = [  (0.091748)*x(5,i-1)+(0.0035586)*x(2,i-1)+(0.0011083)*x(6,i-1)+(2.1977)*x(3,i-1)+(-0.64263)*x(7,i-1)+(-11.686)*x(4,i-1)+(-2.4209)*x(8,i-1)  ];
%     u = subs(controllersol,[q;qd],x(:,i-1));

      u(i) = u_fun(u_coef,x(u_perm,i-1));
%     u(i) = K*x(:,i-1);

    [x(:,i)] = r.update(t(i-1),x(:,i-1),u(i));
%     [~,zz(:,i)] = r.solveLCP(t(i-1),x(:,i-1),u(i));
    
    if norm(x(:,i)) > 1e3
      doBreak = true;
      break
    end
  end
  
%   v = p.constructVisualizer;
%   xtraj = PPTrajectory(foh(t,x));
%   xtraj = xtraj.setOutputFrame(r.getStateFrame);
%   v.playback(xtraj);

  x0_vec(:,trial_ind) = x0;
  
  if doBreak
    xf_vec(:,trial_ind) = x(:,i);
  else
    xf_vec(:,trial_ind) = x(:,N);
  end
  
  xf = xf_vec(:,trial_ind);
  norm(xf_vec(2:end,trial_ind))
  
  x_vec(:,:,trial_ind) = x;
  
  save torso_tmp_data xf_vec x0_vec x_vec
  trial_ind
end