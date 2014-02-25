function [H,C,B,phi,phidot,psi,J,J_f,K,S,U] = reactionEOM_mss(q,v,s_vec,c_vec)
warning('off','Drake:RigidBody:SimplifiedCollisionGeometry');
warning('off', 'Drake:RigidBodyManipulator:UnsupportedContactPoints');
p = PlanarRigidBodyManipulator('ReactionRimless.urdf');

x = q(1);
z =  q(2);

xd = v(1);
zd = v(2);
pitchd = v(3);
thetad = v(4);

s = s_vec(3);
c = c_vec(3);

s_th = s_vec(4);
c_th = c_vec(4);

q_trig=TrigPoly(q,s_vec,c_vec);

[H,C,B] = p.manipulatorDynamics(q_trig,v);
H = getmsspoly(H);
C = getmsspoly(C);

for i=1:4,
  for j=1:4,
    H(i,j) = trigExprReduction(H(i,j),[s;s_th],[c;c_th]);
  end
  C(i) = trigExprReduction(C(i),[s;s_th],[c;c_th]);
end


phi = [z - (8321567036706119*c)/9007199254740992 + (215431620425035*s)/562949953421312 + 1040195879588265/1125899906842624;
       z - (2^(1/2)*((8321567036706119*c)/9007199254740992 - (215431620425035*s)/562949953421312))/2 - (2^(1/2)*((215431620425035*c)/562949953421312 + (8321567036706119*s)/9007199254740992))/2 + 1040195879588265/1125899906842624];

phidot =  [zd - (8321567036706119*-s*pitchd)/9007199254740992 + (215431620425035*c*pitchd)/562949953421312;
           zd - (2^(1/2)*((8321567036706119*-s*pitchd)/9007199254740992 - (215431620425035*c*pitchd)/562949953421312))/2 - (2^(1/2)*((215431620425035*-s*pitchd)/562949953421312 + (8321567036706119*c*pitchd)/9007199254740992))/2];

psi = [xd + (8321567036706119*pitchd*c)/9007199254740992 - (215431620425035*pitchd*s)/562949953421312;
       xd + (2^(1/2)*((8321567036706119*pitchd*c)/9007199254740992 - (215431620425035*pitchd*s)/562949953421312))/2 + (2^(1/2)*((215431620425035*pitchd*c)/562949953421312 + (8321567036706119*pitchd*s)/9007199254740992))/2];

phi_perp = [x + (8321567036706119*s)/9007199254740992 + (215431620425035*c)/562949953421312;
            x + (8321567036706119*s)/9007199254740992 - (215431620425035*c)/562949953421312];

% J = [0 1 (8321567036706119*s)/9007199254740992 -(215431620425035*c)/562949953421312 0;
%   0 1 -(8321567036706119*s)/9007199254740992 -(215431620425035*c)/562949953421312 0];
% 
% J_f = [1 0 ((8321567036706119*c)/9007199254740992 + (215431620425035*s)/562949953421312) 0;
%   1 0 ((8321567036706119*c)/9007199254740992 - (215431620425035*s)/562949953421312) 0];

J = [0 1 (8321567036706119*s)/9007199254740992 + (215431620425035*c)/562949953421312 0;
  0 1 (8321567036706119*s)/9007199254740992 - (215431620425035*c)/562949953421312 0];

J_f = [1 0 ((8321567036706119*c)/9007199254740992 - (215431620425035*s)/562949953421312) 0;
  1 0 ((8321567036706119*c)/9007199254740992 + (215431620425035*s)/562949953421312) 0];


% J_f = J_f([2;1],:);


%potential energy
kinsol = p.doKinematics(q_trig);
% doKinematics(p.model,q_trig);
g = -p.gravity;
U = 0;
for i=1:length(p.body),
  m = p.body(i).mass;
  if m > 0
    com = p.forwardKin(kinsol,i,p.body(i).com);
    U = U + com'*g*m;
  end
end
U = getmsspoly(U);
U = U - subs(U,[q;s_vec;c_vec],[zeros(size(q));zeros(size(s_vec));ones(size(c_vec))]);


% Generate a linear model
[A,B_lin] = p.linearize(0,zeros(8,1),0);
A_sub = A([4;8],[4;8]);
B_sub = B_lin([4;8]);
Q = eye(2);
R = 1;
[K_sub,S_sub] = lqr(A_sub,B_sub,Q,R);
K = K_sub;
S = S_sub;
% K = [zeros(1,3) K_sub(1) zeros(1,3) K_sub(2)];
% S = zeros(8);
% S([4;8],[4;8]) = S_sub;
end




