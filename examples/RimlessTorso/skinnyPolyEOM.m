function [phi,phidot,psi,f_free,f_impact_1,f_impact_2,E,invH,invHJ_1,invHJ_2] = skinnyPolyEOM(degree)
%need phi, phidot, psi, H^-1(Bu-C), H^-1*J'*lambda, 
% lambda is ordered xp;xm;z;xp;xm;z;...

p = PlanarRigidBodyManipulator('WeightedSkinnyTorsoBalance.urdf');

xu = TaylorVar.init(zeros(11,1),degree);
q=xu(1:4);
qd=xu(5:8);
u=xu(9);
lx=xu(10:11);
lz = [1;1];
[H,C,B] = p.manipulatorDynamics(q,qd);
invH = inv(H);
kinsol = p.doKinematics(q);

[phi,n,D] = p.contactConstraints(kinsol);

f_free = invH*(B*u - C);
phidot = n*qd;
psi = D{1}*qd;
invHJ_1 = invH*[D{1}(1,:);n(1,:)]';
invHJ_2 = invH*[D{1}(2,:);n(2,:)]';

f_impact_1 = invHJ_1*[lx(1);lz(1)];
f_impact_2 = invHJ_2*[lx(2);lz(2)];

K = .5*qd'*H*qd;
g = 9.81;

U = 0;
for i=1:p.getNumBodies,
  m_i = p.getBody(i).mass;
  if m_i > 0
    com_i = p.forwardKin(kinsol,i,p.getBody(i).com);
    U = U + m_i*com_i(3)*g;
  end
end
U = U - double(U);
E = K + U;