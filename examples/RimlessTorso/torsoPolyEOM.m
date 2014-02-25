function [phi_poly,phidot_poly,psi_poly,f_free_poly,f_impact_poly,E_poly,q,qd,lambda_x,lambda_z,u,invHJ,invH] = torsoPolyEOM(degree)
%need phi, phidot, psi, H^-1(Bu-C), H^-1*J'*lambda, 
% lambda is ordered xp;xm;z;xp;xm;z;...

p = PlanarRigidBodyManipulator('TorsoBalance.urdf');

xu = TaylorVar.init(zeros(9,1),degree);


q=msspoly('q',p.num_q);
qd=msspoly('v',p.num_q);
u=msspoly('u',1);
lambda_x=msspoly('lx',2);
lambda_z=msspoly('lz',2);
lambda = [lambda_x(1);lambda_z(1);lambda_x(2);lambda_z(2)];
[H,C,B] = p.manipulatorDynamics(xu(1:4),xu(5:8))
% B = p.model.B;

% H=p.getH(xu(1:4));
invH = inv(H);
% C=p.getC(xu(1:4),xu(5:8));

kinsol = p.doKinematics(xu(1:4));
% doKinematicsAndVelocities(p.model,xu(1:4),xu(5:8));
[phi, psi, dPhi, dPsi, J, dJ,psi_full]  = p.contactPositionsAndVelocities(xu(1:4),xu(5:8));
% J_exp = [1 0 0 0;0 1 0 0; 0 -1 0 0; 0 0 1 0; 0 0 0 1; 0 0 0 -1] *J;
% J_exp = [1 0 0 0;-1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 -1 0; 0 0 0 1] *J;
J_exp = J;
plant_dyn = invH*(B*xu(9) - C);
f_free_poly = getmsspoly(plant_dyn,[q;qd;u]);
% keyboard
phi_poly = getmsspoly(phi,[q;qd;u]);
phi_poly = phi_poly([2;4]);

psi_poly = getmsspoly(psi,[q;qd;u]);
phidot_poly = getmsspoly(psi_full([2 4]),[q;qd;u]);

% qt2 = TaylorVar.init(zeros(4,1),degree-1);
% H2 = p.getH(qt2);
% doKinematicsAndVelocities(p.model,qt2,zeros(4,1));
% [~,~,~,~,J2,~] = p.contactPositionsAndVelocities(qt2,zeros(4,1));
% J2_exp = [1 0 0 0;0 1 0 0; 0 -1 0 0; 0 0 1 0; 0 0 0 1; 0 0 0 -1] *J2;
% f_impact_poly = getmsspoly(inv(H2)*J2_exp',q)*lambda;
invHJ = getmsspoly(invH*J_exp',[q;qd;u]);
invH = getmsspoly(invH,[q;qd;u]);
f_impact_poly = invHJ*lambda;

K = .5*xu(5:8)'*H*xu(5:8);
g = 9.81;

U = 0;
for i=1:length(p.model.body),
  T=p.model.body(i).T;
  I = p.model.body(i).I;
  if I(3,3) ~= 0
    tmp = T*[I(1,3)/I(3,3);-I(1,2)/I(3,3);1];
    U = U + tmp(2)*g*I(3,3);
  end
end
U = U - double(U);
E_poly = getmsspoly(K+U,[q;qd;u]);