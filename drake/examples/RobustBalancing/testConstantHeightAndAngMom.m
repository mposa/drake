options.floating = true;
options.terrain = RigidBodyFlatTerrain;
p = PlanarRigidBodyManipulator('PlanarWalker.urdf',options);

%%
load controller
r = TimeSteppingRigidBodyManipulator(p,1e-3);
c = BalancingController(p,coeffs,pows);
sys = r.feedback(c);
q0 = [0;1.35;.137;-1.2;1;-.5;.8];
phi = p.contactConstraints(q0)
q0(2) = q0(2) - phi(2);
x0 = [q0;zeros(7,1)];
traj = sys.simulate([0 .2],[q0;0*q0]);