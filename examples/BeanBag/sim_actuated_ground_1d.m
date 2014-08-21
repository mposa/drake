p = PlanarRigidBodyManipulator('actuated_ground.urdf');
p = p.addLinksToCollisionFilterGroup('ground','ignore',0);
r = TimeSteppingRigidBodyManipulator(p,1e-3);
controller = AffineSystem([],[],[],[],[],[],[],[0 -10 0 -1], p.getMass*norm(p.getGravity));
controller = controller.setOutputFrame(p.getInputFrame);
controller = controller.setInputFrame(p.getStateFrame);
sys_cl = feedback(r,controller);
traj = sys_cl.simulate([0 5],[1;.5;-5;-2]);

v = r.constructVisualizer;
v.playback(traj,struct('slider',true))