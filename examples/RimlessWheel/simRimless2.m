% clear all
h = 1e-2;
options.floating = true;
options.terrain = RigidBodyFlatTerrain;
options.twoD = true;
options.ignore_self_collisions = true;
p = TimeSteppingRigidBodyManipulator('RimlessWheel.urdf',h,options);

% x0 = [0;.2+.076 + cos(pi/8);.11 - pi/8;zeros(3,1)];
% x0 = [0;.15556+cos(pi/8);-.50758 - pi/8;0;0;0];
% x0 = [0;.135+cos(pi/8);.053 - pi/8;0;0;2];

z0 = .13131;
theta0 = -.674;

qd0 = [1.2533;.3639;1.3079];

x0 = [z0*sin(.1);z0*cos(.1) + cos(pi/8);theta0 + pi/8 + .1;qd0];

% x0 = [x0(1);0;x0(2);0;x0(3);0;x0(4);0;x0(5);0;x0(6);0];



v = p.constructVisualizer;
v.drawWrapper(0,x0);
v.axis = [-3 10 -1 2];

traj = p.simulate([0 10],x0);
v.playback(traj)


%%
Vsol = sol.eval(V);
Vdot_freesol = sol.eval(Vdot_free);
ground_angle = .1;
tt = traj.tt;
ttt = tt/tt(end);
xx = traj.xx;
xx_th = xx(3,:) - pi/8 - ground_angle;
xx_th = mod(xx_th + pi/8,pi/4) - pi/8;
xxx = [xx(2,:)*cos(ground_angle) + xx(1,:)*sin(ground_angle) - cos(pi/8); sin(xx_th);cos(xx_th); ... 
  xx(4,:)*cos(ground_angle) - xx(5,:)*sin(ground_angle);xx(5,:)*cos(ground_angle) + xx(4,:)*sin(ground_angle);xx(6,:)];
Vval = msubs(Vsol,[t;q(2:end);qd],[ttt;xxx]);

h_Bval = msubs(h_B,[t;q(2:end);qd],[ttt;xxx]);
h_BTval = msubs(h_BT,[t;q(2:end);qd],[ttt;xxx]);

Vdotval = msubs(Vdot_freesol,[t;q(2:end);qd],[ttt;xxx]);

Vdot_impact_1_sol = sol.eval(Vdot_impact_1);
Vdot_impact_2_sol = sol.eval(Vdot_impact_2);

psi1val = msubs(psi{1},[t;q(2:end);qd],[ttt;xxx]);
psi2val = msubs(psi{2},[t;q(2:end);qd],[ttt;xxx]);

Vdoti1val = msubs(Vdot_impact_1_sol,[t;q(2:end);qd;lx(1)],[ttt;xxx;-sign(psi1val)]);
Vdoti2val = msubs(Vdot_impact_2_sol,[t;q(2:end);qd;lx(2)],[ttt;xxx;-sign(psi2val)]);