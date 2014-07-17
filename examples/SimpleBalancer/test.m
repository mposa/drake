clear all
rbmoptions.twoD = true;
rbmoptions.ignore_self_collisions = true;
p = PlanarRigidBodyManipulator('SimpleBalancer.urdf',rbmoptions);
q_rel = msspoly('qr',3);
v_rel = msspoly('vr',3);
s_rel = msspoly('sr',3);
c_rel = msspoly('cr',3);

q_trig=TrigPoly(q_rel,s_rel,c_rel);

[H,C,B] = p.manipulatorDynamics(q_trig,v_rel);
H = getmsspoly(H);
C = getmsspoly(C);

q = msspoly('q',3);
v = msspoly('v',3);
s = msspoly('s',3);
c = msspoly('c',3);
prog = spotsosprog;

H_bkp = prog.trigExprReduction(H,s_rel,c_rel);

%% substitute relative angles
H = subs(H,c_rel(1),c(1));
H = subs(H,s_rel(1),s(1));
H = subs(H,v_rel(1),v(1));

H = subs(H,c_rel(2),c(1)*c(2) + s(1)*s(2));
H = subs(H,s_rel(2),s(1)*c(2) - s(2)*c(1));
% H = subs(H,v_rel(2),v(2) - v(1));

H = subs(H,c_rel(3),c(2)*c(3) + s(2)*s(3));
H = subs(H,s_rel(3),s(2)*c(3) - s(3)*c(2));
% H = subs(H,v_rel(3),v(3) - v(2));

H = prog.trigExprReduction(H,s,c)