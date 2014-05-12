q = msspoly('q',3);
qd = msspoly('qd',3);
x = q(1);
z = q(2);
theta = q(3);
xd = qd(1);
zd = qd(2);
thetad = qd(3);
c = msspoly('c',1);
s = msspoly('s',1);
lx = msspoly('lx',1);
lz = msspoly('lz',1);

l = 1;
m = 1;
I = m*l^2/50;
mu = 1;
th0 = 1;
thd0 = 0;
g = 10;

1/m - l^2/4/I*cos(th0)*(mu*sin(th0) - cos(th0))

H = diag([m m I]);
C = [0;m*g;0];

phi = z - l/2*s;
J = [0;1;-l/2*c]';
Jf = [1;0;l/2*s]';

psi = Jf*qd;
phid = J*qd;

Jdot = [0 0 .5*s*thetad];

qdd = inv(H)*(-C + J'*lz + Jf'*lx);

phidd = Jdot*qd + J*qdd;

subs(phidd,[s;c;thetad],[sin(th0);cos(th0);0])