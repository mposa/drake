megaclear
degree = 6;

R = 1; % don't change this
R_diag = [.2 2 2 4];
RT_diag = R_diag;
% R_diag = [1 1 1 1];

RT = .05;
theta_bound = .75;
thetaT_bound = .01;
A = 1/R*diag(1./(R_diag.^2));
AT = 1/RT*diag(1./(RT_diag.^2));
% AT = 1/RT*diag([1;1;1;1]);

g = 9.81;

prog = spotsosprog();

%% Add indeterminate variables
q = msspoly('q',4);
qd = msspoly('qd',3);
lx = msspoly('lx',2);
lzsq = [1;1];

x = q(1);
z = q(2);
s = q(3);
c = q(4);

xd = qd(1);
zd = qd(2);
pitchd = qd(3);

x_vars = [q(2:end);qd];

prog = prog.withIndeterminate(q(2:end));
prog = prog.withIndeterminate(qd);
prog = prog.withIndeterminate(lx);

%% Functions V, W
[prog,W,coefw]= prog.newFreePoly(monomials([x_vars],0:degree));
[prog,V,coefv]= prog.newFreePoly(monomials([x_vars],0:degree));

% V = g*z + .5*zd^2 + .5*xd^2 + 1/8*pitchd^2;


%% Dynamics
cpi8 = cos(pi/8);
spi8 = sin(pi/8);
rt2 = sqrt(2);

f_impact=[0, 0;
    0, 0;
    0, 0;
    0, 0;
    lx(1), lx(2);
    lzsq(1)^2, lzsq(2)^2;
    4*lzsq(1)^2*(spi8*c + s*cpi8)- 4*(-lx(1))*(c*cpi8 - s*spi8), 4*lzsq(2)^2*(-spi8*c + s*cpi8)- 4*(-lx(2))*(c*cpi8 + s*spi8)];

phi{1} = z + s*spi8 - c*cpi8 + cpi8;
phi{2} = z - s*spi8 - c*cpi8 + cpi8;

phidot{1} = zd + c*pitchd*spi8 + s*pitchd*cpi8;
phidot{2} = zd - c*pitchd*spi8 + s*pitchd*cpi8;

psi{1} = xd + c*pitchd*cpi8 - s*pitchd*spi8;
psi{2} = xd + c*pitchd*cpi8 + s*pitchd*spi8;

f_free=[xd;
    zd;
    c*pitchd;
    -s*pitchd;
    0;
    -g
    0];

Vdot_free = diff(V,[x;z;s;c;xd;zd;pitchd])*f_free;
Vdot_impact_1 = diff(V,[x;z;s;c;xd;zd;pitchd])*f_impact(:,1);
Vdot_impact_2 = diff(V,[x;z;s;c;xd;zd;pitchd])*f_impact(:,2);

%% SOS functions
% (1) -Vdot_free(x) >= 0 for x admissable and x in B
% (2) -Vdot_impact_1(x,l) >= 0 for (x,l) admissable and x in B
% (3) -Vdot_impact_2(x,l) >= 0 for (x,l) admissable and x in B
% (4) W(x) - V(0,x) - 1 >= 0 for (x) admissable and in B
% (5) V(x) >= 0 for x in BT and admissable
% (6) w >= 0 for x in B

sos_1 = -Vdot_free;
sos_2 = -Vdot_impact_1;
sos_3 = -Vdot_impact_2;
sos_4 = W - V - 1;
sos_5 = V;
sos_6 = W;



%% Add in constraints
const_deg = degree;
sig = {};
coefsig = {};

% non-penetration (admissability of x)
[prog, sos_1, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_1, phi{1}, [x_vars], const_deg);
[prog, sos_1, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_1, phi{2}, [x_vars], const_deg);

[prog, sos_2, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_2, phi{2}, [x_vars;lx(1)], const_deg);

[prog, sos_3, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_3, phi{1}, [x_vars;lx(2)], const_deg);

[prog, sos_4, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_4, phi{1}, x_vars, const_deg);
[prog, sos_4, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_4, phi{2}, x_vars, const_deg);

[prog, sos_5, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_5, phi{1}, x_vars, const_deg);
[prog, sos_5, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_5, phi{2}, x_vars, const_deg);

% % trigonometric constraint
% [prog, sos_1, sig{end+1}, coefsig{end+1}] = spotless_add_eq_sprocedure(prog, sos_1, 1 - s^2 - c^2, [x_vars], const_deg);
% [prog, sos_2, sig{end+1}, coefsig{end+1}] = spotless_add_eq_sprocedure(prog, sos_2, 1 - s^2 - c^2, [x_vars;lx(1)], const_deg);
% [prog, sos_3, sig{end+1}, coefsig{end+1}] = spotless_add_eq_sprocedure(prog, sos_3, 1 - s^2 - c^2, [x_vars;lx(2)], const_deg);
% [prog, sos_4, sig{end+1}, coefsig{end+1}] = spotless_add_eq_sprocedure(prog, sos_4, 1 - s^2 - c^2, x_vars, const_deg);
% [prog, sos_5, sig{end+1}, coefsig{end+1}] = spotless_add_eq_sprocedure(prog, sos_5, 1 - s^2 - c^2, x_vars, const_deg);
% [prog, sos_6, sig{end+1}, coefsig{end+1}] = spotless_add_eq_sprocedure(prog, sos_6, 1 - s^2 - c^2, x_vars, const_deg);

% Ball constraints
h_B = 1 - [z;qd]'*A*[z;qd];
h_BT = 1 - [z;qd]'*AT*[z;qd];

[prog, sos_1, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_1, h_B, [x_vars], const_deg);
[prog, sos_2, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_2, h_B, [x_vars;lx(1)], const_deg);
[prog, sos_3, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_3, h_B, [x_vars;lx(2)], const_deg);
[prog, sos_4, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_4, h_B, x_vars, const_deg);
[prog, sos_5, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_5, h_BT, x_vars, const_deg);
[prog, sos_6, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_6, h_B, x_vars, const_deg);

% Angle constraints
[prog, sos_1, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_1, c - cos(theta_bound), [x_vars], const_deg);
[prog, sos_2, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_2, c - cos(theta_bound), [x_vars;lx(1)], const_deg);
[prog, sos_3, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_3, c - cos(theta_bound), [x_vars;lx(2)], const_deg);
[prog, sos_4, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_4, c - cos(theta_bound), x_vars, const_deg);
[prog, sos_5, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_5, c - cos(thetaT_bound), x_vars, const_deg);
[prog, sos_6, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_6, c - cos(theta_bound), x_vars, const_deg);

[prog, sos_1, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_1, sin(theta_bound)^2 - s^2, [x_vars], const_deg);
[prog, sos_2, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_2, sin(theta_bound)^2 - s^2, [x_vars;lx(1)], const_deg);
[prog, sos_3, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_3, sin(theta_bound)^2 - s^2, [x_vars;lx(2)], const_deg);
[prog, sos_4, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_4, sin(theta_bound)^2 - s^2, x_vars, const_deg);
[prog, sos_5, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_5, sin(thetaT_bound)^2 - s^2, x_vars, const_deg);
[prog, sos_6, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_6, sin(theta_bound)^2 - s^2, x_vars, const_deg);

% Contact constraints (admissability of lambda)
[prog, sos_2, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_2, -phidot{1}, [x_vars;lx(1)], const_deg);
[prog, sos_2, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_2, -lx(1)*psi{1}, [x_vars;lx(1)], const_deg);
[prog, sos_2, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_2, lzsq(1)^2 - lx(1)^2, [x_vars;lx(1)], const_deg);
[prog, sos_2, sig{end+1}, coefsig{end+1}] = spotless_add_eq_sprocedure(prog, sos_2, phi{1}, [x_vars;lx(1)], const_deg);
[prog, sos_2, sig{end+1}, coefsig{end+1}] = spotless_add_eq_sprocedure(prog, sos_2, (lzsq(1)^2 - lx(1)^2)*psi{1}, [x_vars;lx(1)], const_deg);  %should this be psi^2?

[prog, sos_3, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_3, -phidot{2}, [x_vars;lx(2)], const_deg);
[prog, sos_3, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_3, -lx(2)*psi{2}, [x_vars;lx(2)], const_deg);
[prog, sos_3, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_3, lzsq(2)^2 - lx(2)^2, [x_vars;lx(2)], const_deg);
[prog, sos_3, sig{end+1}, coefsig{end+1}] = spotless_add_eq_sprocedure(prog, sos_3, phi{2}, [x_vars;lx(2)], const_deg);
[prog, sos_3, sig{end+1}, coefsig{end+1}] = spotless_add_eq_sprocedure(prog, sos_3, (lzsq(2)^2 - lx(2)^2)*psi{2}, [x_vars;lx(2)], const_deg);  %should this be psi^2?

%% Setup cost function
[vars,alphas,coeff] = decomp(W, prog.freeVar);
assert(isequal(vars(4),s))
assert(isequal(vars(6),c))
assert(isequal(qd,vars([1 3 5])))
assert(isequal(z,vars(2)))

%first, integrate angle
for i=1:size(alphas,1),
  l_theta(i,1) = trigMonomIntegral(alphas(i,[6 4]),-theta_bound,theta_bound);
end
coeff = coeff.*l_theta';
alphas = alphas(:,[2 1 3 5]);

nX = 4;

% [~,alphas] = monomials(x_vars,0:degree);
betas = 0.5*(alphas + 1);
Ra = (R.^(sum(alphas,2) + nX))./(sum(alphas,2) + nX);
IS = 2*prod(gamma(betas),2)./(gamma(sum(betas,2)));
l = Ra.*IS;
alphaszero = (mod(alphas,2) ~= 0);
alphaszero = any(alphaszero,2);
l(alphaszero) = 0;
l = l.*prod(repmat(R_diag,size(alphas,1),1).^(alphas+1),2);

% coeff = ones(size(coeff));

%% Solve program
prog = prog.withSOS(sos_1);
prog = prog.withSOS(sos_2);
prog = prog.withSOS(sos_3);
prog = prog.withSOS(sos_4);
prog = prog.withSOS(sos_5);
prog = prog.withSOS(sos_6);

options = spotprog.defaultOptions;
options.verbose = 1;
options.trig.enable = true;
options.trig.sin = s;
options.trig.cos = c;
options.do_fr = true;
sol = prog.minimize(coeff*l,@spot_mosek,options);

%% Plotting
figure
hold off
[Z,THETA] = meshgrid(linspace(0,R_diag(1),100),linspace(-theta_bound,theta_bound,100));
L=prod(size(Z));
Wval = msubs(sol.eval(W),[q;qd],[zeros(1,L);Z(:)';sin(THETA(:))';cos(THETA(:))';zeros(3,L)]);
Wval=reshape(Wval,size(Z,1),[]);
[cl,h]=contour(THETA,Z,Wval);
clabel(cl,h);

theta = linspace(-theta_bound,theta_bound,100);
z_phi = -( - (8321567036706119*cos(theta))/9007199254740992 - (215431620425035*sin(abs(theta)))/562949953421312 + 1040195879588265/1125899906842624);

hold on
plot(theta,z_phi)

figure
Vval = msubs(sol.eval(V),[q;qd],[zeros(1,L);Z(:)';sin(THETA(:))';cos(THETA(:))';zeros(3,L)]);
Vval=reshape(Vval,size(Z,1),[]);
[cl,h]=contour(THETA,Z,Vval);
clabel(cl,h);

hold on
plot(theta,z_phi)

figure
Vdfval = msubs(sol.eval(Vdot_free),[q;qd],[zeros(1,L);Z(:)';sin(THETA(:))';cos(THETA(:))';zeros(3,L)]);
Vdfval=reshape(Vdfval,size(Z,1),[]);
[cl,h]=contour(THETA,Z,Vdfval);
clabel(cl,h);

hold on
plot(theta,z_phi)
