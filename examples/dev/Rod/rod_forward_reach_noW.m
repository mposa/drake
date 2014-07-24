megaclear
degree = 4;

T = .1;
impact_time_scale = .01;

R = 1; % don't change this
R_diag = [.2 .2 2 2 2];

theta_bound = .5;
A = 1/R*diag(1./(R_diag.^2));

g = 9.81;

prog = spotsosprog();

rbmoptions.floating = true;
rbmoptions.ignore_self_collisions = true;
rbmoptions.terrain = RigidBodyFlatTerrain();
p = PlanarRigidBodyManipulator('Rod.urdf',rbmoptions);

sos_option = 1;
switch sos_option
  case 1
    withSOS_fun = @withSOS;
    PSD_fun = @newPSD;
  case 2
    withSOS_fun = @withSDSOS;
    PSD_fun = @newSDD;
  case 3
    withSOS_fun = @withDSOS;
    PSD_fun = @newDD;
end

%% Add indeterminate variables
q_mss = msspoly('q',3);
s_mss = msspoly('s',3);
c_mss = msspoly('c',3);
q = [q_mss(1:2);s_mss(3);c_mss(3)];
qd = msspoly('v',3);
lx = msspoly('lx',2);
t = msspoly('t',1);
lzsq = [1;1];

x = q(1);
z = q(2);
s = q(3);
c = q(4);

xd = qd(1);
zd = qd(2);
pitchd = qd(3);

V_vars = [q;qd;t];
V_vars_not = [q;qd];
V_vars_sub = [x;s;c;qd;t];

prog = prog.withIndeterminate(q);
prog = prog.withIndeterminate(qd);
prog = prog.withIndeterminate(lx);
prog = prog.withIndeterminate(t);

%% Functions V, W
[prog,V,coefv]= prog.newFreePoly(monomials(V_vars,0:degree));

% linear constraints on V0
N_sample = 100;
th_range = .1;
x0 = [0;.049;0;0;0;0];
sample_range = [.02;.02;.1;.02;.02;.02];
for i=1:N_sample,
  x_sample = (rand(6,1) - .5)*2.*sample_range + x0;
  sample_sub = [x_sample(1:2); sin(x_sample(3)); cos(x_sample(3)); x_sample(4:6)];
  prog = prog.withPos(subs(V,[t;q;qd],[0;sample_sub]) - 1);
end

% prog = prog.withPos(subs(V,[t;q;qd],[0;0;.04;0;1;0;0;0]) - 1);
% prog = prog.withPos(subs(V,[t;q;qd],[0;0;.06;0;1;0;0;0]) - 1);

%% Dynamics
[H,Hf_free,Hf_impact,phi,phidot,psi] = getTrigEOM(p);
H = double(H);

f_free = [xd;zd;c*pitchd;-s*pitchd;inv(H)*Hf_free];

f_impact = [zeros(4,2); inv(H)*Hf_impact];

Vdot_free = T*diff(V,[x;z;s;c;xd;zd;pitchd])*f_free + diff(V,t);
Vdot_impact_1 = diff(V,[x;z;s;c;xd;zd;pitchd])*f_impact(:,1) + impact_time_scale*diff(V,t);
Vdot_impact_2 = diff(V,[x;z;s;c;xd;zd;pitchd])*f_impact(:,2) + impact_time_scale*diff(V,t);

T = 1; % rescale time


%% SOS functions
% (1) -Vdot_free(x) >= 0 for x admissable and x in B
% (2) -Vdot_impact_1(x,l) >= 0 for (x,l) admissable and x in B
% (3) -Vdot_impact_2(x,l) >= 0 for (x,l) admissable and x in B
% (5) V(0,x) >= 0 for x in B

% (linear constraints) V(0,x) >= mu_0 for x in B and admissable

sos_1 = Vdot_free;
sos_2 = Vdot_impact_1;
sos_3 = Vdot_impact_2;
sos_4 = subs(V,t,0);  % doesn't work. sub in T and it does something...but what? 

%% Add in constraints
const_deg = degree;
sig = {};
coefsig = {};

% non-penetration (admissability of x)
[prog, sos_1, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_1, phi(1), [V_vars], const_deg, sos_option);
[prog, sos_1, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_1, phi(2), [V_vars], const_deg, sos_option);

[prog, sos_2, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_2, phi(2), [V_vars_sub;lx(1)], const_deg, sos_option);

[prog, sos_3, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_3, phi(1), [V_vars_sub;lx(2)], const_deg, sos_option);

% Ball constraints
h_B = 1 - [x;z;qd]'*A*[x;z;qd];

[prog, sos_1, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_1, h_B, [V_vars], const_deg, sos_option);
[prog, sos_2, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_2, h_B, [V_vars_sub;lx(1)], const_deg, sos_option);
[prog, sos_3, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_3, h_B, [V_vars_sub;lx(2)], const_deg, sos_option);
[prog, sos_4, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_4, h_B, V_vars_not, const_deg, sos_option);

% Angle constraints
[prog, sos_1, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_1, c - cos(theta_bound), [V_vars], const_deg, sos_option);
[prog, sos_2, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_2, c - cos(theta_bound), [V_vars_sub;lx(1)], const_deg, sos_option);
[prog, sos_3, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_3, c - cos(theta_bound), [V_vars_sub;lx(2)], const_deg, sos_option);
[prog, sos_4, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_4, c - cos(theta_bound), V_vars_not, const_deg, sos_option);

% Contact constraints (admissability of lambda)
[prog, sos_2, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_2, -phidot(1), [V_vars_sub;lx(1)], const_deg, sos_option);
[prog, sos_2, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_2, -lx(1)*psi(1), [V_vars_sub;lx(1)], const_deg-2, sos_option);
[prog, sos_2, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_2, lzsq(1)^2 - lx(1)^2, [V_vars_sub;lx(1)], const_deg, sos_option);
% [prog, sos_2, sig{end+1}, coefsig{end+1}] = spotless_add_eq_sprocedure(prog, sos_2, phi(1), [V_vars_sub;lx(1)], const_deg+1);
[prog, sos_2, sig{end+1}, coefsig{end+1}] = spotless_add_eq_sprocedure(prog, sos_2, (lzsq(1)^2 - lx(1)^2)*psi(1), [V_vars_sub;lx(1)], const_deg-2);  %should this be psi^2?

[prog, sos_3, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_3, -phidot(2), [V_vars_sub;lx(2)], const_deg, sos_option);
[prog, sos_3, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_3, -lx(2)*psi(2), [V_vars_sub;lx(2)], const_deg-2, sos_option);
[prog, sos_3, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_3, lzsq(2)^2 - lx(2)^2, [V_vars_sub;lx(2)], const_deg, sos_option);
% [prog, sos_3, sig{end+1}, coefsig{end+1}] = spotless_add_eq_sprocedure(prog, sos_3, phi(2), [V_vars_sub;lx(2)], const_deg+1);
[prog, sos_3, sig{end+1}, coefsig{end+1}] = spotless_add_eq_sprocedure(prog, sos_3, (lzsq(2)^2 - lx(2)^2)*psi(2), [V_vars_sub;lx(2)], const_deg-2);  %should this be psi^2?

sos_2 = subs(sos_2,z,z-phi(1));
sos_3 = subs(sos_3,z,z-phi(2));

%% Setup cost function
[vars,alphas,coeff] = decomp(subs(V,t,T), prog.freeVar);
assert(isequal(vars(6),s))
assert(isequal(vars(5),c))
assert(isequal(qd,vars([2 4 7])))
assert(isequal(z,vars(3)))
assert(isequal(x,vars(1)))

%first, integrate angle
for i=1:size(alphas,1),
  l_theta(i,1) = trigMonomIntegral(alphas(i,[5 6]),-theta_bound,theta_bound);
end
coeff = coeff.*l_theta';
alphas = alphas(:,[1 3 2 4 7]);

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
prog = withSOS_fun(prog,sos_1);
prog = withSOS_fun(prog,sos_2);
prog = withSOS_fun(prog,sos_3);
prog = withSOS_fun(prog,sos_4);

options = spotprog.defaultOptions;
options.verbose = 1;
options.trig.enable = true;
options.trig.sin = s;
options.trig.cos = c;
options.do_fr = false;
options.solver_options.mosek.MSK_IPAR_INTPNT_BASIS = 'MSK_BI_NEVER';
sol = prog.minimize(coeff*l,@spot_mosek,options);


%% Plotting
close all
% figure
% hold off
[Z,THETA] = meshgrid(linspace(-R_diag(1),R_diag(1),201),linspace(-theta_bound,theta_bound,201));
L=prod(size(Z));

% qd_data = repmat(zeros(3,1),1,L);
% Wval = msubs(sol.eval(W),[q;qd],[zeros(1,L);Z(:)';sin(THETA(:))';cos(THETA(:))';qd_data]);
% Wval=reshape(Wval,size(Z,1),[]);
% [cl,h]=contour(THETA,Z,Wval);
% clabel(cl,h);
% xlabel('theta');
% ylabel('z');
% title('W')

C = cos(THETA);
S = sin(THETA);
phi1val = msubs(subs(phi(1),x,0),[z;s;c],[Z(:) S(:) C(:)]');
phi2val = msubs(subs(phi(2),x,0),[z;s;c],[Z(:) S(:) C(:)]');
phival = min(phi1val,phi2val);
phival = reshape(phival,size(C,1),[]);

% theta = linspace(-theta_bound,theta_bound,100);
% z_phi = -( - (8321567036706119*cos(theta))/9007199254740992 - (215431620425035*sin(abs(theta)))/562949953421312 + 1040195879588265/1125899906842624);

hold on
% plot(theta,z_phi)
contour(THETA,Z,phival,[0 0],'k','Linewidth',3)

figure
Vval = msubs(sol.eval(subs(V,t,0)),[q;qd],[zeros(1,L);Z(:)';sin(THETA(:))';cos(THETA(:))';qd_data]);
Vval=reshape(Vval,size(Z,1),[]);
[cl,h]=contour(THETA,Z,Vval);
clabel(cl,h);

hold on
[cl,h]=contour(THETA,Z,Vval,[1 1],'b','Linewidth',2);
% plot(theta,z_phi)
contour(THETA,Z,phival,[0 0],'k','Linewidth',3)
xlabel('theta');
ylabel('z');
title('V(0,x)')

figure
Vval = msubs(sol.eval(subs(V,t,T)),[q;qd],[zeros(1,L);Z(:)';sin(THETA(:))';cos(THETA(:))';qd_data]);
Vval=reshape(Vval,size(Z,1),[]);
[cl,h]=contour(THETA,Z,Vval);
clabel(cl,h);

hold on
[cl,h]=contour(THETA,Z,Vval,[1 1],'b','Linewidth',2);

% plot(theta,z_phi)
contour(THETA,Z,phival,[0 0],'k','Linewidth',3)
xlabel('theta');
ylabel('z');
title('V(T,x)')

figure
[Z,ZD] = meshgrid(linspace(-R_diag(2)*0,R_diag(2),201),linspace(-R_diag(4),R_diag(4),201));

Vval = msubs(sol.eval(subs(V,t,T)),[q;qd],[zeros(1,L); Z(:)'; zeros(1,L); ones(1,L);zeros(1,L); ZD(:)'; zeros(1,L)]);
Vval=reshape(Vval,size(Z,1),[]);
[cl,h]=contour(Z,ZD,Vval);
clabel(cl,h);

hold on
[cl,h]=contour(Z,ZD,Vval,[1 1],'b','Linewidth',2);

% plot(theta,z_phi)
% contour(THETA,Z,phival,[0 0],'k','Linewidth',3)
xlabel('z');
ylabel('zd');
title('V(T,x)')


figure
[XD,ZD] = meshgrid(linspace(-R_diag(3),R_diag(3),201),linspace(-R_diag(4),R_diag(4),201));

Vval = msubs(sol.eval(subs(V,t,T)),[q;qd],[zeros(3,L);ones(1,L);XD(:)'; ZD(:)'; zeros(1,L)]);
Vval=reshape(Vval,size(Z,1),[]);
[cl,h]=contour(XD,ZD,Vval);
clabel(cl,h);

hold on
[cl,h]=contour(XD,ZD,Vval,[1 1],'b','Linewidth',2);

% plot(theta,z_phi)
% contour(THETA,Z,phival,[0 0],'k','Linewidth',3)
xlabel('xd');
ylabel('zd');
title('V(T,x)')

figure
ZD = linspace(-R_diag(4),R_diag(4),201)';
LL = length(ZD);
Vval = msubs(sol.eval(subs(V,t,T)),[q;qd],[zeros(3,LL);ones(1,LL);zeros(1,LL); ZD(:)'; zeros(1,LL)]);
plot(ZD,Vval)
xlabel('zd')
ylabel('V')
title('z=0')
ylim([0 2])

% figure
% Vdfval = msubs(sol.eval(subs(Vdot_free,t,0)),[q;qd],[zeros(1,L);Z(:)';sin(THETA(:))';cos(THETA(:))';qd_data]);
% Vdfval=reshape(Vdfval,size(Z,1),[]);
% [cl,h]=contour(THETA,Z,Vdfval);
% clabel(cl,h);
% 
% hold on
% % plot(theta,z_phi)
% contour(THETA,Z,phival,[0 0],'k','Linewidth',3)
