megaclear

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
sos_fun= @spot_mosek;
degree = 4;

ground_angle = 0;

finite_time = true;
T = 10;
impact_time_scale = .1;

R = 1; % don't change this
R_diag = [.2 2 1 1];
RT_diag = R_diag;
% R_diag = [1 1 1 1];

RT = .05;
theta_bound = .75;%pi/8;%.75;
thetaT_bound = .01;
A = 1/R*diag(1./(R_diag.^2));
AT = 1/RT*diag(1./(RT_diag.^2));
% AT = 1/RT*diag([1;1;1;1]);

g = 9.81;

prog = spotsosprog();

% xf_fb = [10.391;-0.0830;.1-pi/8;1.2533;.3639;1.3079]; % from sim
% xf_fb = [0;0;0;1.2854;-.7264;1.6393];
% xf = [xf_fb(1);xf_fb(2) - cos(pi/8);sin(xf_fb(3) + pi/8);cos(xf_fb(3) + pi/8);xf_fb(4:6)];

% aT = [0;cos(ground_angle)*xf_fb(4) - sin(ground_angle)*xf_fb(5);sin(ground_angle)*xf_fb(4) + cos(ground_angle)*xf_fb(5);xf_fb(6)]
% aT = zeros(4,1);
aT = zeros(4,1);
%...figure out these conditions
%% Add indeterminate variables
q = msspoly('q',4);
qd = msspoly('qd',3);
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

W_vars = [q(2:end);qd];
V_vars = [q(2:end);qd];

if finite_time
  V_vars = [V_vars;t];
end

prog = prog.withIndeterminate(q(2:end));
prog = prog.withIndeterminate(qd);
prog = prog.withIndeterminate(lx);

if finite_time
  prog = prog.withIndeterminate(t);
end

%% Functions V, W
[prog,W,coefw]= prog.newFreePoly(monomials(W_vars,0:degree));
[prog,V,coefv]= prog.newFreePoly(monomials(V_vars,0:degree));

% V = g*z + .5*zd^2 + .5*xd^2 + 1/8*pitchd^2;

% Make V and W cyclic
% prog = prog.withEqs(decomp(subs(V,[s;c],[sin(pi/8);cos(pi/8)]) - subs(V,[s;c],[sin(-pi/8);cos(-pi/8)]),[t;q;qd]));
% prog = prog.withEqs(decomp(subs(W,[s;c],[sin(pi/8);cos(pi/8)]) - subs(W,[s;c],[sin(-pi/8);cos(-pi/8)]),[t;q;qd]));


%% Dynamics

% do dynamics the right way
cpi8 = cos(pi/8);
spi8 = sin(pi/8);
rt2 = sqrt(2);
f_free=[xd;
    zd;
    c*pitchd;
    -s*pitchd;
    g*sin(ground_angle);
    -g*cos(ground_angle);
    0];
  
fk_foot{1} = [x + s*cpi8 - c*spi8;z + s*spi8 - c*cpi8 + cpi8];
fk_foot{2} = [x + s*cpi8 + c*spi8;z - s*spi8 - c*cpi8 + cpi8];

phi{1} = fk_foot{1}(2,:);
phi{2} = fk_foot{2}(2,:);

J_fk{1} = diff(fk_foot{1},q);
J_fk{2} = diff(fk_foot{2},q);

vk_foot{1} = J_fk{1}*f_free(1:4);
vk_foot{2} = J_fk{2}*f_free(1:4);

phidot{1} = vk_foot{1}(2);
phidot{2} = vk_foot{2}(2);

psi{1} = vk_foot{1}(1);
psi{2} = vk_foot{2}(1);

f_impact = [zeros(4,2); diag([1;1;4])*[1 0 0 0;0 1 0 0; 0 0 c -s]*[J_fk{1}'*[lx(1);lzsq(1)^2] J_fk{2}'*[lx(2);lzsq(2)^2]]]; 


Vdot_free = T*diff(V,[x;z;s;c;xd;zd;pitchd])*f_free + diff(V,t);
Vdot_impact_1 = diff(V,[x;z;s;c;xd;zd;pitchd])*f_impact(:,1) + impact_time_scale*diff(V,t);
Vdot_impact_2 = diff(V,[x;z;s;c;xd;zd;pitchd])*f_impact(:,2) + impact_time_scale*diff(V,t);

T = 1; % rescale time

%% SOS functions
% (1) -Vdot_free(x) >= 0 for x admissable and x in B
% (2) -Vdot_impact_1(x,l) >= 0 for (x,l) admissable and x in B
% (3) -Vdot_impact_2(x,l) >= 0 for (x,l) admissable and x in B
% (4) W(x) - V(0,x) - 1 >= 0 for (x) admissable and in B
% (5) V(x,T) >= 0 for x in BT and admissable
% (6) w >= 0 for x in B

sos_1 = -Vdot_free;
sos_2 = -Vdot_impact_1;
sos_3 = -Vdot_impact_2;
sos_4 = W - subs(V,t,0) - 1;
sos_5 = V; %subs(V,t,T);  %use V for flexible finite time
sos_6 = W;

V_vars_c = V_vars(2:end);



%% Add in constraints
const_deg = degree;
sig = {};
coefsig = {};

% non-penetration (admissability of x)
[prog, sos_1, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_1, phi{1}, [V_vars], const_deg, sos_option);
[prog, sos_1, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_1, phi{2}, [V_vars], const_deg, sos_option);

[prog, sos_2, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_2, phi{2}, [V_vars_c;lx(1)], const_deg, sos_option);

[prog, sos_3, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_3, phi{1}, [V_vars_c;lx(2)], const_deg, sos_option);

[prog, sos_4, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_4, phi{1}, W_vars, const_deg, sos_option);
[prog, sos_4, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_4, phi{2}, W_vars, const_deg, sos_option);

[prog, sos_5, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_5, phi{1}, W_vars, const_deg, sos_option);
[prog, sos_5, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_5, phi{2}, W_vars, const_deg, sos_option);

% % trigonometric constraint
% [prog, sos_1, sig{end+1}, coefsig{end+1}] = spotless_add_eq_sprocedure(prog, sos_1, 1 - s^2 - c^2, [x_vars], const_deg);
% [prog, sos_2, sig{end+1}, coefsig{end+1}] = spotless_add_eq_sprocedure(prog, sos_2, 1 - s^2 - c^2, [x_vars;lx(1)], const_deg);
% [prog, sos_3, sig{end+1}, coefsig{end+1}] = spotless_add_eq_sprocedure(prog, sos_3, 1 - s^2 - c^2, [x_vars;lx(2)], const_deg);
% [prog, sos_4, sig{end+1}, coefsig{end+1}] = spotless_add_eq_sprocedure(prog, sos_4, 1 - s^2 - c^2, x_vars, const_deg);
% [prog, sos_5, sig{end+1}, coefsig{end+1}] = spotless_add_eq_sprocedure(prog, sos_5, 1 - s^2 - c^2, x_vars, const_deg);
% [prog, sos_6, sig{end+1}, coefsig{end+1}] = spotless_add_eq_sprocedure(prog, sos_6, 1 - s^2 - c^2, x_vars, const_deg);

% Ball constraints
h_B = 1 - [z;qd]'*A*[z;qd];
h_BT = 1 - ([z;qd]-aT)'*AT*([z;qd]-aT);

[prog, sos_1, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_1, h_B, [V_vars], const_deg, sos_option);
[prog, sos_2, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_2, h_B, [V_vars_c;lx(1)], const_deg, sos_option);
[prog, sos_3, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_3, h_B, [V_vars_c;lx(2)], const_deg, sos_option);
[prog, sos_4, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_4, h_B, W_vars, const_deg, sos_option);
[prog, sos_5, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_5, h_BT, W_vars, const_deg, sos_option);
[prog, sos_6, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_6, h_B, W_vars, const_deg, sos_option);

% Angle constraints
if theta_bound ~= pi
[prog, sos_1, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_1, c - cos(theta_bound), [V_vars], const_deg, sos_option);
[prog, sos_2, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_2, c - cos(theta_bound), [V_vars_c;lx(1)], const_deg, sos_option);
[prog, sos_3, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_3, c - cos(theta_bound), [V_vars_c;lx(2)], const_deg, sos_option);
[prog, sos_4, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_4, c - cos(theta_bound), W_vars, const_deg, sos_option);
[prog, sos_5, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_5, c - cos(thetaT_bound), W_vars, const_deg, sos_option);
[prog, sos_6, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_6, c - cos(theta_bound), W_vars, const_deg, sos_option);
end

% [prog, sos_1, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_1, sin(theta_bound)^2 - s^2, [x_vars], const_deg);
% [prog, sos_2, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_2, sin(theta_bound)^2 - s^2, [x_vars;lx(1)], const_deg);
% [prog, sos_3, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_3, sin(theta_bound)^2 - s^2, [x_vars;lx(2)], const_deg);
% [prog, sos_4, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_4, sin(theta_bound)^2 - s^2, x_vars, const_deg);
% [prog, sos_5, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_5, sin(thetaT_bound)^2 - s^2, x_vars, const_deg);
% [prog, sos_6, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_6, sin(theta_bound)^2 - s^2, x_vars, const_deg);

% Contact constraints (admissability of lambda)
[prog, sos_2, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_2, -phidot{1}, [V_vars_c;lx(1)], const_deg, sos_option);
[prog, sos_2, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_2, -lx(1)*psi{1}, [V_vars_c;lx(1)], const_deg-1, sos_option);
[prog, sos_2, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_2, lzsq(1)^2 - lx(1)^2, [V_vars_c;lx(1)], const_deg, sos_option);
% [prog, sos_2, sig{end+1}, coefsig{end+1}] = spotless_add_eq_sprocedure(prog, sos_2, phi{1}, [V_vars_c;lx(1)], const_deg+1);
[prog, sos_2, sig{end+1}, coefsig{end+1}] = spotless_add_eq_sprocedure(prog, sos_2, (lzsq(1)^2 - lx(1)^2)*psi{1}, [V_vars_c;lx(1)], const_deg-2);  %should this be psi^2?

[prog, sos_3, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_3, -phidot{2}, [V_vars_c;lx(2)], const_deg, sos_option);
[prog, sos_3, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_3, -lx(2)*psi{2}, [V_vars_c;lx(2)], const_deg-1, sos_option);
[prog, sos_3, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_3, lzsq(2)^2 - lx(2)^2, [V_vars_c;lx(2)], const_deg, sos_option);
% [prog, sos_3, sig{end+1}, coefsig{end+1}] = spotless_add_eq_sprocedure(prog, sos_3, phi{2}, [V_vars_c;lx(2)], const_deg+1);
[prog, sos_3, sig{end+1}, coefsig{end+1}] = spotless_add_eq_sprocedure(prog, sos_3, (lzsq(2)^2 - lx(2)^2)*psi{2}, [V_vars_c;lx(2)], const_deg-2);  %should this be psi^2?

% time constraints
[prog, sos_1, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_1, T*t - t^2, [V_vars], const_deg, sos_option);
[prog, sos_2, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_2, T*t - t^2, [V_vars_c;lx(1)], const_deg, sos_option);
[prog, sos_3, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_3, T*t - t^2, [V_vars_c;lx(2)], const_deg, sos_option);
[prog, sos_4, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_4, T*t - t^2, W_vars, const_deg, sos_option);
[prog, sos_5, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_5, T*t - t^2, W_vars, const_deg, sos_option);
[prog, sos_6, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_6, T*t - t^2, W_vars, const_deg, sos_option);

sos_2 = subs(sos_2,z,z-phi{1});
sos_3 = subs(sos_3,z,z-phi{2});
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
prog = withSOS_fun(prog,sos_1);
prog = withSOS_fun(prog,sos_2);
prog = withSOS_fun(prog,sos_3);
prog = withSOS_fun(prog,sos_4);
prog = withSOS_fun(prog,sos_5);
prog = withSOS_fun(prog,sos_6);

options = spotprog.defaultOptions;
options.verbose = 1;
options.trig.enable = true;
options.trig.sin = s;
options.trig.cos = c;
options.do_fr = true;
options.solveroptions.MSK_IPAR_INTPNT_BASIS = 'MSK_BI_NEVER';
tic
sol = prog.minimize(coeff*l,sos_fun,options);
toc

%% Plotting
close all
figure
hold off
[Z,THETA] = meshgrid(linspace(-R_diag(1),R_diag(1),100),linspace(-theta_bound,theta_bound,100));
L=prod(size(Z));

qd_data = repmat(aT(2:4),1,L);
Wval = msubs(sol.eval(W),[q;qd],[zeros(1,L);Z(:)';sin(THETA(:))';cos(THETA(:))';qd_data]);
Wval=reshape(Wval,size(Z,1),[]);
[cl,h]=contour(THETA,Z,Wval);
clabel(cl,h);

C = cos(THETA);
S = sin(THETA);
phi1val = msubs(subs(phi{1},x,0),[z;s;c],[Z(:) S(:) C(:)]');
phi2val = msubs(subs(phi{2},x,0),[z;s;c],[Z(:) S(:) C(:)]');
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
% plot(theta,z_phi)
contour(THETA,Z,phival,[0 0],'k','Linewidth',3)


figure
Vdfval = msubs(sol.eval(subs(Vdot_free,t,0)),[q;qd],[zeros(1,L);Z(:)';sin(THETA(:))';cos(THETA(:))';qd_data]);
Vdfval=reshape(Vdfval,size(Z,1),[]);
[cl,h]=contour(THETA,Z,Vdfval);
clabel(cl,h);

hold on
% plot(theta,z_phi)
contour(THETA,Z,phival,[0 0],'k','Linewidth',3)
