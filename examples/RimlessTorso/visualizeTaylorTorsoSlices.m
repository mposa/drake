% load(datapath('torso_cubic_controller_taylor_iter_39'))
% load(datapath('zscale_torso_cubic_controller_taylor_iter_4'))
load(datapath('zscale_skinny_cubic_controller_taylor_iter_70'));
load skinny_taylor_eom
% load torso_taylor_eom_fix
z_scale = .1;

q = msspoly('q',4);
qd = msspoly('qd',4);
u = msspoly('u',1);
lx = msspoly('lx',2);
x = q(1);
z = q(2);
pitch = q(3);
theta = q(4);
thetad = qd(4);
xd = qd(1);
zd = qd(2);

taylor_degree = 3;
taylor_vars = [q;qd;u;lx];

f_impact_1 = clean(getmsspoly(f_impact_1,taylor_vars,taylor_degree));
f_impact_2 = clean(getmsspoly(f_impact_2,taylor_vars,taylor_degree));
f_free = clean(getmsspoly(f_free,taylor_vars,taylor_degree));
phi = clean(getmsspoly(phi,taylor_vars,taylor_degree-1));
phidot = clean(getmsspoly(phidot,taylor_vars,taylor_degree-1));
psi = clean(getmsspoly(psi,taylor_vars,taylor_degree));
E = clean(getmsspoly(E,taylor_vars,taylor_degree));

% phi = subs(phi,[z;zd],[z;zd]*z_scale);
Vsol = subs(Vsol,[z;zd],[z;zd]/z_scale);
TA = diag([1/z_scale;ones(6,1)]);
AO = TA'*AO*TA;
AI = TA'*AI*TA;

K = [10 1];

f_free = subs(f_free,u,-K*[theta;thetad]);
Vdot_free = diff(Vsol,q)*qd + diff(Vsol,qd)*f_free;
Vdot_impact_1 = diff(Vsol,qd)*f_impact_1;
Vdot_impact_2 = diff(Vsol,qd)*f_impact_2;


pitch_val =  -.5:.05:.5;
z_range = 0:.01:.2;
theta_val =  -.5:.05:.5;
%%
% load iter_4
% rho_i = R;
% rho_o = 1;
% load(datapath('torso_taylor_iter_28'))

% load torso_ff_03
% load skinny_ff_06



rho_i = R;
rho_o = D;
Ao2 = AO;
Ai = AI;

ball_vec = [z;pitch;0;zeros(4,1)];

h_Bo = ball_vec'*Ao2*ball_vec;
% h_Bo2 = .1 - (2-c-c_th) - z^2 - .05*.5*qd'*H*qd;

%searched for V with .01, worked
h_Bi = ball_vec'*Ai*ball_vec; %worked with .01 and E, but failed sdsos
Vsub = subs(Vsol,[theta;qd],[zeros(5,1)]);
phisub = subs(phi,[theta;qd],[zeros(5,1)]);
Vdot_freesub = subs(Vdot_free,[theta;qd],[zeros(5,1)]);


[PITCH,Z] = meshgrid(pitch_val,z_range);
C = cos(PITCH);
S = sin(PITCH);

Vval = dmsubs(Vsub,[z;pitch],[Z(:) PITCH(:)]');
Vval = reshape(Vval,size(Z,1),[]);

BIval = dmsubs(h_Bi,[z;pitch],[Z(:) PITCH(:)]');
BIval = reshape(BIval,size(Z,1),[]);

BOval = dmsubs(h_Bo,[z;pitch],[Z(:) PITCH(:)]');
BOval = reshape(BOval,size(Z,1),[]);

PHIval = min(dmsubs(phisub,[z;pitch],[Z(:) PITCH(:)]'));
PHIval = reshape(PHIval,size(Z,1),[]);

Vdotval = dmsubs(Vdot_freesub,[z;pitch],[Z(:) PITCH(:)]');
Vdotval = reshape(Vdotval,size(Z,1),[]);

figure(1)
hold off
[cl, h] = contour(PITCH,Z,Vval,[1 1]);
clabel(cl,h);
hold on
[cl, h] = contour(PITCH,Z,BIval,[rho_i rho_i]);
clabel(cl,h);

[cl, h] = contour(PITCH,Z,BOval,[rho_o rho_o]);
clabel(cl,h);

[cl, h] = contour(PITCH,Z,PHIval,[0 0]);
clabel(cl,h);

% figure(4)
% hold off
% [cl, h] = contour(PITCH,Z,Vdotval);
% clabel(cl,h);
% hold on
% [cl, h] = contour(PITCH,Z,PHIval,[0 0]);
% clabel(cl,h);


%%
ball_vec = [z;0;theta;zeros(4,1)];
% Ao2(4,4) = 1;
% Ao2(5,5) = 1;

h_Bo = ball_vec'*Ao2*ball_vec;
% h_Bo2 = .1 - (2-c-c_th) - z^2 - .05*.5*qd'*H*qd;

%searched for V with .01, worked
h_Bi = ball_vec'*Ai*ball_vec; %worked with .01 and E, but failed sdsos
Vsub = subs(Vsol,[pitch;qd],[0;zeros(4,1)]);

[THETA,Z] = meshgrid(theta_val,z_range);

Vval = dmsubs(Vsub,[z;theta],[Z(:) THETA(:)]');
Vval = reshape(Vval,size(Z,1),[]);

BIval = dmsubs(h_Bi,[z;theta],[Z(:) THETA(:)]');
BIval = reshape(BIval,size(Z,1),[]);

BOval = dmsubs(h_Bo,[z;theta],[Z(:) THETA(:)]');
BOval = reshape(BOval,size(Z,1),[]);

figure(2)
hold off
[cl, h] = contour(THETA,Z,Vval,[1 1]);
clabel(cl,h);
hold on
[cl, h] = contour(THETA,Z,BIval,[rho_i rho_i]);
clabel(cl,h);

[cl, h] = contour(THETA,Z,BOval,[rho_o rho_o]);
clabel(cl,h);

%%
ball_vec = [z;pitch;theta;zeros(4,1)];

h_Bo = ball_vec'*Ao2*ball_vec;
h_Bi = ball_vec'*Ai*ball_vec;
Vsub = subs(Vsol,qd,zeros(4,1));
Vdot_freesub = subs(Vdot_free,qd,zeros(4,1));

% theta =  -.5:.05:.5;
% pitch =  -.5:.05:.5;
[PITCH,THETA,Z] = meshgrid(pitch_val,theta_val,z_range);

Vval = dmsubs(Vsub,[z;pitch;theta],[Z(:) PITCH(:) THETA(:)]');
Vval = reshape(full(Vval),size(Z,1),size(Z,2),[]);

Vdotval = dmsubs(Vdot_freesub,[z;pitch;theta],[Z(:) PITCH(:) THETA(:)]');
Vdotval = reshape(full(Vdotval),size(Z,1),size(Z,2),[]);

BIval = dmsubs(h_Bi,[z;pitch;theta],[Z(:) PITCH(:) THETA(:)]');
BIval = reshape(full(BIval),size(Z,1),size(Z,2),[]);

BOval = dmsubs(h_Bo,[z;pitch;theta],[Z(:) PITCH(:) THETA(:)]');
BOval = reshape(full(BOval),size(Z,1),size(Z,2),[]);

PHIval = min(dmsubs(phi,[z;pitch;theta],[Z(:) PITCH(:) THETA(:)]'));
PHIval = reshape(full(PHIval),size(Z,1),size(Z,2),[]);

figure(3)
close(3)
figure(3)
hold off
surf=isosurface(PITCH,THETA,Z,Vval,1);
p=patch(surf);
set(p,'facecolor','blue')
% set(p,'edgealpha',0)
alpha(.5)
xlabel('Pitch')
ylabel('Theta')
zlabel('z')
hold on
% isosurface(PITCH,THETA,Z,BIval,rho_i);
% alpha(.5)

surf=isosurface(PITCH,THETA,Z,PHIval,0);
p = patch(surf);
set(p,'facecolor','red')
set(p,'edgealpha',0)
alpha(.3);

% [PITCH_PHI,THETA_PHI] = meshgrid(pitch,theta);
% % z_phi = -( - (8321567036706119*cos(pitch))/9007199254740992 - (215431620425035*sin(abs(pitch)))/562949953421312 + 1040195879588265/1125899906842624);
% z_phi = max(-double(subs(dmsubs(phi,[s;c],[sin(pitch);cos(pitch)]),z,0)));
% 
% Z_PHI= repmat(z_phi',1,length(theta));
% 
% surf(PITCH_PHI,THETA_PHI,Z_PHI')
% alpha(.3);

%%
doSampleVdot = false;
if doSampleVdot
  ball_vec = [z;pitch;theta;qd];
  h_Bo = ball_vec'*Ao2*ball_vec;
  
  AI = double(AI);
  U = inv(chol(AI));
  N = 100000;

  x_val = zeros(7,N);
  rad_i = zeros(N,1);
  for i=1:N,
    x_val(:,i) = U*randn(7,1)*sqrt(rho_i/7);
    rad_i(i) = x_val(:,i)'*AI*x_val(:,i);
  end
   Vdotval = dmsubs(Vdot_free,[z;pitch;theta;qd],x_val);
   PHIval = dmsubs(phi,[z;pitch;theta;qd],x_val);
   BOval = dmsubs(h_Bo,[z;pitch;theta;qd],x_val);
   psival = dmsubs(psi,[z;pitch;theta;qd],x_val);
   phidotval = dmsubs(phidot,[z;pitch;theta;qd],x_val);
   lxval = -sign(psival);
   Vdotimpact1val = dmsubs(Vdot_impact_1,[z;pitch;theta;qd;lx(1)],[x_val; lxval(1,:)]);
   Vdotimpact2val = dmsubs(Vdot_impact_2,[z;pitch;theta;qd;lx(2)],[x_val; lxval(2,:)]);
   
   maxVdot = max(Vdotval(BOval < rho_o & min(PHIval) > 0))
   max(Vdotimpact1val(abs(PHIval(1,:)) < 1e-3 & BOval < rho_o & PHIval(2,:) > 0 & phidotval(1,:) < 0))
   max(Vdotimpact2val(abs(PHIval(2,:)) < 1e-3 & BOval < rho_o & PHIval(1,:) > 0 & phidotval(2,:) < 0))
   
end
