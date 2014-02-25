
q = msspoly('q',4);
qd = msspoly('qd',4);
s_vec = msspoly('s',4);
c_vec = msspoly('c',4);
[H,C,B,phi,phidot,psi,J,J_f,K,S,U] = reactionEOM_mss(q,qd,s_vec,c_vec);

x = q(1);
z = q(2);
theta = q(4);
s = s_vec(3);

c = c_vec(3);

%%
load flex2_iter_2
rho_i = R;
rho_o = 1;
% Ai = subs(AI,[s;c],[0;0;1;1]);
Ai = double(AI);
% rho_o = .15;
% rho = .125;
% load torso_data_fixed_4_2
% load torso_data_sd_0_04
ball_vec = [z;s;1-c;0;zeros(4,1)];
% Ao2 = zeros(9)*z;
% Ao2(1,1) = 10;
% Ao2(2,2) = .5;
% Ao2(3,3) = .5;
% Ao2(4,4) = 2.5/5;
% Ao2(5,5) = 2.5/5;
% Ao2(6:9,6:9) = .05*H;

h_Bo = ball_vec'*Ao2*ball_vec;

%searched for V with .01, worked
h_Bi = ball_vec'*Ai*ball_vec; %worked with .01 and E, but failed sdsos
Vsub = subs(Vsol,[theta;qd],[0;zeros(4,1)]);

pitch =  -.5:.02:.5;
[PITCH,Z] = meshgrid(pitch,0:.001:.15);
C = cos(PITCH);
S = sin(PITCH);

Vval = msubs(Vsub,[z;s;c],[Z(:) S(:) C(:)]');
Vval = reshape(Vval,size(C,1),[]);

BIval = msubs(h_Bi,[z;s;c],[Z(:) S(:) C(:)]');
BIval = reshape(BIval,size(C,1),[]);

BOval = msubs(h_Bo,[z;s;c],[Z(:) S(:) C(:)]');
BOval = reshape(BOval,size(C,1),[]);

figure(1)
hold off
[cl, h] = contour(PITCH,Z,Vval,[1 1]);
clabel(cl,h);
hold on
[cl, h] = contour(PITCH,Z,BIval,[rho_i rho_i]);
clabel(cl,h);

[cl, h] = contour(PITCH,Z,BOval,[rho_o rho_o]);
clabel(cl,h);

z_phi = -( - (8321567036706119*cos(pitch))/9007199254740992 - (215431620425035*sin(abs(pitch)))/562949953421312 + 1040195879588265/1125899906842624);
plot(pitch,z_phi,'r','Linewidth',3)

%%
ball_vec = [z;0;0;theta;zeros(4,1)];
% Ao2(4,4) = 1;
% Ao2(5,5) = 1;

h_Bo = ball_vec'*Ao2*ball_vec;

%searched for V with .01, worked
h_Bi = ball_vec'*Ai*ball_vec; %worked with .01 and E, but failed sdsos
Vsub = subs(Vsol,[s;c;qd],[0;1;zeros(4,1)]);

theta_val =  -1:.02:1;
[THETA,Z] = meshgrid(theta_val,0:.001:.15);

Vval = msubs(Vsub,[z;theta],[Z(:) THETA(:)]');
Vval = reshape(Vval,size(THETA,1),[]);

BIval = msubs(h_Bi,[z;theta],[Z(:) THETA(:)]');
BIval = reshape(BIval,size(THETA,1),[]);

BOval = msubs(h_Bo,[z;theta],[Z(:) THETA(:)]');
BOval = reshape(BOval,size(THETA,1),[]);

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
ball_vec = [z;s;1-c;theta;zeros(4,1)];

h_Bo = ball_vec'*Ao2*ball_vec;
h_Bi = ball_vec'*Ai*ball_vec;
Vsub = subs(Vsol,qd,zeros(4,1));

% theta =  -.5:.05:.5;
% pitch =  -.5:.05:.5;
[PITCH,THETA,Z] = meshgrid(pitch,theta_val,0:.005:.15);
C = cos(PITCH);
S = sin(PITCH);

Vval = msubs(Vsub,[z;s;c;theta],[Z(:) S(:) C(:) THETA(:)]');
Vval = reshape(full(Vval),size(C,1),size(C,2),[]);

BIval = msubs(h_Bi,[z;s;c;theta],[Z(:) S(:) C(:) THETA(:)]');
BIval = reshape(full(BIval),size(C,1),size(C,2),[]);

figure(3)
close(3)
figure(3)
hold off
isosurface(PITCH,THETA,Z,Vval,1);
alpha(.5)
xlabel('Pitch')
ylabel('Theta')
zlabel('z')
hold on
isosurface(PITCH,THETA,Z,BIval,rho_i);
alpha(.5)


[PITCH_PHI,THETA_PHI] = meshgrid(pitch,theta_val);
z_phi = -( - (8321567036706119*cos(pitch))/9007199254740992 - (215431620425035*sin(abs(pitch)))/562949953421312 + 1040195879588265/1125899906842624);
Z_PHI= repmat(z_phi',1,length(theta_val));

surf(PITCH_PHI,THETA_PHI,Z_PHI')
alpha(.3);