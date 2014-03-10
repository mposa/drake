
q = msspoly('q',3);
qd = msspoly('qd',3);
x = q(1);
z = q(2);
theta = q(3);

%%

load taylor_iter_10

% D = 0;
% AO = eye(6);
% AI = eye(6);
rho_i = R
rho_o = D
ball_vec = [z;theta;zeros(3,1)];

h_Bo = ball_vec'*AO*ball_vec;
% h_Bo2 = .1 - (2-c-c_th) - z^2 - .05*.5*qd'*H*qd;

%searched for V with .01, worked
h_Bi = ball_vec'*AI*ball_vec; %worked with .01 and E, but failed sdsos
Vsub = subs(Vsol,qd,zeros(3,1));

pitch =  -.5:.01:.5;
[THETA,Z] = meshgrid(pitch,0:.001:.15);

Vval = msubs(Vsub,[z;theta],[Z(:) THETA(:)]');
Vval = reshape(Vval,size(Z,1),[]);

BIval = msubs(h_Bi,[z;theta],[Z(:) THETA(:)]');
BIval = reshape(BIval,size(Z,1),[]);

BOval = msubs(h_Bo,[z;theta],[Z(:) THETA(:)]');
BOval = reshape(BOval,size(Z,1),[]);

figure(1)
hold off
[cl, h] = contour(THETA,Z,Vval,[1 1],'b','Linewidth',3);
xlabel('{\theta} (rad)','FontSize',14)
ylabel('z (m)','FontSize',14)
clabel(cl,h);
hold on
[cl, h] = contour(THETA,Z,BIval,[rho_i rho_i]);
clabel(cl,h);

[cl, h] = contour(THETA,Z,BOval,[rho_o rho_o]);
clabel(cl,h);

cpi8 = cos(pi/8);
spi8 = sin(pi/8);
rt2 = sqrt(2);
s = sin(pitch);
c = cos(pitch);
z_phi = -[s*spi8 - c*cpi8 + cpi8; - s*spi8 - c*cpi8 + cpi8];

Z_PHI = max(z_phi);

% z_phi = -( - (8321567036706119*cos(pitch))/9007199254740992 - (215431620425035*sin(abs(pitch)))/562949953421312 + 1040195879588265/1125899906842624);
plot(pitch,Z_PHI,'r','Linewidth',3)
