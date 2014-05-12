
q = msspoly('q',3);
qd = msspoly('qd',3);
s_vec = msspoly('s',3);
c_vec = msspoly('c',3);
x = q(1);
z = q(2);

s = s_vec(3);

c = c_vec(3);

%%

load(datapath('flex_mult_afix_iter_48'))
load(datapath('no_help_sq_lowmu_flex_mult_afix_iter_13'))
% load(datapath('sq_lowmu_flex_mult_afix_iter_14'))

% D = 0;
% AO = eye(6);
% AI = eye(6);
rho_i = R
rho_o = D
ball_vec = [z;s;1-c;zeros(3,1)];

h_Bo = ball_vec'*AO*ball_vec;
% h_Bo2 = .1 - (2-c-c_th) - z^2 - .05*.5*qd'*H*qd;

%searched for V with .01, worked
h_Bi = ball_vec'*AI*ball_vec; %worked with .01 and E, but failed sdsos
Vsub = subs(Vsol,qd,zeros(3,1));

pitch =  -.5:.01:.5;
pitch = -.2:.01:.8;
[PITCH,Z] = meshgrid(pitch,0:.001:.13);
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
[cl, h] = contour(PITCH,Z,Vval,[1 1],'b','Linewidth',3);
xlabel('{\theta} (rad)','FontSize',24)
ylabel('z (m)','FontSize',24)
% clabel(cl,h);
hold on
[cl, h] = contour(PITCH,Z,BIval,[rho_i rho_i],'k','Linewidth',3);
% clabel(cl,h);

[cl, h] = contour(PITCH,Z,BOval,[rho_o rho_o],'m','Linewidth',3);
% clabel(cl,h);

cpi8 = cos(pi/8);
spi8 = sin(pi/8);
rt2 = sqrt(2);
z_phi = -[s*spi8 - c*cpi8 + cpi8; - s*spi8 - c*cpi8 + cpi8];

Z_PHI = max(msubs(z_phi,[s;c],[sin(pitch); cos(pitch)]));

% z_phi = -( - (8321567036706119*cos(pitch))/9007199254740992 - (215431620425035*sin(abs(pitch)))/562949953421312 + 1040195879588265/1125899906842624);
plot(pitch,Z_PHI,'r','Linewidth',3)
set(gca,'FontSize',24)