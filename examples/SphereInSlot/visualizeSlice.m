
q = msspoly('q',2);
qd = msspoly('qd',3);
s = msspoly('s',1);
c = msspoly('c',1);

xd = qd(1);
zd = qd(2);
thetad = qd(3);

%%
Vsub = subs(Vsol,qd,zeros(3,1));
% Vsub = subs(sol.eval(Vdot_free),qd,zeros(3,1));
% Vsub = subs(sol.eval(Vdot_impact(1)),qd,zeros(3,1));
x_range = linspace(-.5,.5,1000);
z_range = linspace(0,1,1000);
[X,Z] = meshgrid(x_range,z_range);
Vval = msubs(Vsub,[x;z],[X(:) Z(:)]');
Vval = reshape(Vval,size(X,1),[]);

figure(1)
hold off
[cl, h] = contour(X,Z,Vval);
% xlabel('{\theta} (rad)','FontSize',14)
% ylabel('z (m)','FontSize',14)
clabel(cl,h);

% %%
% 
% % load(datapath('flex_mult_afix_iter_48'))
% load(datapath('no_help_sq_lowmu_flex_mult_afix_iter_17'))
% % load(datapath('sq_lowmu_flex_mult_afix_iter_14'))
% 
% % D = 0;
% % AO = eye(6);
% % AI = eye(6);
% rho_i = R
% rho_o = D
% ball_vec = [z;s;1-c;zeros(3,1)];
% 
% h_Bo = ball_vec'*AO*ball_vec;
% % h_Bo2 = .1 - (2-c-c_th) - z^2 - .05*.5*qd'*H*qd;
% 
% %searched for V with .01, worked
% h_Bi = ball_vec'*AI*ball_vec; %worked with .01 and E, but failed sdsos
% Vsub = subs(Vsol,qd,zeros(3,1));
% 
% pitch =  -.5:.01:.5;
% [PITCH,Z] = meshgrid(pitch,0:.001:.15);
% C = cos(PITCH);
% S = sin(PITCH);
% 
% Vval = msubs(Vsub,[z;s;c],[Z(:) S(:) C(:)]');
% Vval = reshape(Vval,size(C,1),[]);
% 
% BIval = msubs(h_Bi,[z;s;c],[Z(:) S(:) C(:)]');
% BIval = reshape(BIval,size(C,1),[]);
% 
% BOval = msubs(h_Bo,[z;s;c],[Z(:) S(:) C(:)]');
% BOval = reshape(BOval,size(C,1),[]);
% 
% figure(1)
% hold off
% [cl, h] = contour(PITCH,Z,Vval,[1 1],'b','Linewidth',3);
% xlabel('{\theta} (rad)','FontSize',14)
% ylabel('z (m)','FontSize',14)
% clabel(cl,h);
% hold on
% [cl, h] = contour(PITCH,Z,BIval,[rho_i rho_i]);
% clabel(cl,h);
% 
% [cl, h] = contour(PITCH,Z,BOval,[rho_o rho_o]);
% clabel(cl,h);
% 
% cpi8 = cos(pi/8);
% spi8 = sin(pi/8);
% rt2 = sqrt(2);
% z_phi = -[s*spi8 - c*cpi8 + cpi8; - s*spi8 - c*cpi8 + cpi8];
% 
% Z_PHI = max(msubs(z_phi,[s;c],[sin(pitch); cos(pitch)]));
% 
% % z_phi = -( - (8321567036706119*cos(pitch))/9007199254740992 - (215431620425035*sin(abs(pitch)))/562949953421312 + 1040195879588265/1125899906842624);
% plot(pitch,Z_PHI,'r','Linewidth',3)
