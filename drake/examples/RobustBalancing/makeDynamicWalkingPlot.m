clear all
% load lipmVariationWorkspace_7
load lipmVariation_switching_6
data_outer = load('V0_LIPMVariation2D');
data_outer_1 = load('V1_LIPMVariation2D');
data_samples = load('lipmVariation_switching_samples');


% zero step region
% figure(1)


dN = lipmCaptureLimit(model.T,model.foot_radius,model.step_max,model.z_nom,model.gravity,0)

r_ic = x(1) + x(2)*sqrt(model.z_nom / model.gravity);

figure(1)
colormap default

hold off

offset = .0;
scale = 1;
h=contourSpotless(scale*r_ic'*r_ic,x(1),x(2),[-1 1],[-3 3],[t; x(3:end)],zeros(model.num_states-1,1),[scale*dN^2]+offset,{'k'},500,500,offset,true);
% set(h,'Fill','On')
% set(h','Alpha',.5)
set(h,'LineWidth',2)
hold on
% caxis([-1 1])
offset = 0;
scale = 1;
h=contourSpotless(V*scale,x(1),x(4),[-1 1],[-3 3],[t;x([2;3;5;6])],zeros(model.num_states-1,1),1*scale,{'k'},500,500,0,false);
set(h,'LineWidth',2)

h=contourSpotless(data_outer.Vsol,x(1),x(4),[-1 1],[-3 3],[t;x([2;3;5;6])],zeros(model.num_states-1,1),0,{'k'},500,500,0,false);
set(h,'LineWidth',2)

% h=contour(data_samples.X,data_samples.XD,reshape(data_samples.V0,size(data_samples.X)),[1 1])

% h=contourSpotless([subs(V,x(1:2),[sqrt(model.g/model.z_nom)*x(1) - x(2);sqrt(model.g/model.z_nom)*x(1) + x(2)])],x(1),x(2),[-.2 .2],[-.5 .5],[t;x(3:end)],zeros(model.num_states-1,1),1+offset,{'k'},500,500,offset,true);
% set(h,'Fill','On')
set(gca,'LooseInset',get(gca,'TightInset'))

xlabel('x_c_m','FontSize',24)
ylabel('xdot_c_m','FontSize',24)
title('0-Step Regions','FontSize',24)
hold off
axis([-.5 .5 -1.5 1.5])
% colorbar


%%
figure(2)
hold off
dN = lipmCaptureLimit(model.T,model.foot_radius,model.step_max,model.z_nom,model.gravity,1)
h=contourSpotless(scale*r_ic'*r_ic,x(1),x(2),[-1 1],[-3 3],[t;x(3:end)],zeros(model.num_states-1,1),[scale*dN^2]+offset,{'k'},500,500,offset,true);
set(h,'LineWidth',2)

hold on
h=contourSpotless(V2,x(1),x(4),[-1 1],[-3 3],[t;x([2;3;5;6])],zeros(model.num_states-1,1),dmsubs(rho2,t,0),{'k'});
set(h,'LineWidth',2)

h=contourSpotless(data_outer_1.Vsol,x(1),x(4),[-1 1],[-3 3],[t;x([2;3;5;6])],zeros(model.num_states-1,1),0,{'k'},500,500,0,false);
set(h,'LineWidth',2)

% h=contourSpotless(V*scale,x(1),x(4),[-1 1],[-3 3],[t;x([2;3;5;6])],zeros(model.num_states-1,1),1*scale,{'k'},500,500,0,false);
% set(h,'LineWidth',2)

[cl,h]=contour(X,XD,reshape(V0,NS,NS),[1 1]*1)
set(h,'LineWidth',2)
set(h','LineColor','k')

xlabel('x_c_m','FontSize',24)
ylabel('xdot_c_m','FontSize',24)
title('1-Step Regions','FontSize',24)
axis([-.5 .5 -1.5 1.5])
% hold off

 %% simulate
% NS = 10;
% [X,XD] = meshgrid(linspace(-.5,.5,NS),linspace(-1.5,1.5,NS));
% xf = zeros(6,length(X));
% for i=1:numel(X),
%   xfi = simLIPMVariationSwitching(model,Bu,false,[X(i);0;0;XD(i);0;0]);
%   xf(:,i) = xfi;
%   i
% end
% 
% s = msspoly('s',model.num_reset_inputs);
% r = model.reset(t,x,s);
% for i=1:numel(X)
%   Vfi = subs(V,x,xf(:,i) + x - r);
%   V0(i) = double(min(dmsubs(Vfi,s,linspace(-1,1,100))));
% %   i
% end
% V0 = double(V0)
% [cl,h]=contour(X,XD,reshape(V0,NS,NS),[1 1]*1)
% set(h,'LineWidth',2)
% set(h','LineColor','r')

 %% simulate S1
NS = 20;
[X,XD] = meshgrid(linspace(-.5,.5,NS),linspace(-1.5,1.5,NS));
xf = zeros(6,length(X));
for i=1:numel(X),
  xfi = simLIPMVariationSwitching(model,B2,false,[X(i);0;0;XD(i);0;0]);
  xf(:,i) = xfi;
  i
end

s = msspoly('s',model.num_reset_inputs);
r = model.reset(t,x,s);
for i=1:numel(X)
  Vfi = subs(V,x,xf(:,i) + x - r);
  V0(i) = double(min(dmsubs(Vfi,s,linspace(-1,1,100))));
%   i
end
%%
figure
% V0 = double(V0)N
NS = sqrt(length(data_samples.V0));
[cl,h]=contour(data_samples.X,data_samples.XD,reshape(data_samples.V0,NS,NS),[1 1]*1)
set(h,'LineWidth',2)
set(h','LineColor','k')

%%
figure(3)
hold off
dN = lipmCaptureLimit(model.T,model.foot_radius,model.step_max,model.z_nom,model.gravity,0);
r_ic = x(1) + x(2)*sqrt(model.z_nom / model.gravity);
h=contourSpotless(r_ic'*r_ic,x(1),x(2),[-1 1],[-3 3],[t;x(3:end)],zeros(model.num_states-1,1),[dN^2],{'k'},500,500);
hold on
set(h,'LineWidth',2);

dN = lipmCaptureLimit(model.T,model.foot_radius,model.step_max,model.z_nom,model.gravity,1);
h=contourSpotless(r_ic'*r_ic,x(1),x(2),[-1 1],[-3 3],[t;x(3:end)],zeros(model.num_states-1,1),[dN^2],{'k'},500,500);
set(h,'LineWidth',2);

dN = lipmCaptureLimit(model.T,model.foot_radius,model.step_max,model.z_nom,model.gravity,2);
h=contourSpotless(r_ic'*r_ic,x(1),x(2),[-1 1],[-3 3],[t;x(3:end)],zeros(model.num_states-1,1),[dN^2],{'k'},500,500);
set(h,'LineWidth',2);

xlabel('x_c_m','FontSize',24)
ylabel('xdot_c_m','FontSize',24)
title('LIPM Backwards Reachable Sets','FontSize',24)
axis([-.5 .5 -1.5 1.5])

%%
close all
contour_inds = [1 4 2];
sub_inds = setdiff(1:6,contour_inds);
figure
contourSpotless3D(subs(Bu(2)*1e3,  [t;x(sub_inds)], zeros(1+length(sub_inds),1)), x(contour_inds), 0, [.5 1.5 .5]);
xlabel('x_c_m'); ylabel('xdot_c_m'); zlabel('z_c_m');

figure
h=contourSpotless(Bu(2),x(1),x(4),[-1 1],[-3 3],[t;x([2;3;5;6])],zeros(model.num_states-1,1),0,{'k'},500,500);
set(h,'LineWidth',2);
hold on
h=contourSpotless(V,x(1),x(4),[-1 1],[-3 3],[t;x([2;3;5;6])],zeros(model.num_states-1,1),1,{'k'},500,500,0,false);
set(h,'LineWidth',2)
axis([-.5 .5 -1.5 1.5])

xlabel('x_c_m','FontSize',24)
ylabel('xdot_c_m','FontSize',24)
title('Control Input: z-acceleration','FontSize',24)
