
load V0_TransformedFull2DModel_inner_2.mat
x = msspoly('x',model.num_states);
t = msspoly('t',1);
r_ic = x(1) + x(2)*sqrt(model.z_nom / model.g);
dN=lipmCaptureLimit(model.T,model.foot_radius, model.step_max, model.z_nom, model.g, 0)
figure(1)
colormap winter
hold off

offset = 0;
h=contourSpotless([subs(V_inner,x(1:2),[sqrt(model.g/model.z_nom)*x(1) - x(2);sqrt(model.g/model.z_nom)*x(1) + x(2)])],x(1),x(2),[-.2 .2],[-.5 .5],[t;x(3:end)],zeros(model.num_states-1,1),1+offset,{'k'},500,500,offset,true);
set(h,'Fill','On')
hold on

scale = 1
h=contourSpotless(scale*r_ic'*r_ic,x(1),x(2),[-.2 .2],[-.5 .5],[t;x(3:end)],zeros(model.num_states-1,1),[scale*dN^2],{'k'},500,500,0,false);
% set(h,'Fill','On')
% % set(h','Alpha',.5)
  set(h,'LineWidth',5)

  set(gca,'LooseInset',get(gca,'TightInset'))

xlabel('x_c_m','FontSize',24)
ylabel('xdot_c_m','FontSize',24)
  title('0-Step Regions','FontSize',24)
hold off