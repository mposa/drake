[r_inv,reset_constraint_inv] = model_swing.inverse_reset(t,x,s);
% reset_constraint_inv = s;
V3 = x'*Q*x*100;
rho3 = msspoly(1) + 10*t;

%%
for i=1:10,
[V3,rho3] = switchingControlHybridApproxTest(x,s,V3,rho3,step_time,r_inv,reset_constraint_inv,V);
% [V3,rho3] = switchingControlHybridApproxTestNormalReset(x,s,V3,rho3,step_time,r,reset_constraint,V);  % achieves a lower cost
rho3 = msspoly(rho3);
figure(2)
  hold off
  contourSpotless(V,x(1),x(3),[-1 1],[-3 3],[t;x([2;4;5])],zeros(model_swing.num_states-1,1),1,{'g'});
  hold on
  contourSpotless(V3,x(1),x(3),[-1 1],[-3 3],[t;x([2;4;5])],[zeros(model_swing.num_states-2,1);0],dmsubs(rho3,t,0),{'k'});
  contourSpotless(V3,x(1),x(3),[-1 1],[-3 3],[t;x([2;4;5])],[step_time;zeros(model_swing.num_states-3,1);0],dmsubs(rho3,t,step_time),{'r'});
  hold off
  
  figure(3)
  hold off
  contourSpotless(V,x(1),x(3),[-1 1],[-3 3],[t;x([2;4;5])],zeros(model_swing.num_states-1,1),1,{'g'});
  hold on
  contourSpotless(V3,x(1),x(3),[-1 1],[-3 3],[t;x([2;4;5])],[zeros(model_swing.num_states-2,1);.7],dmsubs(rho3,t,0),{'k'});
  contourSpotless(V3,x(1),x(3),[-1 1],[-3 3],[t;x([2;4;5])],[step_time;zeros(model_swing.num_states-3,1);.7],dmsubs(rho3,t,step_time),{'r'});
  hold off
end