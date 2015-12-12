function simVariableHeightPointMass2D
model_name = 'VariableHeightPointMass2D';
n = 1;
sample_slice = true;

V_data = load(solutionFileName(model_name,n));
V0_data = load(solutionFileName(model_name,n-1));

model = V_data.model;
step_time = V_data.T;
Vsol = V_data.Vsol;

t = msspoly('t',1);
x = msspoly('x',model.num_states);
u = msspoly('u',model.num_inputs);
s = msspoly('s',model.num_reset_inputs);
xdot = model.dynamics(t,x,u);

Vdot = diff(Vsol,x)*xdot + diff(Vsol,t);
dVdotdu = diff(Vdot,u);

[pows,coeffs] = decomp_ordered(dVdotdu,[t;x]);

plant = NStepCapturabilityPlant(model);
controller = NStepCapturabilityController(plant,coeffs,pows);


sys_cl = plant.feedback(controller);


xp = model.reset(t, x, s);

V0p = subs(V0_data.Vsol,[x;t],[xp;0]);


sys_cl = sys_cl.setSimulinkParam('Solver','ode4');
sys_cl = sys_cl.setSimulinkParam('FixedStep','.01');

%%
if sample_slice
  sub_inds = [2;4];
  sub_val = 0*sub_inds;
  sample_inds = [1;3];
  N_sample = 100;
  [XS_1,XS_2] = meshgrid(linspace(-V_data.R_diag(sample_inds(1)),V_data.R_diag(sample_inds(1)),N_sample),...
    linspace(-V_data.R_diag(sample_inds(2)),V_data.R_diag(sample_inds(2)),N_sample));
  XS_1 = XS_1(:);
  XS_2 = XS_2(:);
  VOPT = -inf(size(XS_1));
  for i=1:length(XS_1),
    x0 = zeros(model.num_states,1);
    x0(sample_inds) = [XS_1(i);XS_2(i)];
    if double(subs(Vsol,[t;x],[0;x0])) > 0 && norm(x0./V_data.R_diag') < 1
      % do simulation
      traj = sys_cl.simulate([0 step_time],x0);
      
      t_sim = traj.pp.breaks;
      x_sim = traj.eval(t_sim);
      
      V0p_sim = subs(V0p,x,traj.eval(step_time));
      S = linspace(-model.step_max,model.step_max,100);
      V0p_vals = msubs(V0p_sim,s,S);
      [V0_opt,i_opt] = max(V0p_vals);
      V0_opt
      VOPT(i) = V0_opt;
      s_opt = S(i_opt);
    end
  end
  Ip = find(VOPT >= 0);
  In = (VOPT < 0 & VOPT > -inf);
  figure(n*10+2)
  hold on
  plot(XS_1(Ip),XS_2(Ip),'go',XS_1(In),XS_2(In),'rx')
  hold off
else
x0 = [0;0;.9;0];


t_sim = traj.pp.breaks;
x_sim = traj.eval(t_sim);
V_sim = msubs(Vsol,[t;x],[t_sim;x_sim]);
figure(1)
subplot(2,1,1)
plot(t_sim,x_sim)
legend('x_1','x_2','x_3','x_4')


V0p_sim = subs(V0p,x,x_sim(:,end));
S = linspace(-step_max,step_max,1000);
V0p_vals = msubs(V0p_sim,s,S);
[V0_opt,i_opt] = max(V0p_vals);
V0_opt
s_opt = S(i_opt)

end

%% 
% x0_0step = double(subs(xp,[x;s],[x_sim(:,end);s_opt]));
% % x0_0step = [0-.5;0;1.623;0];
% Vdot_0step = diff(V0.Vsol,x)*xdot + diff(V0.Vsol,t);
% dVdotdu_0step = diff(Vdot_0step,u);
% 
% [pows_0step,coeffs_0step] = decomp_ordered(dVdotdu_0step,[t;x]);
% controller_0step =  NStepCapturabilityController(plant,coeffs_0step,pows_0step);
% 
% sys_cl_0step = plant.feedback(controller_0step);
% sys_cl_0step = sys_cl_0step.setSimulinkParam('Solver','ode4');
% sys_cl_0step = sys_cl_0step.setSimulinkParam('FixedStep','.001');
% traj_0step = sys_cl_0step.simulate([0 1],x0_0step);
% 
% t_sim_0step = traj_0step.pp.breaks;
% x_sim_0step = traj_0step.eval(t_sim_0step);
% 
% V0_sim_0step = msubs(V0.Vsol,[t;x],[t_sim_0step;x_sim_0step]);
% 
% figure(1)
% subplot(2,1,2)
% plot(t_sim_0step,x_sim_0step)
% legend('x_1','x_2','x_3','x_4')
end


function filename = solutionFileName(model_name, n)
filename = sprintf(['V%d_' model_name '.mat'], n);
end