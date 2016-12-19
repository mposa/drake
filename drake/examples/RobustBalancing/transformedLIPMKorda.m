function transformedLIPMKorda()

g = 10;
z_nom = 1;
step_max = .7;
step_time = 0.3;
cop_max = .5; % set to 0 to get point foot model with no continuous inputs

model = TransformedLIPM2D(g, z_nom, step_max, step_time, cop_max);
% model = LIPM2D(g, z_nom, step_max, step_time, cop_max);

R_diag = [1 cop_max*2*sqrt(g/z_nom)];
% R_diag = [1 1];

options.split_inputs = false;
options.beta = 1;
options.R_diag = R_diag;
options.ubar = 1;
options.l_x = @(x) x(1)^2 + x(2)^2;
options.l_u = @(x) 1e-1;
options.degree = 10;
options.M = 1.01 * (options.l_x([1; 0]) + options.l_u([1; 0]) * options.ubar) / options.beta;
sol = korda2015(model, options);

figure(1);
clf;
R_diag = options.R_diag;
xlim = [-R_diag(1) R_diag(1)];
ylim = [-R_diag(2) R_diag(2)];
axis([xlim, ylim]);
hold on;

while true
  x0 = zeros(model.num_states, 1);
  [x0(1), x0(2)] = ginput(1);
  lastwarn('');
  T = 10;
  [~, x_traj] = ode45(@(t, x) sol.fbar(x), [0 T], x0);
  if strcmp(lastwarn, '');
    plot(x_traj(:, 1), x_traj(:, 2), 'k-');
    axis([xlim, ylim]);
  end
end

end
