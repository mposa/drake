function variableHeightPointMass2DNStepCapturability(n)

g = 10;
step_max = .7;
step_time = 0.3;
f_max = 1.5 * g;
z_nom = 1;

model = VariableHeightPointMass2D(g, z_nom, step_max, step_time, f_max);

if n > 0
  T = step_time;
else
  T = 2;
end
options.degree = 6;

% R_diag = 2 * ones(1, model.num_states);
R_diag = [1, 0.5, 1, 1];

% goal_radius = 0.1;
% target = @(x) goal_radius^2 - x'*x;
target = zeros(model.num_states, 1);

nStepCapturabilitySOS(model, T, R_diag, target, n, options)

end
