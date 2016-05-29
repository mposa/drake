function variableHeightAndPitchNStepCapturability(n)
if nargin < 1
n = 0;
end
g = 1;
z_nom = 1;
R_diag = [2, 1, 1, 2, 2, 2]*2;
% R_diag = [2, 1, 1, 4, 2, 2];
% f_min = .5;
% f_max = 1.5;
f_min = .9;
f_max = 1.1;
div_max = .1;
% inertia_ratio = .3^2/2; 
inertia_ratio = .6^2/2; 
step_max = .7;
step_time = .3;

model = VariableHeightandPitch2D(g, z_nom, step_max, step_time, inertia_ratio, f_max, f_min, div_max);

options.degree = 4;
options.do_backoff = false;
options.backoff_ratio = 1.02;
options.free_final_time = false;
options.scale = 1/3;
options.scale_input = 1;
options.control_design = false;
options.korda_control_design = true;

target = [];
% goal_radius = 0.01;
% target = @(x) goal_radius^2 - x'*x;

if n > 0
  T = step_time;
else
  T = 1;
end

nStepCapturabilitySOS(model, T, R_diag, target, n, options)

end