function scaledFullNStepCapturability(n)
if nargin < 1
n = 0;
end
g = 1;
z_nom = 1;
R_diag = [2, 1, 1, 2, 2, 2];
% f_min = .5;
% f_max = 1.5;
fx_range = 1;
fz_range = .5;
% inertia_ratio = .3^2/2; 
inertia_ratio = .6^2/2; 
step_max = .7;
step_time = .3;
foot_size = .4;

model = ScaledFull2D(g, z_nom, step_max, step_time,  inertia_ratio, fx_range, fz_range, foot_size);

options.degree = 4;
options.do_backoff = false;
options.backoff_ratio = 1.02;
options.free_final_time = false;
options.scale = 1;
options.scale_input = 1;
options.control_design = false;

target = [];
% goal_radius = 0.01;
% target = @(x) goal_radius^2 - x'*x;

if n > 0
  T = step_time;
else
  T = 5;
end

nStepCapturabilitySOS(model, T, R_diag, target, n, options)

end