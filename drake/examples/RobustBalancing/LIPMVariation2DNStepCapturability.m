function LIPMVariation2DNStepCapturability(n)
if nargin < 1
n = 0;
end
g = 10;
z_nom = 1;
R_diag = [1, .5, pi/2, 2, 2, 2];
% R_diag = [2, 1, 1, 4, 2, 2];
% f_min = .5;
% f_max = 1.5;
uz_bnd = .3;
ux_bnd = .1;
foot_radius = .05;
% inertia_ratio = .3^2/2; 
inertia_ratio = .6^2/2; 
step_max = .7;
step_time = .3;
z_inv_degree = 1;

model = LIPMVariation2D(g, z_nom, step_max, step_time, inertia_ratio, ux_bnd, uz_bnd, foot_radius, z_inv_degree);

options.degree = 4;
options.do_backoff = false;
options.backoff_ratio = 1.02;
options.free_final_time = false;
options.scale = 1;
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