function lipmHeightVariationCapturability(n)
if nargin < 1
n = 0;
end
g = 10;
z_nom = 1;
uz_bnd = .5;
foot_radius = .05;
step_max = .7;
step_time = .3;
z_inv_degree = 3;

model = LIPMHeightVariation2D(g, z_nom, step_max, step_time, uz_bnd, foot_radius, z_inv_degree);


R_diag = [.5, .5, 1, 2];

options.degree = 6 ;
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