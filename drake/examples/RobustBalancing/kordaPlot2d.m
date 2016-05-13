function kordaPlot2d(model,plant, x,g_X,g_X_target, V_outer, V_inner, R_diag)
clf;
legend_strings = {};
xlim = [-R_diag(1) R_diag(1)];
ylim = [-R_diag(2) R_diag(2)];
hold on;

handles = [];

handles = [handles; contourSpotless(g_X, x(1), x(2), xlim, ylim, [], [], 0, {'k'})];
legend_strings{end + 1} = 'g_X';

handles = [handles; contourSpotless(g_X_target, x(1), x(2), xlim, ylim, [], [], 0, {'g'})];
legend_strings{end + 1} = 'g_X^T';

handles = [handles; contourSpotless(V_outer, x(1), x(2), xlim, ylim, [], [], 0, {'b'})];
legend_strings{end + 1} = 'v (outer)';

handles = [handles; contourSpotless(V_inner, x(1), x(2), xlim, ylim, [], [], 0, {'r'})];
legend_strings{end + 1} = 'v (inner)';
% doesn't work because of Matlab bug:
%       handles = [handles; inner_contour_handles];
%       for i = 1 : nbeta
%         legend_strings{end + 1} = ['v (inner), \beta = ' num2str(betas(i))];
%       end

xlabel('x_1')
ylabel('x_2')
legend(gca(), handles, legend_strings);

% set(gcf, 'ButtonDownFcn',@ImageClickCallback)
% 
%   function ImageClickCallback(objectHandle, eventData)
%     coordinates = get(axesHandle,'CurrentPoint');
%     coordinates = coordinates(1,1:2);
%     message     = sprintf('x: %.1f , y: %.1f',coordinates (1) ,coordinates (2));
%     helpdlg(message);
%   end

while true
  x0 = zeros(model.num_states, 1);
  mouse = ginput(1);
  if isempty(mouse)
    break
  else
    x0 = mouse';
  end
  lastwarn('');
  
  try
    traj = plant.simulate([0 10],[1;x0;0;0;0]);
    x_traj = traj.eval(traj.pp.breaks)';
    if strcmp(lastwarn, '');
      plot(x_traj(:, 1), x_traj(:, 2), 'k-');
    end
  catch e
    display(e.message)
  end
%   [~, x_traj] = ode45(@(t, x) msubs(sol.fbar, x, x), [0 10], x0);

end

hold off;
end