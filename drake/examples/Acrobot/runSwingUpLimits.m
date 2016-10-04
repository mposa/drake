% run optimization once with slack variable set to 1
[p,xtraj,utraj,ljltraj,info,z,traj_opt] = swingUpLimits(1);

% run optimization once with slack variable set to 0.1
[p,xtraj,utraj,ljltraj,info,z,traj_opt] = swingUpLimits(0.1,xtraj,utraj,ljltraj);

% run again with no slack
[p,xtraj,utraj,ljltraj,info,z,traj_opt] = swingUpLimits(0,xtraj,utraj,ljltraj);

% playback the trajectory
v = p.constructVisualizer;
v.axis = [-3 3 -4 4];
v.playback(xtraj);

%%
figure(1)
plot(z(traj_opt.ljl_inds)')
xlabel('time')
ylabel('Joint limit forces')

%%
figure(2)
plot(z(traj_opt.x_inds(2,:))')
xlabel('time')
ylabel('elbow position (rad)')