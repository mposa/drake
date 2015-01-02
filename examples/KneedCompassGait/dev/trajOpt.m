
%%
tic
[p,xtraj,utraj,ltraj,ljltraj,z,F,info,traj_opt] = testNewTrajOpt();
toc
v = p.constructVisualizer; v.playback(xtraj);
figure(1)
t = xtraj.pp.breaks;
l = ltraj.eval(t);
plot(t,l(1:4:end,:))
display('iter 1')
%%
tic
[p,xtraj2,utraj2,ltraj2,ljltraj2,z2,F2,info2,traj_opt2] = testNewTrajOpt(xtraj,utraj,ltraj,ljltraj,.1);
toc
v = p.constructVisualizer; v.playback(xtraj2);
figure(1)
t = xtraj2.pp.breaks;
l = ltraj2.eval(t);
plot(t,l(1:4:end,:))
display('iter 2')
%%
tic
[p,xtraj3,utraj3,ltraj3,ljltraj3,z3,F3,info3,traj_opt3] = testNewTrajOpt(xtraj2,utraj2,ltraj2,ljltraj2,.01);
toc
v = p.constructVisualizer; v.playback(xtraj3);
figure(1)
t = xtraj3.pp.breaks;
l = ltraj3.eval(t);
plot(t,l(1:4:end,:))
display('iter 3')
%%
tic
[p,xtraj4,utraj4,ltraj4,ljltraj4,z4,F4,info4,traj_opt4] = testNewTrajOpt(xtraj3,utraj3,ltraj3,ljltraj3,1e-3);
toc
v = p.constructVisualizer; v.playback(xtraj4);
figure(1)
t = xtraj4.pp.breaks;
l = ltraj4.eval(t);
plot(t,l(1:4:end,:))
display('iter 4')
%%
tic
[p,xtraj5,utraj5,ltraj5,ljltraj5,z5,F5,info5,traj_opt5] = testNewTrajOpt(xtraj4,utraj4,ltraj4,ljltraj4,0);
toc
figure(1)
t = xtraj5.pp.breaks;
l = ltraj5.eval(t);
plot(t,l(1:4:end,:))
v = p.constructVisualizer; v.playback(xtraj5);
display('iter 5')
%%
if info5 ~= 1 || true
  tic
[p,xtraj5,utraj5,ltraj5,ljltraj5,z5,F5,info5,traj_opt5] = testNewTrajOpt(xtraj5,utraj5,ltraj5,ljltraj5,0);
toc
v = p.constructVisualizer; v.playback(xtraj5);
display('iter 5.2')
end

