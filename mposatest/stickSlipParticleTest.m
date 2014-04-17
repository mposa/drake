%% generate data
N = 10000;
t = linspace(0,1,N);
dt = t(2) - t(1);
u0 = cos(10*t);
F_stick = .7;
F = .6;
Q = .2;
u = u0 + Q*randn(size(u0));
x = zeros(2,N);
mode = zeros(1,N);
%%
for i=2:N,
  v0 = x(2,i-1) + dt*u(i);
  
  if abs(v0) < dt*F_stick,
    v1 = 0;
    mode(i) = 0;
    t(i)
  elseif v0 > 0
    v1 = v0 - dt*F;
    mode(i) = 1;
  else
    v1 = v0 + dt*F;
    mode(i) = -1;
  end
  x(:,i) = [x(1,i-1) + dt*x(2,i-1); v1];
end

figure(1)
subplot(2,1,1)
plot(t,x)
subplot(2,1,2)
plot(t,u,t,mode)