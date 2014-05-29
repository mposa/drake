% function xdot = obstest(t,x)
%   % q
%   % qdot
%   % qhat
%   % qhatdot
%   y = x(1) + .1*randn;
%   
%   xdot = zeros(4,1);
%   xdot(1) = x(2);
%   xdot(2) = cos(t) + 1*randn - .1*sign(x(2));
%   xdot(3) = x(4) + 10*(y - x(3));
%   xdot(4) = cos(t) + 10*(y - x(3));
% end

h = 1e-3;
T = 20;
t = 0:h:T;
N = length(t);

x0 = .2*randn(2,1);

x = zeros(2,N);
x(:,1) = x0;

xhat = zeros(2,N);

F = .8;
Fhat = .5*ones(1,N);
fric = zeros(1,N);

for i=2:N,
  u = cos(t(i));
  xi = x(:,i-1);
  xn = zeros(2,1);
  xn(1) = xi(1) + h*xi(2);
  if xi(2) > 0
    vnom = xi(2) + h*(u - F);
    if vnom > 0
      xn(2) = vnom;
      fric(i) = -F;
    elseif vnom < -2*h*F
      xn(2) = xi(2) + h*(u + F);
      fric(i) = F;
    else
      xn(2) = 0;
      fric(i) = -xi(2)/h - u;
    end
  elseif xi(2) < 0
    vnom = xi(2) + h*(u + F);
    if vnom < 0
      xn(2) = vnom;
      fric(i) = F;
    elseif vnom > 2*h*F
      xn(2) = xi(2) + h*(u - F);
      fric(i) = -F;
    else
      xn(2) = 0;
      fric(i) = -xi(2)/h - u;
    end
  else
    if abs(u) > F,
      xn(2) = xi(2) + h*(u - sign(u)*F);
      fric(i) = -F*sign(u);
    else
      xn(2) = 0;
      fric(i) = -xi(2)/h - u;
    end
  end
  x(:,i) = xn;
  
  K = [10;10;10];
  y = xi(1);
  xhati = xhat(:,i-1);
  xhatn = zeros(2,1);
  xhatn(1) = xhati(1) + h*xhati(2) + h*K(1)*(y - xhati(1));
  uhat = u+K(2)*(y - xhati(1));
  if xhati(2) > 0
    vnom = xhati(2) + h*(uhat - Fhat(i-1));
    if vnom > 0
      xhatn(2) = vnom;
      fric(i) = -Fhat(i-1);
    elseif vnom < -2*h*Fhat(i-1)
      xhatn(2) = xhati(2) + h*(uhat + Fhat(i-1));
      fric(i) = Fhat(i-1);
    else
      xhatn(2) = 0;
      fric(i) = -xhati(2)/h - uhat;
    end
  elseif xhati(2) < 0
    vnom = xhati(2) + h*(uhat + Fhat(i-1));
    if vnom < 0
      xhatn(2) = vnom;
      fric(i) = Fhat(i-1);
    elseif vnom > 2*h*Fhat(i-1)
      xhatn(2) = xhati(2) + h*(uhat - Fhat(i-1));
      fric(i) = -Fhat(i-1);
    else
      xhatn(2) = 0;
      fric(i) = -xhati(2)/h - uhat;
    end
  else
    if abs(uhat) > Fhat(i-1),
      xhatn(2) = xhati(2) + h*(uhat - sign(uhat)*Fhat(i-1));
      fric(i) = -Fhat(i-1)*sign(uhat);
    else
      xhatn(2) = 0;
      fric(i) = -xhati(2)/h - uhat;
    end
  end
  
%   xhatn(2) = xhati(2) + h*(u+K(2)*(y - xhati(1)));
  
  xhat(:,i) = xhatn;
  Fhat(i) = Fhat(i-1) + h*K(3)*(y - xhati(1));
% Fhat(i) = Fhat(i-1) + h*K(3)*(xi(2) - xhati(2));
end

figure(1)
plot(t,xhat,t,x)
legend('qhat','vhat','q','v')
figure(2)
plot(t,Fhat)