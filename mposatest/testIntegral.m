% A = diag(rand(2,1))
%
% th = linspace(0,2*pi,1000);
% f = exp(-A(1,1)*sin(th).^2 - A(2,2)*cos(th).^2);
%
% I0 = mean(f)*2*pi
%
% C = exp(-1/2*(A(1,1)+A(2,2)));
%
% I1 = 2*pi*C*besseli(0,1/2*(A(1,1)-A(2,2)))
%
% I2 = C*2*pi*mean(exp(1/2*(A(1,1)-A(2,2))*cos(2*th)))
%
%
% F0 = mean(exp(-A(2,2)*cos(th)))*2*pi
%
% F1 = 2*pi*besseli(0,-A(2,2))

A = randn(2,2);
A = (A+A')/2;
while min(eig(A)) < 0
  A = randn(2,2);
  A = (A+A')/2;
end

x = randn(2,1);
x = x/norm(x);

1 - .5*x'*A*x

exp(-x'*A*x)