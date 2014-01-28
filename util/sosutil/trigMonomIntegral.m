function I = trigMonomIntegral(alpha, xmin, xmax)
% function I = trigMonomIntegral(alpha, xmin, xmax)
%   for f = cos(theta)^alpha(1) * sin(theta)^alpha(2)
%     compute the definite integral I = int(f,xmin,xmax)
%   where alpha must be 1x2

  % check size
  assert(isequal(size(alpha),[1 2]));
  
  % check non-negativity
  assert(all(alpha >= 0))
  
  % check integer
  assert(all(rem(alpha,1) == 0))
  
  % constant term
  if ~any(alpha)
    I = xmax - xmin;
    return
  end
  
  % sin only
  if alpha(1) == 0
    if alpha(2) == 1
      %sin
      I = -cos(xmax) + cos(xmin);
    else
      %sin^n
      n = alpha(2);
      I = -sin(xmax)^(n-1)*cos(xmax)/n + sin(xmin)^(n-1)*cos(xmin)/n + (n-1)/n*trigMonomIntegral([0, n-2], xmin, xmax);
    end
    return
  end
  
  % cos only
  if alpha(2) == 0
    if alpha(1) == 1
      %cos
      I = sin(xmax) - sin(xmin);
    else
      n = alpha(1);
      I = cos(xmax)^(n-1)*sin(xmax)/n - cos(xmin)^(n-1)*sin(xmin)/n + (n-1)/n*trigMonomIntegral([n-2, 0], xmin, xmax);
    end
    return
  end
  
  % mixed sin and cos
  n = alpha(2);
  m = alpha(1);
  if m > 1
    I = sin(xmax)^(n+1)*cos(xmax)^(m-1)/(n+m) - sin(xmin)^(n+1)*cos(xmin)^(m-1)/(n+m) + (m-1)/(n+m)*trigMonomIntegral([m-2,n], xmin, xmax);
  elseif n > 1
    I = -sin(xmax)^(n-1)*cos(xmax)^(m+1)/(n+m) + sin(xmin)^(n-1)*cos(xmin)^(m+1)/(n+m) + (n-1)/(n+m)*trigMonomIntegral([m,n-2], xmin, xmax);
  else
    %sin*cos
    I = -.5*cos(xmax)^2 + .5*cos(xmin)^2;
  end
end