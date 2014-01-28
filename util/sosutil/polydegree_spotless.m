function degree = polydegree_spotless(p,y)
%POLYDEGREE
% degree = polydegree(p,y)
%  Computes the degree of the polynomial p w.r.t the variables y
%  example:
%   p = x^2*y + z;
%   polydegree(p,[x;y;z]) = 3
%   polydegree(p,[x;z]) = 2
if isa(p,'numeric')
  degree = 0;
else
  if nargin > 1
    degree = deg(p,y);
  else
    degree = deg(p);
  end
end
end

