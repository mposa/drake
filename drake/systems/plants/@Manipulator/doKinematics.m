function kinsol = doKinematics(obj, q, v)
kinsol.q = q;
if nargin > 2
  kinsol.v = v;
end
kinsol.mex = false;
end

