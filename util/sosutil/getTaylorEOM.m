function [f_free,f_impact,phi,phidot,psi,mss_vars] = getTaylorEOM(p,degree,x0,u0)
  nq = p.getNumPositions;
  nv = p.getNumVelocities;
  nu = p.getNumInputs;
  nC = p.getNumContactPairs;
  
  q_mss = msspoly('q',nq);
  v_mss = msspoly('v',nq);
  u_mss = msspoly('u',nu);
  lx_mss = msspoly('lx',nC);
  
  mss_vars = [q_mss;v_mss;u_mss;lx_mss];
  
  vars_init = zeros(nq+nv+nu+nC,1);
  if nargin > 2
    vars_init(1:nq+nv) = x0;
  end
  if nargin > 3
    vars_init(nq+nv+1:nq+nv+nu) = u0;
  end
  
  tvars = TaylorVar.init(vars_init,degree);
  q=tvars(1:nq);
  v=tvars(nq+1:nq+nv);
  u=tvars(nq+nv+1:nq+nv+nu);
  lx=tvars(nq+nv+nu+1:nq+nv+nu+nC);
  lz = ones(nC,1);
  
  [H,C,B] = p.manipulatorDynamics(q,v);
  
  invH = inv(H);
  
  f_free = invH*(B*u - C);
  
  [phi,normal,d,xA,xB,idxA,idxB,mu,n,D] = p.contactConstraints(q);
  phidot = n*v;
  psi = D{1}*v;
  
  f_impact = invH*(n'*diag(lz) + D{1}'*diag(lx));
end