function [H,Hf_free,Hf_impact,phi,phidot,psi] = getTrigEOM(p)
  prog = spotsosprog;
  nq = p.getNumPositions;
  nv = p.getNumVelocities;
  nu = p.getNumInputs;
  nC = p.getNumContactPairs;
  
  q = msspoly('q',nq);
  s = msspoly('s',nq);
  c = msspoly('c',nq);
  v = msspoly('v',nv);
  u = msspoly('u',nu);
  lx = msspoly('lx',nC);
  lz = ones(nC,1);
  
  q_trig = TrigPoly(q,s,c);
  
  [H,C,B] = p.manipulatorDynamics(q_trig,v); 
  
  Hf_free = clean(getmsspoly(B*u - C));
  
  [phi,normal,d,xA,xB,idxA,idxB,mu,n,D] = p.contactConstraints(q_trig);
  phi = clean(getmsspoly(phi));
  phidot = clean(getmsspoly(n*v));
  psi = clean(getmsspoly(D{1}*v));
  
  Hf_impact = clean(getmsspoly(n'*diag(lz) + D{1}'*diag(lx)));
  
  H = clean(getmsspoly(H));
  
  H = prog.trigExprReduction(H,s,c);
  Hf_free = prog.trigExprReduction(Hf_free,s,c);
  Hf_impact = prog.trigExprReduction(Hf_impact,s,c);
  phi = prog.trigExprReduction(phi,s,c);
  phidot = prog.trigExprReduction(phidot,s,c);
  psi = prog.trigExprReduction(psi,s,c);
end