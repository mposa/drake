function [H,C,B,phi,phidot,psi,J,J_f,K,S,U] = torsoEOM_mss(q,v,s_vec,c_vec)
p = PlanarRigidBodyManipulator('TorsoBalance.urdf');

x = q(1);
z =  q(2);

xd = v(1);
zd = v(2);
pitchd = v(3);
thetad = v(4);

s = s_vec(3);
c = c_vec(3);

s_th = s_vec(4);
c_th = c_vec(4);

q_trig=TrigPoly(q,s_vec,c_vec);

[H,C,B] = p.manipulatorDynamics(q_trig,v);
H = getmsspoly(H);
C = getmsspoly(C);

%crude cleanup of H--NEED TO CHECK FOR S^2 IN THE COEFFICIENTS RET VAL
for i=1:4,
  for j=1:4,
    H(i,j) = trimTrigPoly(H(i,j),c_th,s_th);
    H(i,j) = trimTrigPoly(H(i,j),c,s);   
  end
  C(i) = trimTrigPoly(C(i),c_th,s_th);
  C(i) = trimTrigPoly(C(i),c,s);
end

function f = trimTrigPoly(f,s_var,c_var)
%   [Rs,ps] = pdecomp(f,s_var)
%   [Rc,pc] = pdecomp(f,c_var);
%   
%   s_ind = find(ps == 2);
%   c_ind = find(pc == 2);
%   
%   if ~isempty(s_ind) && ~isempty(c_ind)
%     
%   end
  [vars,pows,coeffs]=decomp(f);
  c_ind = find(msseq(vars,c_var));
  s_ind = find(msseq(vars,s_var));
  
  if ~isempty(s_ind) && ~isempty(c_ind)
    %terms quadratic in s
    s_quad_ind = find(pows(:,s_ind) >= 2);
    s_quad_coeff = coeffs(s_quad_ind);
    s_quad_pows = pows(s_quad_ind,:);
    s_quad_pows(:,s_ind) = s_quad_pows(:,s_ind) - 2;
    
    s_mat = [s_quad_pows s_quad_coeff'];
    
    c_quad_ind = find(pows(:,c_ind) >= 2);
    c_quad_coeff = coeffs(c_quad_ind);
    c_quad_pows = pows(c_quad_ind,:);
    c_quad_pows(:,c_ind) = c_quad_pows(:,c_ind) - 2;
    
    c_mat = [c_quad_pows c_quad_coeff'];
    
    [match,indx] = ismember(s_mat,c_mat,'rows');
    
    inds = find(match,1);
    
    if isempty(inds)
      return
    end
    
    elems = vars'.^s_mat(inds,1:end-1);
    m = s_mat(inds,end);
    for loopv=1:length(elems),
      m = m*elems(loopv);
    end
    f = f + m - m*s_var^2 - m*c_var^2;
    f = trimTrigPoly(f,s_var,c_var);
  end
end


phi = [z - (8321567036706119*c)/9007199254740992 + (215431620425035*s)/562949953421312 + 1040195879588265/1125899906842624;
       z - (2^(1/2)*((8321567036706119*c)/9007199254740992 - (215431620425035*s)/562949953421312))/2 - (2^(1/2)*((215431620425035*c)/562949953421312 + (8321567036706119*s)/9007199254740992))/2 + 1040195879588265/1125899906842624];

phidot =  [zd - (8321567036706119*-s*pitchd)/9007199254740992 + (215431620425035*c*pitchd)/562949953421312;
           zd - (2^(1/2)*((8321567036706119*-s*pitchd)/9007199254740992 - (215431620425035*c*pitchd)/562949953421312))/2 - (2^(1/2)*((215431620425035*-s*pitchd)/562949953421312 + (8321567036706119*c*pitchd)/9007199254740992))/2];

psi = [xd + (8321567036706119*pitchd*c)/9007199254740992 - (215431620425035*pitchd*s)/562949953421312;
       xd + (2^(1/2)*((8321567036706119*pitchd*c)/9007199254740992 - (215431620425035*pitchd*s)/562949953421312))/2 + (2^(1/2)*((215431620425035*pitchd*c)/562949953421312 + (8321567036706119*pitchd*s)/9007199254740992))/2];

phi_perp = [x + (8321567036706119*s)/9007199254740992 + (215431620425035*c)/562949953421312;
            x + (8321567036706119*s)/9007199254740992 - (215431620425035*c)/562949953421312];

% J = [0 1 (8321567036706119*s)/9007199254740992 -(215431620425035*c)/562949953421312 0;
%   0 1 -(8321567036706119*s)/9007199254740992 -(215431620425035*c)/562949953421312 0];
% 
% J_f = [1 0 ((8321567036706119*c)/9007199254740992 + (215431620425035*s)/562949953421312) 0;
%   1 0 ((8321567036706119*c)/9007199254740992 - (215431620425035*s)/562949953421312) 0];

J = [0 1 (8321567036706119*s)/9007199254740992 + (215431620425035*c)/562949953421312 0;
  0 1 (8321567036706119*s)/9007199254740992 - (215431620425035*c)/562949953421312 0];

J_f = [1 0 ((8321567036706119*c)/9007199254740992 - (215431620425035*s)/562949953421312) 0;
  1 0 ((8321567036706119*c)/9007199254740992 + (215431620425035*s)/562949953421312) 0];


% J_f = J_f([2;1],:);


%potential energy
doKinematics(p.model,q_trig);
g = -p.model.gravity(2);
U = 0;
for i=1:length(p.model.body),
  T=p.model.body(i).T;
  I = p.model.body(i).I;
  if I(3,3) ~= 0
    tmp = T*[I(1,3)/I(3,3);-I(1,2)/I(3,3);1];
    U = U + tmp(2)*g*I(3,3);
  end
end
U = getmsspoly(U);
U = U - subs(U,[q;s_vec;c_vec],[zeros(size(q));zeros(size(s_vec));ones(size(c_vec))]);


% Generate a linear model
[A,B_lin] = p.linearize(0,zeros(8,1),0);
A_sub = A([4;8],[4;8]);
B_sub = B_lin([4;8]);
Q = eye(2);
R = 1;
[K_sub,S_sub] = lqr(A_sub,B_sub,Q,R);
K = K_sub;
S = S_sub;
% K = [zeros(1,3) K_sub(1) zeros(1,3) K_sub(2)];
% S = zeros(8);
% S([4;8],[4;8]) = S_sub;
end




