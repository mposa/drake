function [ Y,y,Z,sigma,R,Rt] = mprw2( A, b, C, max_iter, sigma, tol, Y, Z, R, Rt)
% simple "regularization method" : 
% (based on primal prox, or dual augmented Lagrangien)
% corresponds to boundary point method, and to block-coordinate descent

% solves: min <C,X> s.t. A*X(:) = b; X psd
%         max b'y s.t. C - A^t(y) = Z psd
% max_iter: max number of iterations
% sigma: penalty parameter for augmented Lagrangian 
%        a typical value (b and C normalized) is between .1 and 10
% tol: stopping condition
%      we stop once primal and dual relative infeasibility < tol
% 
% simplest call: [ X, y, Z] = mprw( A, b, C);
% call: [ X,y,Z,sigma,R,Rt] = mprw( A, b, C, max_iter, sigma, tol, X, Z,R, Rt);
% optional input (for restart):  sigma, tol, X, Z

% version 2 : 
% * scaling (it should help polynomial optimization)
% * slight change of sigma update

tstart = cputime; 
m=min(size(A));
n2=max(size(A));
n = sqrt(n2);
if length(b)~=m; error('Size A mismatches b.'); end
if (min(size(C)) == 1)&&(max(size(C)) == n2); C=reshape(C,n,n); end
if (size(C,1)~=n)||(size(C,2)~=n); error('Size C mismatches A.'); end

% rescale data
normA0 = full(min(1e12,max(1,norm(A,'fro'))));%full(min(1e12,max(1,sqrt(A(:)'*A(:)))));
normb0 = max(1,norm(b));
normC0 = full(min(1e12,max(1,sqrt(C(:)'*C(:)))));
A = A/normA0;
C = C/normC0; 
b = b/normb0; 
%R = R/normA0;
%R = R/sqrt(normA0);
%Rt = Rt/normA0;
%Rt = Rt/sqrt(normA0); 

% form A*A' (and cholesky, if not part of input) 
if size(A,1) == m
  At = A';
else
  At = A; A = At';
end
AAT = A* At;           % form A*A'
if nargin <= 6;
   [R,p] = chol( AAT);
   Rt = R';    % only for speed-up
   secs = cputime - tstart;
   fprintf(' secs after chol:   %12.5f \n', secs);
   if p>0; 
         fprintf(' rows of A linearly dependent. p = %4.0d, m = %4.0d\n',p,m);
          error(' nothing done');
   end
end


% sinon
normb = norm(b);
normC = sqrt(C(:)'*C(:));

% initialize
if nargin == 3; max_iter=100; sigma=1; tol=1e-5; end
if nargin == 4; sigma = 1; tol = 1e-5; end
if nargin == 5; tol = 1e-5; end
if nargin <= 6; Y = zeros(n); Z = Y; end;
  
% outer stopping condition
iter = 1;              % iteration count
g = b-A*Y(:);          % primal residue
done = 0;              % stopping condition 
update_it = 10;        % update sigma every update_it iterations 
fprintf(' it     secs       dual         primal     lg(rrd)   lg(rrp)   sigma\n');

% start outer iterations
while done == 0;       % while stopping conditions not satisfies

% given Y,Z  and sigma, solve for y
  rhs = A*(C(:) - Z(:)) + g/sigma;
  y_tmp = Rt\rhs; y = R\y_tmp;
%  y = AAT\rhs;         % solve sigma*AA^Ty = sigma*A(L+Z) + A(Y)-b

% now form W= W(y,Y,sigma)  
  Aty = reshape(At*y, n, n);    % A^T(y)
  M = Aty - C + Y/sigma;
  M = (M + M')/2;      % should be symmetric
  
% now compute eigenvalue decomposition of W  
  [ev, lam] = eig(M); lam = diag(lam);

% compute projection Mp and Mn of M onto PSD and NSD
  I = find(lam> 0); j = length(I);
  if j < n/2;
     evp = zeros( n,j);
     for r=1:j;
         ic = I(r); evp(:,r) = ev(:,ic)*sqrt(lam(ic)); 
     end
     if j==0; evp = zeros(n,1); end
     evpt= evp';
     Mp = evp*evpt;     
     Mn = M - Mp;
  else
     I = find(lam < 0); j = length( I); % should be <= n/2
     evn = zeros( n,j);
     for r=1:j;
        ic = I(r); evn(:,r) = ev(:,ic)*sqrt( -lam(ic)); 
     end
     if j==0; evn = zeros(n,1); end
     evnt= evn';
     Mn = -evn*evnt;
     Mp = M - Mn; 
  end

  Mn = (Mn+Mn')/2; Mp = (Mp+Mp')/2;  % these should be symm.

  Z = -Mn; X = sigma * Mp;
  g = b-A * X(:);

% update Y
  Y = X;
  G = -Z + C - Aty;

% some output
err_d = norm( G,'fro'); dual = b'*y*normC0*normb0/normA0; 
secs = cputime-tstart;
err_p = norm(g);   primal = C(:)'*Y(:)*normC0*normb0/normA0;
rel_err_p = err_p/(1+normb); rel_err_d = err_d/(1+normC);
iter = iter+ 1;

if (mod(iter,50)==0);
fprintf( '%3.0d %8.2f %13.5e %13.5e %8.3f %8.3f  %9.6f\n', iter, ...
   secs,  dual, primal,log10(rel_err_d), log10(rel_err_p), sigma );
end

% check stopping conditions
if max(rel_err_d, rel_err_p) < tol || iter > max_iter ; 
  fprintf( '%3.0d %8.2f %13.5e %13.5e %8.3f %8.3f  %9.6f\n', iter, ...
  secs, dual, primal, log10(rel_err_d), log10(rel_err_p), sigma );
  fprintf('total time: %10.3f \n', secs); 
  if (iter>max_iter);  
      fprintf('max outer iterations reached. \n');
  end
  done=1;
end

% check for reduction of sigma
ratio = rel_err_p/rel_err_d;  
const=1.2; % 1.2
if (ratio < 0.2) && (mod(iter,update_it)==0) %0.2
   sigma = min(1e3,const*sigma);      
elseif (ratio > 5)  % 5
   sigma = max(1e-2,sigma/const); 
end

end

%descale data
Y=Y*normb0/normA0;
y=y*normC0/normA0;
Z=normC0*Z;
