[X,Y] = meshgrid(linspace(-5,5,200),linspace(-4,4,200));
XMAT = [X(:)';Y(:)';zeros(4,numel(X))];
s = linspace(-1,1,100);
L = chol(Q);
vecs=kron(XMAT,ones(1,numel(s))) + repmat(c*s,1,numel(X));
Z=reshape(min(reshape(sum(vecs.*(Q*vecs)),numel(s),[])),size(X,1),[]);
figure(1)
hold off
contour(X,Y,Z,[1 1])
hold on

rad = b^2 - 4*a*d;


% RAD = reshape(msubs(rad,x,XMAT),size(X,1),[]);
% contour(X,Y,RAD,[0 0])

contourSpotless(x'*Q*x,x(1),x(2),[-5 5],[-4 4],[t;x(3:end)],zeros(model.num_states-1,1),1,{'c'});
contourSpotless(b^2 - 4*a*d,x(1),x(2),[-5 5],[-4 4],[t;x(3:end)],zeros(model.num_states-1,1),0,{'b'});
contourSpotless(b,x(1),x(2),[-5 5],[-4 4],[t;x(3:end)],zeros(model.num_states-1,1),0,{'g'});
contourSpotless(2*a-b,x(1),x(2),[-5 5],[-4 4],[t;x(3:end)],zeros(model.num_states-1,1),0,{'r'});
contourSpotless(b^2-4*a*d-(b-2*a)^2,x(1),x(2),[-5 5],[-4 4],[t;x(3:end)],zeros(model.num_states-1,1),0,{'k'});
contourSpotless((2*a+b)^2 - b^2 + 4*a*d,x(1),x(2),[-5 5],[-4 4],[t;x(3:end)],zeros(model.num_states-1,1),0,{'k'});

% %%
% for i=1:size(X,1),
%   for j=1:size(X,2),
%     ZZ(i,j) = min(dmsubs(x'*Q*x,x,repmat([X(i,j);Y(i,j);zeros(4,1)],1,numel(s)) + c*s));
%   end
% end
% %%
% contour(X,Y,ZZ,[1 1])

contourSpotless(x'*S*x,x(1),x(2),[-5 5],[-4 4],[t;x(3:end)],zeros(model.num_states-1,1),1,{'o'});