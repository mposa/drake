x0 = linspace(-1,1,1000);
% y0 = x0.^2;
y0 = exp(-1./(x0.^2));

theta = linspace(-.5,.5,1000);

pts = zeros(2,length(x0),length(theta));
min_ind = zeros(length(theta),1);
min_val  = zeros(length(theta),1);
for i=1:length(theta),
  R = [cos(theta(i)) -sin(theta(i)); sin(theta(i)) cos(theta(i))];
  for j=1:length(x0),
    pts(:,j,i) = R*[x0(j);y0(j)];
  end
  [min_val(i),min_ind(i)] = min(pts(2,:,i));
%   tmp = pts(:,:,i) - repmat([0;-1],1,length(x0));
%   [min_val(i),min_ind(i)] = min(sum(tmp.*tmp));
end

l = sqrt(x0(min_ind).^2 + y0(min_ind).^2);
figure(1)
% plot(theta,l,theta,x0(min_ind).^2./l)
plot(theta,l)

figure(2)
plot(x0,x0.^2./sqrt(x0.^2 + x0.^4))

figure(3)
plot(x0,y0)