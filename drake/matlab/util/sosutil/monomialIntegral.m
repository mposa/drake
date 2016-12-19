function y = monomialIntegral(alphas,sphere_inds,A_diag,box_inds,box_lims)
if isempty(box_inds)
  y = 1;
else
  box_pows = alphas(box_inds) + 1;
  y = prod(1./box_pows .*(box_lims(:,2)'.^box_pows - box_lims(:,1)'.^box_pows));
end

if ~isempty(sphere_inds)
  n_sphere = length(sphere_inds);
  sphere_alphas = alphas(sphere_inds);
  
  sphere_betas = 0.5*(sphere_alphas + 1);
  Ra = (1.^(sum(sphere_alphas,2) + n_sphere))./(sum(sphere_alphas,2) + n_sphere);
  IS = 2*prod(gamma(sphere_betas),2)./(gamma(sum(sphere_betas,2)));
  l = Ra.*IS;
  alphaszero = (mod(sphere_alphas,2) ~= 0);
  alphaszero = any(alphaszero,2);
  l(alphaszero) = 0;
  
  y = y*l*prod(repmat(A_diag,size(sphere_alphas,1),1).^(sphere_alphas+1),2);  
  
end
end