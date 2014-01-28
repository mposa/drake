function [p] = reconstructPoly_spotless(vars,degree,coef)
%reconstructPoly Reconstruct the sprocedure multiplier from the coefs
%  @input vars array of sdpvars
%  @input degree the degree of the multiplier
%  @input coef The coefficients for the multiplier
%  @return p The polynomial multiplier

m = monomials(vars,0:degree);
p = coef'*m;

end