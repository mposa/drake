function [p] = reconstructDualPoly(x1,dmax1,dmin1,x2,dmax2,dmin2,coef)
%reconstructDualPoly Reconstruct the sprocedure multiplier from the coefs

mono1 = monolist(x1,dmax1);
if dmin1 <= dmax1 & dmin1>0
    s = nchoosek(length(x1) + dmin1-1,dmin1-1);
    mono1 = extsubsref(mono1,s+1:length(mono1));
end

mono2 = monolist(x2,dmax2);
if dmin2 <= dmax2 & dmin2>0
    s = nchoosek(length(x2) + dmin2-1,dmin2-1);
    mono2 = extsubsref(mono2,s+1:length(mono2));
end
monoms = mono1*mono2';
monoms = unique(monoms(:));
p = coef'*monoms;

end