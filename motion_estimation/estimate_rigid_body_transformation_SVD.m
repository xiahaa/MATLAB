function [Ropt,topt] = estimate_rigid_body_transformation_SVD(ptsrc, ptdst)
%% SVD 

ptsrcmean = mean(ptsrc,2);
ptdstmean = mean(ptdst,2);

ptsrcrefine = ptsrc - repmat(ptsrcmean, 1, size(ptsrc,2));
ptdstrefine = ptdst - repmat(ptdstmean, 1, size(ptsrc,2));

Y = ptdstrefine';
X = ptsrcrefine;
S = X*Y;
[U,~,V] = svd(S);

D = V*U';
if det(D) < 0
    Ropt = V*[1 0 0;0 1 0;0 0 -1]*U';
else
    Ropt = V*U';
end
topt = ptdstmean - Ropt * ptsrcmean;


end