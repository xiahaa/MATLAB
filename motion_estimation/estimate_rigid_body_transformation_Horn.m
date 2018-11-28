function [Ropt,topt] = estimate_rigid_body_transformation_Horn(ptsrc, ptdst)

ptsrcmean = mean(ptsrc')';
ptdstmean = mean(ptdst')';

ptsrcrefine = ptsrc - repmat(ptsrcmean, 1, size(ptsrc,2));
ptdstrefine = ptdst - repmat(ptdstmean, 1, size(ptsrc,2));
M = [ptsrcrefine(1,:)*ptdstrefine(1,:)' ptsrcrefine(1,:)*ptdstrefine(2,:)' ptsrcrefine(1,:)*ptdstrefine(3,:)';...
     ptsrcrefine(2,:)*ptdstrefine(1,:)' ptsrcrefine(2,:)*ptdstrefine(2,:)' ptsrcrefine(2,:)*ptdstrefine(3,:)';...
     ptsrcrefine(3,:)*ptdstrefine(1,:)' ptsrcrefine(3,:)*ptdstrefine(2,:)' ptsrcrefine(3,:)*ptdstrefine(3,:)'];

% 3 N
N = [M(1,1)+M(2,2)+M(3,3) M(2,3)-M(3,2) M(3,1)-M(1,3) M(1,2)-M(2,1);...
     M(2,3)-M(3,2) M(1,1)-M(2,2)-M(3,3) M(1,2)+M(2,1) M(3,1)+M(1,3);...
     M(3,1)-M(1,3) M(1,2)+M(2,1) -M(1,1)+M(2,2)-M(3,3) M(2,3)+M(3,2);...
     M(1,2)-M(2,1) M(3,1)+M(1,3) M(2,3)+M(3,2) -M(1,1)-M(2,2)+M(3,3)];
 
[evec, eval] = eig(N);
evalvec = [eval(1,1);eval(2,2);eval(3,3);eval(4,4)];
[maxeval, maxid] = max(evalvec);
q = evec(:,maxid);

Ropt = q2r(q);

topt = ptdstmean - Ropt * ptsrcmean; 
 
end