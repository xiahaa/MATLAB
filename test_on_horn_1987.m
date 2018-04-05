function test_on_horn_1987

ptsrc = rand(3,1000).*10 - 5;
yaw = 57 * pi / 180.0;
pitch = 60 * pi / 180.0;
roll = -30 * pi / 180.0;
R = rot(roll,pitch,yaw);
t = [-3;-2;5];
ptdst = R * ptsrc;
ptdst(1,:) = ptdst(1,:)+t(1);
ptdst(2,:) = ptdst(2,:)+t(2);
ptdst(3,:) = ptdst(3,:)+t(3);

% 1 centroid
ptsrcmean = mean(ptsrc')';
ptdstmean = mean(ptdst')';

% 2 M
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

%% SVD 
Y = ptdstrefine';
X = ptsrcrefine;
S = X*Y;
[U,S,V] = svd(S);

D = V*U';
if det(D) < 0
    Rsvd = V*[1 0 0;0 1 0;0 0 -1]*U';
else
    Rsvd = V*U';
end
tsvd = ptdstmean - Rsvd * ptsrcmean;


return

function R = q2r(q)

R = [q(1)*q(1)+q(2)*q(2)-q(3)*q(3)-q(4)*q(4) 2*(q(2)*q(3)-q(1)*q(4)) 2*(q(2)*q(4)+q(1)*q(3));...
     2*(q(2)*q(3)+q(1)*q(4)) q(1)*q(1)-q(2)*q(2)+q(3)*q(3)-q(4)*q(4) 2*(q(3)*q(4)-q(1)*q(2));...
     2*(q(2)*q(4)-q(1)*q(3)) 2*(q(3)*q(4)+q(1)*q(2)) q(1)*q(1)-q(2)*q(2)-q(3)*q(3)+q(4)*q(4)];

return

function R = rot(roll, pitch, yaw)

R1 = [cos(yaw) sin(yaw) 0;...
     -sin(yaw) cos(yaw) 0;...
     0 0 1];
R2 = [cos(pitch) 0 -sin(pitch);0 1 0;sin(pitch) 0 cos(pitch)];

R3 = [1 0 0;0 cos(roll) sin(roll);0 -sin(roll) cos(roll)];

R = R3*R2*R1;

return