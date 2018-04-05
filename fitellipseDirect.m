function A = fitellipseDirect()

fid = fopen('C:/Users/xiahaa/workspace/data/webcam/imgmoment/sample.txt','r');
pts = fscanf(fid,'%f %f\n',[2 inf]);

if size(pts,2) == 2
    x = pts(:,1);
    y = pts(:,2);
else
    x = pts(1,:)';
    y = pts(2,:)';
end

%%
D  = [x.*x x.*y y.*y x y ones(size(x,1),1)];
S = D'*D;
C(6,6)=0;C(1,3)=-2;C(3,1)=-2;C(2,2) = 1;

[gevec,geval]=eig(S,C);
[negr,negc]=find(geval<0 & ~isinf(geval));
a = gevec(:,negc);

D1 = [x.*x x.*y y.*y];
D2 = [x y ones(size(x,1),1)];
S1 = D1'*D1;
S2 = D1'*D2;
S3 = D2'*D2;
T = -inv(S3)*S2';
M = S1+S2*T;
M = [M(3,:)./2;-M(2,:);M(1,:)./2];
[evec,eval]=eig(M);
cond = 4*evec(1,:).*evec(3,:)-evec(2,:).^2;
a1 = evec(:,find(cond>0));
astar = [a1;T*a1];

A = a(1);
B = a(2);
C = a(3);
D = a(4);
E = a(5);
F = a(6);

xc = (2*C*D-B*E)/(B*B-4*A*C);yc = (2*A*E-B*D)/(B*B-4*A*C);
a = -sqrt(2*(A*E*E+C*D*D-B*D*E+(B*B-4*A*C)*F)*(A+C+(sqrt((A-C)*(A-C)+B*B))))/(B*B-4*A*C);
b = -sqrt(2*(A*E*E+C*D*D-B*D*E+(B*B-4*A*C)*F)*(A+C-(sqrt((A-C)*(A-C)+B*B))))/(B*B-4*A*C);

if (B==0 && A<C)
    theta = 0;
elseif (B==0 && A>C)
    theta = pi/2;
else
    theta = atan2(C-A-sqrt((A-C)*(A-C)+B*B),B);
end

figure
plot(x,y,'.')
hold on;

t=-pi:0.01:pi;
xe=a*cos(t);
ye=b*sin(t);
xr = xc + cos(theta).*xe+sin(theta).*ye;
yr = yc + -1*sin(theta).*xe+cos(theta).*ye;

plot(xr,yr)

A = astar(1);
B = astar(2);
C = astar(3);
D = astar(4);
E = astar(5);
F = astar(6);

xc = (2*C*D-B*E)/(B*B-4*A*C);yc = (2*A*E-B*D)/(B*B-4*A*C);
a = -sqrt(2*(A*E*E+C*D*D-B*D*E+(B*B-4*A*C)*F)*(A+C+(sqrt((A-C)*(A-C)+B*B))))/(B*B-4*A*C);
b = -sqrt(2*(A*E*E+C*D*D-B*D*E+(B*B-4*A*C)*F)*(A+C-(sqrt((A-C)*(A-C)+B*B))))/(B*B-4*A*C);

if (B==0 && A<C)
    theta = 0;
elseif (B==0 && A>C)
    theta = pi/2;
else
    theta = atan2(C-A-sqrt((A-C)*(A-C)+B*B),B);
end
t=-pi:0.01:pi;
xe=a*cos(t);
ye=b*sin(t);
xr = xc + cos(theta).*xe+sin(theta).*ye;
yr = yc + -1*sin(theta).*xe+cos(theta).*ye;

plot(xr,yr,'k+')

return