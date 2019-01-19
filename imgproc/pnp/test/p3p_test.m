%
% Verification of Kneip's p3p code
%
%Per-Erik Forss?n, 2016
addpath('../3rdparty/kneipp3p/');
%% simulation of homography decomposition
addpath('../MatrixLieGroup');
addpath('../quaternion');
addpath('../beautiful_plot');
addpath('../');
% Generate three random 3D points
X=rand(3,3)+[0;0;0.5]*[1 1 1];

% Generate a random rotation
R = rand(3,3);
[E,D] = eig(R*R');
R = E*det(E); % ensure det(R)=1
% Generate a random translation
t=rand(3,1)+1.0;

% Get 3D points in global frame by applying the rotation and translation
Xw = [R t]*[X; ones(1,3)];

% Get image points by projection and normalization.
% By "unitary" Kneip means that the points should be normalized to unit vectors.
xp = X*diag(1./sqrt(sum(X.^2,1)));
% Xw = [0.0784    0.5682   -0.0190; ...
%     1.0668    1.1543    1.2836; ...
%     1.3042    0.8841    1.4226];
% xp = [0.4843    0.0534    0.6199; ...
%     0.0863    0.3408    0.1738; ...
%     0.8707    0.9386    0.7652];
% Call p3p to find the four solutions
tic
P = kneip_p3p(Xw,xp);
toc
tic
[R1,t1] = p3p_kneip(Xw,xp,eye(3));
toc
% P
% R1
% t1

% for k=1:4,
%     Pk=P(:,1+(k-1)*4:4+(k-1)*4);
%     if isreal(Pk),
%         fprintf('solution %d:',k);
%         Pk
%     else
%         fprintf('solution %d is complex\n',k);
%     end
% end

fprintf('Ground truth:');
[t R] % Solutions are packed in this way