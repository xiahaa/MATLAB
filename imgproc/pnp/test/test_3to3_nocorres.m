clc;
clear all;
close all;

%% simulation of homography decomposition
addpath('../../../MatrixLieGroup');
addpath('../../../quaternion');
addpath('../../../beautiful_plot');
addpath('../');
addpath('./eurasip/');
T1 = fakeRT();
tform = affine3d(T1');
% ptCloud = pcread('teapot.ply');
% ptCloud = pcread('C:\Users\xiahaa\3rdparty\opencv-3.1.0\samples\cpp\tutorial_code\viz\bunny.ply');
% gridStep = 0.02;
% ptCloud = pcdownsample(ptCloud,'gridAverage',gridStep);
% ptCloudTformed = pctransform(ptCloud,tform);
% figure
% pcshow(ptCloud); hold on;
% pcshow(ptCloudTformed); 
% tic
% tform1 = pcregistericp(ptCloudTformed,ptCloud,'Extrapolate',true);
% toc
% tform2 = invert(tform1);
% disp(tform2.T);
% p = ptCloud.Location;
% q = ptCloudTformed.Location;
% p = p';
% q = q';

N = 50;
p = rand([3,N]) * 5 - 2.5;
% p(3,:) = 5;%p(3,:);
p(1,:) = p(1,:);
p(2,:) = p(2,:);
s = 1;
q = s.*(T1(1:3,1:3)*p + repmat(T1(1:3,4),1,N));
q = q(:,1:N) + rand([3,N]).*0.0;

font_size = 14;
blue_color = [0 116 186]/255;
orange_color = [223 80 35]/255;
plot3(q(1,:),q(2,:),q(3,:),'o','MarkerSize',5,'MarkerEdgeColor',orange_color,'MarkerFaceColor',orange_color)
hold on;
plot3(p(1,:),p(2,:),p(3,:),'o','MarkerSize',5,'MarkerEdgeColor',blue_color,'MarkerFaceColor',blue_color)
xlabel('$x$','Interpreter','latex')
ylabel('$y$','Interpreter','latex')
zlabel('$z$','Interpreter','latex')

disp(T1)

tic    
%     nq = size(q,2);
%     np = size(p,2);
%     %% step 1: forme SO(3) for qn and P
%     set_q = nchoosek(1:nq,2);
%     set_p = nchoosek(1:np,2);
%     
%     SO3_q = formeSO3(q, set_q);
%     SO3_p = formeSO3(p, set_p);
%    
% %     SO3_q = rankingSO3(SO3_q);
% %     SO3_p = rankingSO3(SO3_p);
%     
%     R1 = eye(3);
% 
%     SO3_p1 = SO3_p;
%     
%     Mq1 = mean_1st_order(SO3_q);
%     Mq2 = mean_iterative_kron(SO3_q,Mq1);
%     Mq3 = FNS_iterative(SO3_q,Mq1);
% 
%     Mp1 = mean_1st_order(SO3_p1);
%     Mp2 = mean_iterative_kron(SO3_p1,Mp1);
%     Mp3 = FNS_iterative(SO3_p1,Mp1);
%     
%     R3 = Mq3*inv(Mp3)
%     R2 = Mq2*inv(Mp2)
%     R = Mq1*inv(Mp1)
    [R1,t1,R2,t2] = pose_estimation_correspondence_free(p, q);
toc
    [R1 t1;[0 0 0 1]]
    [R2 t2;[0 0 0 1]]
%     norm(logm(R2'*T1(1:3,1:3)),'fro')*57.3

    
function T = fakeRT()
    euler(1) = (rand(1)*pi/2 - pi/4);
    euler(2) = (rand(1)*pi/2 - pi/4);
    euler(3) = (rand(1)*2*pi - pi);
    R1 = euler2rot(euler(1),euler(2),euler(3));
    t1 = rand([1,3]) * 2 + 3;
    t1 = t1';
    t1(1) = t1(1);
    t1(2) = t1(2);
    t1(3) = 1;%t1(3);
    T = [R1 t1;[0 0 0 1]];
end

function goodSO3s = rankingSO3(SO3s)
    angles = zeros(1,size(SO3s,3));
    for i = 1:size(SO3s,3)
        angles(i) = abs(norm(rot2vec(SO3s(:,:,i))));
    end
    [angless,id] = sort(angles);
    num = 100;
    if size(SO3s,3) < 100
        num = size(SO3s,3);
    end
    goodSO3s = SO3s(:,:,id(1:num));
end







