function main_homography_decomposition_test
close all;
    %% simulation of homography decomposition
    addpath('../MatrixLieGroup');
    addpath('../quaternion');
    addpath('../beautiful_plot');
    addpath('./homography')

    T1 = fakeRT();
    T2 = fakeRT();

    figure()
    DrawAxis(T1, 0.5, 'r', 'g', 'b');
    hold on;
    DrawAxis(T2, 0.5, 'm', 'c', 'y');

    N = 200;
    p = rand([3,N]) * 5 - 2.5;
    p(3,:) = 5;
    p(1,:) = p(1,:);
    p(2,:) = p(2,:);
    for i = 1:size(p,2)
        plot3(p(1,:),p(2,:),p(3,:),'bo','MarkerSize',2);
    end
    axis equal;
    view(3);
    xlabel('x');ylabel('y');zlabel('z');

    K = [200 0 160;0 200 120;0 0 1];

    [uv1, in1] = proj(T1, p, K);
    im = zeros(240,320);

    [uv2, in2] = proj(T2, p, K);
    im2 = zeros(240,320);
%     figure()
%     for i = 1:size(uv2,2)
%         im2(int32(uv2(2,i)),int32(uv2(1,i))) = 255;
%     end
%     imshow(im2,[]);

    %% matching
    matching = in1 & in2;
    if isempty(matching)
        warning('No overlap');
    end
    Tr = T2*inv(T1);
    n = T1(1:3,1:3) * [0;0;1];
    ps1 = T1(1:3,1:3) * p(:,1) - T1(1:3,4);
%     ps2 = T1(1:3,1:3) * p(:,2) + T1(1:3,4);
    d = dot(ps1,n);
%     d2 = dot(-ps2,n);

    uv1n = inv(K)*uv1(:,matching);
    uv2n = inv(K)*uv2(:,matching);
    H1 = homography_est(uv1n, uv2n);
    H1 = H1./H1(3,3);

%     H1 = [0.857462971105479  -0.514545676478740   0.238814648426845; ...
%           0.514545676478740   0.857462971105479  -0.111976619581788; ...
%           0.000000000000000   0.000000000000000   1.000000000000000];

    H1 = normalizeHomography(H1);
    [R1,t1,n1] = homo_decom_svd_olivier(H1, uv1n);
    [R2,t2,n2] = homo_decom_malis(H1, uv1n);
    [R3,t3,n3] = homo_decom_svd_mayi(H1, uv1n);
    t1
    t2
    t3
%
%
end

function T = fakeRT()
    euler(1) = (rand(1)*pi/2 - pi/4)*0;
    euler(2) = (rand(1)*pi/2 - pi/4)*0;
    euler(3) = (rand(1)*2*pi - pi);
    R1 = euler2rot(euler(1),euler(2),euler(3));
    t1 = rand([1,3]) * 2;
    t1 = t1';
    t1(1) = t1(1);
    t1(2) = t1(2);
    t1(3) = 1;%t1(3);
    T = [R1 t1;[0 0 0 1]];
end

function [uv1, in] = proj(T, p, K)
    P1 = K*([T(1:3,1:3) T(1:3,4)]);
    phomo = [p;ones(1,size(p,2))];
    uv1 = P1*phomo;
    uv1 = uv1./uv1(3,:);
    in = uv1(1,:) > 0 & uv1(1,:) < 321 & uv1(2,:) > 0 & uv1(2,:) < 241;
%     uv1 = uv1(:,in);
end
