    clc;
    close all;
    clear all;

    %% simulation of homography decomposition
    addpath('../MatrixLieGroup');
    addpath('../quaternion');
    addpath('../beautiful_plot');
    addpath('../');
    T1 = fakeRT();
    
    N = 200;
    p = rand([3,N]) * 5 - 2.5;
    p(3,:) = 5;
    p(1,:) = p(1,:);
    p(2,:) = p(2,:);
    
    K = [1 0 0;0 1 0;0 0 1];
    
    [uv1, in1] = proj(T1, p, K);
    im = zeros(240,320);
    
    q1 = uv1(:,in1);
    P = p(:,in1);
    
    pr = T1(1:3,1:3)*p + repmat(T1(1:3,4),1,N);
    pr = pr(:,in1);
    
    id = randperm(size(q1,2),3);
    pr(:,id)
    T1
    
    q1n = K\q1;
    
    [R, t] = p3p_kneip(P(:,:), q1, K);
%     [R, t] = orthogonal_iterative_optimization(P(:,:), q1n);
    
    minerr = 1e6;
    minid = 0;
    P = [P;ones(1,size(P,2))];
    for i = 1:size(R,3)
        P1 = K*([R(1:3,1:3,i) t(1:3,1,i)]);
        uv1rep = P1*P;
        uv1rep = uv1rep./uv1rep(3,:);
        err = uv1rep - q1;
        avgerr = sum(diag(err'*err)) / size(q1,2);
        if avgerr < minerr
            minerr = avgerr;
            minid = i;
        end
    end
    R(:,:,minid)
    t(:,:,minid)
    
    
    

    
    
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