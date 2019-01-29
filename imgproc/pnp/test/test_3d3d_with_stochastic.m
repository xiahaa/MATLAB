    clc;
    close all;
    clear all;

    %% simulation of homography decomposition
    addpath('../../../MatrixLieGroup');
    addpath('../../../quaternion');
    addpath('../../../beautiful_plot');
    addpath('../');
    T1 = fakeRT();
    
    N = 10;
    p = rand([3,N]) * 5 - 2.5;
%     p(3,:) = 5;%p(3,:);
    p(1,:) = p(1,:);
    p(2,:) = p(2,:);
    
    q = T1(1:3,1:3)*p + repmat(T1(1:3,4),1,N);
    q = q(:,1:N);
    T1
    
    [R,t] = pose_est_without_correspondance(p, q, T1(1:3,1:3));

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
    euler(1) = (rand(1)*pi/2 - pi/4);
    euler(2) = (rand(1)*pi/2 - pi/4);
    euler(3) = (rand(1)*2*pi - pi);
    R1 = euler2rot(euler(1),euler(2),euler(3));
    t1 = rand([1,3]) * 2;
    t1 = t1';
    t1(1) = t1(1);
    t1(2) = t1(2);
    t1(3) = 1;%t1(3);
    T = [R1 t1;[0 0 0 1]];
end

function [R,t] = pose_est_without_correspondance(p, q, Rx)
    % Plan: can we obtain PnP solution without correspondance or with
    % less than 3 correspondance but with a large batch of 3D and 2D
    % points. We know they will match with each other. However, we don't
    % know which one should match to which one.
    %
    nq = size(q,2);
    np = size(p,2);
    %% step 1: forme SO(3) for qn and P
    set_q = nchoosek(1:nq,3);
    set_p = nchoosek(1:np,3);
    
    SO3_q = formeSO3(q, set_q);
    SO3_p = formeSO3(p, set_p);
    
    %% step 2: forme relative SO(3) for qn and P
    raw_SO3_q = forme_relative_SO3(SO3_q);
    raw_SO3_p = forme_relative_SO3(SO3_p);
    
    relative_SO3_q = rankingSO3(raw_SO3_q);
    relative_SO3_p = rankingSO3(raw_SO3_p);
%     relative_SO3_q = raw_SO3_q;
%     relative_SO3_p = raw_SO3_p;
    
    %% just for debug
%     for i = 1:1:size(relative_SO3_q,3)
%         sel = randperm(nq,2);
%         ii = sel(1);jj = sel(2);
% %         u1 = q(:,ii)-q(:,jj);u1=u1./norm(u1);
% %         v1 = p(:,ii)-p(:,jj);v1=v1./norm(v1);
% %         norm(u1-Rx*v1)
%         R1 = relative_SO3_q(:,:,i)*Rx - Rx*relative_SO3_p(:,:,i);
%         err(i) = norm(R1,'fro');
%     end
    
    %% tempotary solution: SO3 to SE3
    SE3_q = SO3toSE3(relative_SO3_q);
    SE3_p = SO3toSE3(relative_SO3_p);
    
    %% step 3: now, we know: ideally, for those who can be matched, AX = XB should be maintained.
    %% which equals to a stochastic hand eye calibration problem.
    addpath C:\Users\xiahaa\Documents\MATLAB\hand_eye_calib\solver\
    addpath C:\Users\xiahaa\Documents\MATLAB\hand_eye_calib\stochastics
    addpath C:\Users\xiahaa\Documents\MATLAB\hand_eye_calib\stochastics\mean
    addpath C:\Users\xiahaa\Documents\MATLAB\hand_eye_calib\stochastics\utils
    addpath C:\Users\xiahaa\Documents\MATLAB\MatrixLieGroup\barfoot_tro14\
%     [a1,a2,a3]  = size(SE3_q);
%     A_noise_mex = reshape(SE3_q, a1, a2*a3);
%     [a1,a2,a3]  = size(SE3_p);
%     B_mex = reshape(SE3_p, a1, a2*a3);
    opt = 2;
    [X_solved, MA, MB, SigA, SigB] = batchSolveNew(SE3_q, SE3_p, opt);
%     dT = sol_manifold_opt_SE3(SE3_q, SE3_p, size(SE3_q,1));
%     X_solved
    
    
end

function goodSO3s = rankingSO3(SO3s)
    angles = zeros(1,size(SO3s,3));
    for i = 1:size(SO3s,3)
        angles(i) = abs(norm(rot2vec(SO3s(:,:,i))));
    end
    [~,id] = sort(angles,'descend');
    num = 100;
    if size(SO3s,3) < 50
        num = size(SO3s,3);
    end
    goodSO3s = SO3s(:,:,id(1:num));
end

function SE3s = SO3toSE3(SO3s)
    SE3s = zeros(4,4,size(SO3s,3));
    for i = 1:size(SO3s,3)
        SE3s(:,:,i) = [SO3s(:,:,i) [0;0;0];[0,0,0,1]];
    end
end

function SE3s = SO3toSE32(SO3s)
    SE3s = zeros(size(SO3s,3),4,4);
    for i = 1:size(SO3s,3)
        SE3s(i,:,:) = [SO3s(:,:,i) [0;0;0];[0,0,0,1]];
    end
end

function SO3s = formeSO3(p, setp)
    n_so3 = size(setp,1);
    SO3s = zeros(3,3,n_so3);
    for i=1:n_so3
        ii = setp(i,1);
        jj = setp(i,2);
        kk = setp(i,3);
        v1 = p(:,ii)-p(:,jj); 
        v2 = p(:,ii)-p(:,kk);
        v1 = v1./norm(v1);v2 = v2./norm(v2);
        v3 = cross(v1,v2);v3 = v3./norm(v3);
        v2 = cross(v3,v1);v2 = v2./norm(v2);
        SO3s(:,:,i) = [v1 v2 v3];
    end
end

function relative_SO3s = forme_relative_SO3(SO3s)
    nSO3s = size(SO3s,3);
    set_r_SO3s = nchoosek(1:nSO3s, 2);
    n_r_SO3s = size(set_r_SO3s,1);
    relative_SO3s = zeros(3,3,n_r_SO3s);
    for i = 1:1:n_r_SO3s
        ii = set_r_SO3s(i,1);
        jj = set_r_SO3s(i,2);
        relative_SO3s(:,:,i) = SO3s(:,:,ii)*SO3s(:,:,jj)';
    end
end