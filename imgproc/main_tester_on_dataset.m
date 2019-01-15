function main_tester_on_dataset
    clc;close all;clear all;
    
    tmp = matlab.desktop.editor.getActive;
    cd(fileparts(tmp.Filename));
    
    %% data directory
    addpath C:\Users\xiahaa\Documents\MATLAB\3rdparty\allyourposesrours-master\simulator;
    %% five point implemented by someone from CMU
    addpath C:\Users\xiahaa\Documents\MATLAB\3rdparty\allyourposesrours-master\paper1_impl;
    %% five points
    addpath C:\Users\xiahaa\Documents\MATLAB\3rdparty\E_estim\E_estim\Grobner_five_point
    addpath C:\Users\xiahaa\Documents\MATLAB\3rdparty\E_estim\E_estim\Nister_five_point
    addpath C:\Users\xiahaa\Documents\MATLAB\3rdparty\E_estim\E_estim\Polynomial_Eig_five_point
    addpath C:\Users\xiahaa\Documents\MATLAB\3rdparty\E_estim\E_estim\severn_point
    addpath C:\Users\xiahaa\Documents\MATLAB\3rdparty\E_estim\E_estim\sixp
    addpath C:\Users\xiahaa\Documents\MATLAB\3rdparty\E_estim\E_estim\utils
    %% addpath 
    addpath './essential'
    
    %% Load the camera calibrations
    load('templeSparseRing_data.mat');
    [ P, K, R, t ] = extractKRt(data);
    %% Generate a 3D point cloud
    [ pointList ] = generatePoints( );
    %% Project the 3D points onto each camera's normalized image plane
    [ pts1, pts2 ] = projectPoints( pointList, P{1}, P{2} );
    
    % Homogeneous version of the pointss
    h_pts1 = [pts1 ones(size(pts1,1),1)];
    h_pts2 = [pts2 ones(size(pts2,1),1)];
    
    %% coordinate normalization
    n_pts1 = K{1} \ h_pts1';
    n_pts2 = K{2} \ h_pts2';
    n_pts1 = n_pts1';
    n_pts2 = n_pts2';
    
    %% Perfect backprojection original projection matrices
    % backprojected points do not need to be premultiplied by K1 and K2 since
    % the function backproject makes use of the entire projection matrix which
    % includes camera intrinsics and extrinsics.
    [ idealReprojPoints ] = backproject( h_pts1, h_pts2, P{1}, P{2} );
    figure();
    plot3(idealReprojPoints(:,1), idealReprojPoints(:,2), idealReprojPoints(:,3), '.');
    title('Gnd truth projection matrices');
    axis equal;
    
    %% Backprojection with error-free E and calculated projection matrices
    Etrue = findGndTruthE(R{1}, t{1}, R{2}, t{2});
    reprojPoints = triangulateWithE(h_pts1, h_pts2, Etrue, K{1}, K{2});
    
%     func = {@five_point,          ...  %% impl1
%             @calibrated_fivepoint,          ...  %% grob
%             @E_est_Nister_nongb_warpper,    ...  %% nister
%             @peig_five_point_warpper,       ...  %% peig
%             @sevenp,                        ...  %% 7 point
%             @Essential_est_8_point,         ...  %% 8 point2
%             @Essential_est_5point,          ...  %% self
%             @Essential_est_five_point};        
    func = {@calibrated_fivepoint,          ...  %% grob
            @E_est_Nister_nongb_warpper,    ...  %% nister
            @Essential_est_5point,          ...  %% self
            @Essential_est_8_point,         ...  %% 8 point2
            @sevenp,                        ...  %% 7 point
            @sixp_pizarro,                  ...  %% todo
     };  

    N = [5 5 5 8 7 6];%             

    %% return type: cell, 9x4, 9x4, 9x4, 9x3, 9x6, 3x3 3x3x4 3x3x4
    ts = zeros(6,100);
    res = zeros(9,6,100);
    
    for k = 1:100
        k
        %% Our 5-point implementation
        Ns = randperm(size(n_pts1,1),20);
        for i = 1:size(func,2)
            id = Ns(1:N(i));
            handle = func{i};
%             if i == 1
%                 tic
%                 Elist = handle(n_pts1(id,:), n_pts2(id,:));
%                 t1 = toc;
%             else
                tic
                Elist = handle(n_pts1(id,:)', n_pts2(id,:)');
                t1 = toc;
%             end
            %% find the best one
            if ~iscell(Elist)
                if size(Elist,1) == 9
                    nsol = size(Elist,2);
                    Emat = permute(reshape(Elist, 3, 3, nsol), [2,1,3]);
    %                 Ecell = mat2cell(Emat, [3, 3], ones(1, nsol));
                else
                    Emat = Elist;
                end

                nsol = size(Emat,3);
                for j = 1:nsol
                    Ecell{j} = Emat(:,:,j) ./ Emat(3,3,j);
                end
            else
                Ecell = Elist;
                for j = 1:size(Ecell,2)
                    E1 = Ecell{j};
                    Ecell{j} = E1(:,:) ./ E1(3,3);
                end
            end

            bestE = findbestE(Ecell, n_pts1(:,:), n_pts2(:,:));
            res(:,i,k) = bestE(:);
            ts(i,k) = t1;
        end
    end
    Etrue = Etrue./Etrue(3,3);
    gt = 3;
    err = zeros(size(func,2),1);
    for i = 1:100
        for j = 1:size(func,2)
            v1 = Etrue(:);
            v2 = res(:,j,i);
            err(j) = err(j) + norm(v1-v2);
        end
    end
    fontsize = 12;
    xlabels = {'Grob','Nie','Self','8pt','7pt','6pt'};
    figure(1);plot(err,'-o');
    xticks([1 2 3 4 5 6]);
    xticklabels(xlabels);
    ylabel('$Error: (s)$','Interpreter','latex','FontSize',fontsize);
    title('Performance Validation','Interpreter','latex','FontSize',12)
    set(gca,'TickLabelInterpreter','latex')
    set(gca,'FontSize',12)
    set(gca,'XTickLabelRotation',45);
    figure(2);plot(mean(ts,2),'-o');
    xticks([1 2 3 4 5 6]);
    xticklabels(xlabels);
    ylabel('$Time: (s)$','Interpreter','latex','FontSize',fontsize);
    title('Runtime Evaluation','Interpreter','latex','FontSize',12)
    set(gca,'TickLabelInterpreter','latex')
    set(gca,'FontSize',12)
    set(gca,'XTickLabelRotation',45);
    
    triangulateWithE(h_pts1, h_pts2, bestE, K{1}, K{2});
end

function bestE = findbestE(EList, n_pts1, n_pts2)
    %% Check error
    minSumError = inf;
    for i=1:length(EList)
        E1 = EList{i};
%         i;
        % Check determinant constraint! 
%         error_determinate = det( E1);
        % Check trace constraint
%         error_trace = 2 *E1*transpose(E1)*E1 -trace( E1*transpose(E1))*E1;
        % Check reprojection errors
        error_reprojection = diag( n_pts2 * E1 * n_pts1');
        sum_error_proj = sum(abs(error_reprojection));
        % Find E with the smallest error
        if (sum_error_proj < minSumError)
            minSumError = sum_error_proj;
            bestE = E1;
        end
    end
    % Show the points in 3D backprojected with the best E.
end

function [ points ] = triangulateWithE( pts_1, pts_2, E, K1, K2 )
    %triangulate Finds the 3D points given the 2D correspondences and camera
    %intrinsics.
    M1 = K1 * eye(3,4);
    % camera2 generates 4 possible camera projection matrices. We will display
    % with each of them.
    %M2 = camera2(E', K2);
    % figure();
    %
    % for i = 1:4
    %     points = backproject( pts_1, pts_2, M1, M2{i} );
    %
    %     subplot(2,2,i);
    %     plot3(points(:,1), points(:,2), points(:,3), '.');
    %     axis equal;
    %
    % end
    M2 = K2 * robustCameraRecovery(E, pts_1, pts_2);
    figure();
    points = backproject( pts_1, pts_2, M1, M2);
    plot3(points(:,1), points(:,2), points(:,3), '.');
    axis equal;
end

function E = findGndTruthE(R1,t1,R2,t2)
    T1 = [R1 t1;[0 0 0 1]];
    T2 = [R2 t2;[0 0 0 1]];
    T1to2 = T2 * inv(T1);
    R1to2 = T1to2(1:3,1:3);
    t1to2 = T1to2(1:3,4);
    E = skewm(t1to2)*R1to2;
end

function [ P, K, R, t ] = extractKRt(data)
    numPos = size(data,1);
    K = cell(numPos, 1);
    R = cell(numPos, 1);
    t = cell(numPos, 1);
    P = cell(numPos, 1);
    for i = 1:numPos
        % K and R are transposed because the reshape function places elements
        % into the new matrix column-wise but the original data is stored
        % row-wise.
        K{i} = reshape(data(i,1:9),3,3)';
        R{i} = reshape(data(i,10:18),3,3)';
        t{i} = reshape(data(i,19:21),3,1);
        P{i}  = K{i} * [ R{i} t{i}];
    end
end

function [ pointList ] = generatePoints( )
    %generatePoints Generates a list of points
    %   Detailed explanation goes here
    % A cube:
    % I'm adding the corner points first to ensure that the 5-point used to
    % calculate E are not on a plane. All 5-points on a plane seemed to cause
    % problems earlier.
    pointList = [   -0.5, -0.5, -0.5;
                    -0.5, -0.5,  0.5;
                    -0.5,  0.5, -0.5;
                    -0.5,  0.5,  0.5;
                     0.5, -0.5, -0.5;
                     0.5, -0.5,  0.5;
                     0.5,  0.5, -0.5;
                     0.5,  0.5,  0.5];
    % Add interior points
    [X, Y, Z] = ndgrid(-0.5:0.1:0.5, -0.5:0.1:0.5, -0.5:0.1:0.5);
    pointList = [   pointList;
                    X(:), Y(:), Z(:)];
    % DEBUG: View points
%     {
%     figure, plot3(pointList(:,1), pointList(:,2), pointList(:,3), '.');
%     title('Original 3D points');
%     }
end

function [ points ] = backproject( pts_1, pts_2, M1, M2 )
    %triangulate Finds the 3D points given the 2D correspondences and camera
    % projection matrices.

    % points: n-by-3 matrix of 3D points in the form (x, y, z) with each row
    % representing one point

    % Preallocate
    points = zeros(length(pts_1), 4);

    for i = 1:size(pts_1, 1)
        % Method taken from section 12.2 Linear triangulation methods in
        % Hartley R., Zisserman A., Multiple View Geometry in Computer Vision  
        x1 = pts_1(i,1);
        y1 = pts_1(i,2);
        x2 = pts_2(i,1);
        y2 = pts_2(i,2);
        % A is a 4-by-4 matrix    
        A = [  x1 * M1(3,:) - M1(1,:);
                y1 * M1(3,:) - M1(2,:);
                x2 * M2(3,:) - M2(1,:);
                y2 * M2(3,:) - M2(2,:)];
        % Ah = 0, where h is the 3D homogeneous coordinate of the real world
        % point
        [~, ~, V] = svd(A);
        % Take the last column of V
        h = V(:,end);
        % Scale by the homogeneous component
        h = h / h(end);
        % Place into the point list
        points(i,:) = h';
    end
    % Trim off the homogenous component
    points = points(:,1:3);
end

