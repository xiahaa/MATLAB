
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
    
    func = {@five_point,                    ...  %% impl1   five_point
            @calibrated_fivepoint,          ...  %% grob
            @E_est_Nister_nongb_warpper,    ...  %% nister
            @peig_five_point_warpper,       ...  %% peig
            @Essential_est_5point,          ...  %% self
            @Essential_est_five_point};        
    N = [5 5 5 5 5 5];%             @sixp_pizarro,                  ...  %% todo  

    %% return type: cell, 9x4, 9x4, 9x4, 9x3, 9x6, 3x3 3x3x4 3x3x4
    
    %% Our 5-point implementation
    Q1 = [0.4964 ,1.0577; ...
                    0.3650, -0.0919; ...
                    -0.5412, 0.0159; ...
                    -0.5239, 0.9467; ...
                    0.3467, 0.5301; ...
                    0.2797, 0.0012; ...
                    -0.1986, 0.0460; ...
                    -0.1622, 0.5347; ...
                    0.0796, 0.2379; ...
                    -0.3946, 0.7969];

    Q2 = [0.7570, 2.7340; ...
                    0.3961, 0.6981; ...
                    -0.6014, 0.7110; ...
                    -0.7385, 2.2712; ...
                    0.4177, 1.2132; ...
                    0.3052, 0.4835; ...
                    -0.2171, 0.5057; ...
                    -0.2059, 1.1583; ...
                    0.0946, 0.7013; ...
                    -0.6236, 3.0253];
    
    Q1 = [Q1';ones(1,size(Q1,1))];
    Q2 = [Q2';ones(1,size(Q2,1))];
    
    res = zeros(9,size(func,2),200);
    ts = zeros(size(func,2),200);
    for k = 1:1:100
        k
%         Q1 = rand(3,20);
%         Q2 = rand(3,20);
        for i = 5:size(func,2)
            id = 1:N(i);
            handle = func{i};
            if i == 1
                tic
                Elist = handle(Q1(:,id)', Q2(:,id)');
                ts(i,k) = toc;
            else
                tic
                Elist = handle(Q1(:,id), Q2(:,id));
                ts(i,k) = toc;
            end
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
            bestE = findbestE(Ecell, Q1(:,:)', Q2(:,:)');
            res(:,i,k) = bestE(:);
        end
    end

    gt = 3;
    err = zeros(size(func,2),1);
    for i = 1:200
        for j = 1:size(func,2)
            v1 = res(:,gt,i);
            v2 = res(:,j,i);
            err(j) = err(j) + norm(v1-v2);
        end
    end
    fontsize = 12;
    xlabels = {'Syb1','Grob','Nie','Peig','Self','Syb2'};
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
