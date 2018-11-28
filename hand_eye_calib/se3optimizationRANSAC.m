function dT = se3optimizationRANSAC(TA, TB, N)
%% input should be Nx4x4 3D matrices.    
    addpath('./solver/');
    
    %% 1st, prepare data, transform to cell array
    for j = 1:N
        Ts = TA(j,:,:);Ts = reshape(Ts,4,4);
        A{j} = Ts;invA{j} = inv(A);
        Ts = TB(j,:,:);Ts = reshape(Ts,4,4);
        B{j} = Ts;
    end

    %% RANSAC parameters
    bestCost = -1e8;
    bestInliers = [];
    sampleCNT = 2;
    inlierThreshold = 10 * pi / 180.0;
    bestdists = [];
    maxIter = 50;
    iter = 1;

    while iter < maxIter
        %% random sampling
        is =randperm(N,sampleCNT);
        sA = A{is};
        sB = B{is};
        %% estimation
        dTopts = se3optimization(sA, sB, sampleCNT);
        %% verification
        dists = zeros(N,1);
        inlierflag = ones(N,1);
        cnt = N;
        for i = 1:N
            dist = SE3distance(invA{i}*dTopts*B{i}*dTopts',eye(4));
            dists(i) = dist;
            if dist >= inlierThreshold
                inlierflag(i) = 0;
                cnt = cnt - 1;
            end
        end
        cost = cnt;
        
        if bestCost < cost
            bestCost = cost;
            bestdists = dists;
            bestInliers = inlierflag;
        end 
        iter = iter + 1;
    end
    
    %% refine with all inliers
    ids = 1:N;ids = ids';
    ids(bestInliers<1) = [];
    for j = 1:bestCost
        Ain{j} = A{ids(j)};
        Bin{j} = B{ids(j)};
    end
    dT = se3optimization(Ain, Bin, bestCost);%% estimation
end

function dist = SE3distance(SE31,SE32)
    dist = norm(tran2vec(SE31'*SE32));
end
