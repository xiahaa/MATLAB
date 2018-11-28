function dRopt = so3optimizationRANSAC(Rii, Rci, N)
    %% 1st, prepare data
    for j = 1:N
        Rs = Rii{j};
        Ri{j} = Rs;%reshape(Rs,3,3,1);
        Rs = Rci{j};
        Rc{j} = Rs;%reshape(Rs,3,3,1);
        invRc{j} = inv(Rs);
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
        is =randi(N,sampleCNT,1);
        for i=1:sampleCNT
            Ris{i} = Ri{is(i)};
            Rcs{i} = Rc{is(i)};
        end
        
        dRopts = so3optimization(Ris, Rcs, sampleCNT);%% estimation
        dists = zeros(N,1);
        inlierflag = ones(N,1);
        cnt = N;
        for i = 1:N
            dist = SO3distance(invRc{i}*dRopts*Ri{i}*dRopts',eye(3));
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
        Rin{j} = Ri{ids(j)};%reshape(Rs,3,3,1);
        Rcn{j} = Rc{ids(j)};
    end
    dRopt = so3optimization(Rin, Rcn, bestCost);%% estimation
end

function dist = SO3distance(SO31,SO32)
    dist = norm(rot2vec(SO31'*SO32));
end

function dRopt = so3optimization(Rii, Rci, N)
    maxIter = 20;    
    % Solve for pose using our algorithm
    Rest = eye(3);
    for i=1:maxIter      % Gauss-Newton iterations
        LHS = zeros(3);
        RHS = zeros(3,1);
        
        for k=1:N
            G = JSO3(Rii{k},Rest);
            e = eSO3(Rci{k},Rii{k},Rest);
            LHS = LHS + G'*G;
            RHS = RHS + G'*e;
        end
        xi = -LHS \ RHS;
        Rest = vec2rot( xi ) * Rest;
    end
    dRopt = Rest;
end

function e = eSO3(Rc,Ri,Rest)
    e1 = rot2vec(Rc');
    e2 = rot2vec(Ri);
    e = e1 + Rest * e2;
end

function G = JSO3(Ri,Rest)
    so31 = rot2vec(Ri);
    so32 = Rest*so31;
    invJ = vec2jacInv( so32 );
    C = vec2rot(so32);
    G = invJ * (eye(3) - C);
end