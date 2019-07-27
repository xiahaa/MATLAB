function [ X, MeanA, MeanB, t_error ] = batchSolveSoftUseScrew(A, B)
    %% Mixed version for solving AX = XB
    %%
    
    %% first step, use screw parameters to select similar set
    addpath('C:/Users/xiahaa/Documents/MATLAB/hand_eye_calib/solver');

    [id1,id2] = findMaximumConsistentSet(A,B);
    A = A(:,:,id1(:,1));
    B = B(:,:,id2(:,1));
    
    [a1,a2,a3]  = size(A);
    [b1,b2,b3]  = size(B);
%     A_mex = reshape(A, a1, a2*a3);
%     B_mex = reshape(B, b1, b2*b3);
%     [a_ins, b_ins] = identifyInliers(A_mex, B_mex);
%     Ain = A(:,:,a_ins);
%     Bin = B(:,:,b_ins);
%     [a1,a2,a3]  = size(Ain);
%     [b1,b2,b3]  = size(Bin);
%     A_mex = reshape(Ain, a1, a2*a3);
%     B_mex = reshape(Bin, b1, b2*b3);

    A_mex = reshape(A, a1, a2*a3);
    B_mex = reshape(B, b1, b2*b3);
%     
    w1 = ones(1,a3);
    w2 = ones(1,b3);
%     
    MeanA = mean_Taylor_2nd_adv_recursive3( A_mex,  w1, 50); 
    MeanB = mean_Taylor_2nd_adv_recursive3( B_mex,  w2, 50); 
%     
    [X, t_error] = findX(A_mex, B_mex, MeanA, MeanB, w1, w2);
    
%     Ain = A;
%     Bin = B;
%     T2 = zeros(size(Ain,3),4,4);
%     for i = 1:size(Ain,3)
%         T2(i,:,:) = Ain(:,:,i);
%     end
%     T1 = zeros(size(Bin,3),4,4);
%     for i = 1:size(Bin,3)
%         T1(i,:,:) = Bin(:,:,i);
%     end
%     X = sol_manifold_opt_SE3(T1,T2,size(T1,1));%sol_cvx_sdp
%     MeanA = [];
%     MeanB = [];
%     t_error = 0;
end

function [id1,id2] = findMaximumConsistentSet(A,B)
    n1 = size(A,3);
    n2 = size(B,3);
    
    % select the pivoting set as the compared set
    if n1 > n2
        query = B;
        n3 = n2;
        train = A;
        n4 = n1;
    else
        query = A;
        n3 = n1;
        train = B;
        n4 = n2;
    end

    % collect screw parameters
    screws_query = zeros(n3,2);
    screws_train = zeros(n4,2);
    valid_query = zeros(n3,1);
    valid_train = zeros(n4,1);
    for i = 1:n3
        T1 = query(:,:,i);
        [thetaa, da, la, ma] = screwParams(T1(1:3,1:3),T1(1:3,4));
        
        % in case of NAN or INF
        isvalid = (ma(1) == ma(1)) & (ma(2) == ma(2)) & (ma(3) == ma(3)) ...
                 & (la(1) == la(1)) & (la(2) == la(2)) & (la(3) == la(3));

        if isvalid
            valid_query(i) = 1;
            screws_query(i,:) = [thetaa, da];            
        end
    end
    for i = 1:n4
        T1 = train(:,:,i);
        [thetaa, da, la, ma] = screwParams(T1(1:3,1:3),T1(1:3,4));
        
        % in case of NAN or INF
        isvalid = (ma(1) == ma(1)) & (ma(2) == ma(2)) & (ma(3) == ma(3)) ...
                 & (la(1) == la(1)) & (la(2) == la(2)) & (la(3) == la(3));

        if isvalid
            valid_train(i) = 1;
            screws_train(i,:) = [thetaa, da];            
        end
    end
    
    % 3rd step, similarities matrix
    similarities = zeros(n3,n4);
    w1 = 1;w2 = 1;
    scores = zeros(n3,2);
    matches = zeros(n3,2);
    for i = 1:n3
        thetaa = screws_query(i,1);
        da = screws_query(i,2);
        tmp = repmat([thetaa, da], n4, 1);
        tmp1 = abs(tmp - screws_train)';
        similarity = (w1.*tmp1(1,:) + w2.*tmp1(2,:));
        similarities(i,:) = similarity;
        [val,sid] = sort(similarity,'ascend');
        scores(i,:) = [val(1),val(2)];
        matches(i,:) = [sid(1),sid(2)];
    end
    
    % 4th step, find score for each query
    bestscore = scores(:,1);
    inlierthresh = 0.5;%mean(bestscore);%0.01;
    valid = bestscore < inlierthresh;
    
    % otherway, just voting
    validid1 = zeros(1,n3);
    validid2 = zeros(1,n4);
    for i = 1:n3
        similarity = similarities(i,:);
        [val,sid] = sort(similarity,'ascend');
        vv = val <= inlierthresh;
        sid = sid(vv);
        if ~isempty(sid)
            validid1(i) = 1;
            validid2(sid) = 1;
        end
    end
    validid1 = validid1 == 1;
    validid2 = validid2 == 1;
    
    id1 = (1:n3)';
    id2 = (1:n4)';
    id1 = id1(validid1);
    id2 = id2(validid2);
    
    if n1 > n2
        tmp = id1;
        id1 = id2;
        id2 = tmp;
    else
%         id = [id1 id2];
    end
    
end

function Sigma = cov_SE3_weight(MX, X, weight)
    % compute the covariance of a batch of matrices
    n = size(X,2) / 4;
    Sigma = zeros(6,6);
    
    % do first order approximation
    for i = 1:n
        X_i = X(:,4*(i-1)+1:4*i);
        P = MX\X_i;
        Sigma = Sigma + weight(i).*se3_vec(logm(P))*se3_vec(logm(P))';
    end
    Sigma = Sigma./sum(weight);
end

function [X, t_error] = findX(A_mex, B_mex, MeanA, MeanB, w1, w2)
    % this will serve as the back-end optimization, up-to-this stage,
    % inliers should be identified.
    SigA = cov_SE3_weight(MeanA, A_mex, w1);
    SigB = cov_SE3_weight(MeanB, B_mex, w2);
    
    [ VA, ~ ] = eig( SigA(1:3,1:3) );
    [ VB, ~ ] = eig( SigB(1:3,1:3) );
    Q1 = eye(3);
    Q2 = [-1 0 0; 0 -1 0; 0 0 1];
    Q3 = [-1 0 0; 0 1 0; 0 0 -1];
    Q4 = [1 0 0; 0 -1 0; 0 0 -1];
    
    % 8 solution
    Rx_solved(:,:,1) = VA*Q1*VB';
    Rx_solved(:,:,2) = VA*Q2*VB';
    Rx_solved(:,:,3) = VA*Q3*VB';
    Rx_solved(:,:,4) = VA*Q4*VB';
    Rx_solved(:,:,5) = VA*-Q1*VB';
    Rx_solved(:,:,6) = VA*-Q2*VB';
    Rx_solved(:,:,7) = VA*-Q3*VB';
    Rx_solved(:,:,8) = VA*-Q4*VB';
    
    % in fact, just to extract the screw params
    [~, Na, ~, ~] = param_extract(MeanA);
    [~, Nb, ~, ~] = param_extract(MeanB);
    % rotation axis
    na = so3_vec(Na);
    nb = so3_vec(Nb);
    % traverse the 8 sols and find the minimum
    min = inf;
    Rx = Rx_solved(:,:,1);
    for i = 1:8
        if (abs(det(Rx_solved(:,:,i))-1)<0.001) && (norm(na-Rx_solved(1:3,1:3,i)*nb) < min)
            min = norm(na-Rx_solved(1:3,1:3,i)*nb);
            Rx = Rx_solved(:,:,i);
        else
        end
    end
    % equ 20b, recorver translation
    tx_temp = so3_vec(((Rx'*SigA(1:3,1:3)*Rx)^-1*(SigB(1:3,4:6)-Rx'*SigA(1:3,4:6)*Rx))');
    tx = -Rx*tx_temp;
    X = [Rx tx; [0 0 0] 1];
    t_error = (MeanA(1:3,1:3) - eye(3))*tx - Rx*MeanB(1:3,4) + MeanA(1:3,4);
    t_error = norm(t_error);
end

function [a_ins, b_ins] = identifyInliers(A_mex, B_mex)
    % number of sampels for each SE3
    n1 = round(size(A_mex,2)/4);
    n2 = round(size(B_mex,2)/4);

    w1 = ones(1,n1);
    w2 = ones(1,n2);
    
    maxiter = 30;
    iter = 1;
    inlierthresh = 0.8;
    
    Xslack = eye(4);
    Ts = zeros(4,4,2);
    while iter <= maxiter
        MeanA = mean_Taylor_2nd_adv_recursive3(A_mex, w1, 5);
        MeanB = mean_Taylor_2nd_adv_recursive3(B_mex, w2, 5);    
        % update weight
        [X, ~] = findX(A_mex, B_mex, MeanA, MeanB, w1, w2);
        
        oldw1 = w1;
        oldw2 = w2;
        
        Ts(:,:,1) = Xslack;
        err1 = calcErr(A_mex,B_mex,Xslack);err1 = sum(min(err1));
        reprojerrors(1) = err1;
        Ts(:,:,2) = X;
        err2 = calcErr(A_mex,B_mex,X);err2 = sum(min(err2));
        reprojerrors(2) = err2;
        [Tfused, ~, ~] = poseFusionOnSE3(Ts, reprojerrors, 2);
        Xslack = Tfused;
        
        [w1,w2,dists] = update_weight(A_mex,B_mex,Xslack,inlierthresh);
        
%         w1 = 0.8*oldw1+neww1*0.2;
%         w2 = 0.8*oldw2+neww2*0.2;
%         w1(w1>1) = 1;
%         w2(w2>1) = 1;
%         
        if iter > 1
            if norm(oldw1-w1) < 1e-6 && norm(oldw2-w2) < 1e-6
                break;
            end
        end
        iter = iter + 1;
    end    
    
    % identify inliers according to their weights, here add a function to
    % find the consistent set
    valid1 = zeros(1,n1);
    valid2 = zeros(1,n2);
    for i = 1:n1
        if w1(i) < 0.2  continue; end
        [~,minj] = min(dists(i,:));
        [~,mini] = min(dists(:,minj));
        if i == mini
            valid1(i) = 1; valid2(minj) = 1;
        end
    end
    
    id1 = 1:n1;
    id2 = 1:n2;
    
    a_ins = id1(valid1==1);
    b_ins = id2(valid2==1);
end

function dists = calcErr(A_mex,B_mex,X)
    m = round(size(A_mex,2)/4);
    n = round(size(B_mex,2)/4);
    
    se3a = zeros(6,m);
    se3b = zeros(6,n);
    for i = 1:m
        se3a(:,i) = tran2vec(A_mex(:,(i-1)*4+1:i*4)*X);
    end
    
    for i = 1:n
        se3b(:,i) = tran2vec(X*B_mex(:,(i-1)*4+1:i*4));
    end
    
    % compute inter-distance
    dists = zeros(m,n);
    for i = 1:m
        err = se3b - repmat(se3a(:,i),1,n);
        err = sqrt(diag(err'*err));
        dists(i,:) = err';
    end
end

function [Test, V, Sigma_est] = poseFusionOnSE3(Ts, reprojerrors, N)
    %% 1st, prepare covariance matrix
    Sigma = cell(N,1);
    for j = 1:N
        Sigma{j} = eye(6,6).*(reprojerrors(1,j)*reprojerrors(1,j));
        T{j} = Ts(:,:,j);
    end
    
    % precompute inverses and Cholesky factors
    for k = 1:N
        invSigma{k} = inv( Sigma{k} );
        cholSigma{k} = chol( Sigma{k}, 'lower' );
    end
    
    n = 7;
    maxIter = 20;
    
    % Solve for pose using our algorithm
    Test = eye(4);
    for i=1:maxIter      % Gauss-Newton iterations
        LHS = zeros(6);
        RHS = zeros(6,1);
        for k=1:N
            xik = tran2vec( Test*inv(T{k}) );
            if n <= 6
               invJ = vec2jacInvSeries( xik, n );
            else
               invJ = vec2jacInv( xik );
            end
            invJtS = invJ'*invSigma{k};
            LHS = LHS + invJtS*invJ;
            RHS = RHS + invJtS*xik;
        end
        xi = -LHS \ RHS;
        Test = vec2tran( xi )*Test;
    end

    % How low did the objective function get?
    V = 0;
    for k=1:N
        xik = tran2vec( Test*inv(T{k}) );
        V = V + xik'*invSigma{k}*xik/2;
    end

    % How big is the covariance?
    Sigma_est = inv(LHS);
end

function [w1,w2,dists] = update_weight(A_mex,B_mex,X,inlierthresh)
    m = round(size(A_mex,2)/4);
    n = round(size(B_mex,2)/4);
    dists = calcErr(A_mex,B_mex,X);
    
    % update strategy
    min1 = min(dists,[],2)';
    min2 = min(dists,[],1);
    gamma = 0.0001;
    min1(abs(min1) < 1e-6) = gamma;
    min2(abs(min2) < 1e-6) = gamma;
    
    w1 = inlierthresh ./ min1;
    w2 = inlierthresh ./ min2;
%     w1 = w1./sum(w1).*m;
%     w2 = w2./sum(w2).*n;
    w1(w1>1) = 1;
    w2(w2>1) = 1;
    
%     w1 = w1;
%     w2 = w2./n;
    
%     l1 = 1/m;
%     l2 = 1/n;
%     beta = 0.001;
%     alpha = 3*pi/180.0;
%     wp = l1.*exp(-beta.*(min1-alpha))';    
%     wq = l2.*exp(-beta.*(min2-alpha));
%     wp = wp ./ norm(wp);
%     wq = wq ./ norm(wq);
end