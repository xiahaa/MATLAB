function [R, t] = epnp_re_ransac2(p,q)

    sol_iter = 1; %indicates if the initial solution must be optimized
    dims = 4;     %kernel dimensions
    
    minerror = 0.01;
    
%     qn1 = p(:,1) ./ norm(p(:,1));qn1 = qn1';
%     qn2 = [-qn1(2) qn1(1) 0];qn2 = qn2./norm(qn2);
%     qn3 = cross(qn1,qn2);qn3 = qn3 ./ norm(qn3);
%     qn1 = [1 0 0];
%     qn2 = [-qn1(2) qn1(1) 0];qn2 = qn2./norm(qn2);
%     qn3 = cross(qn1,qn2);qn3 = qn3 ./ norm(qn3);
    
    % implementation of two point ransac for pnp
    minset = 1;
    maxiter = 1e6;
    iter = 0;
    bestcost = 0;
    inliers = [];
    
    % total number
    n = size(p,2);
    id = 1:n;
    pd = 0.99;         % Desired probability of choosing at least one sample
                      % free from outliers (probably should be a parameter)
    k = 1;
    while iter < maxiter
        if k > n
            break;
        end
        % random pick two
        ss = id(k);%randperm(n, minset);
        k = k + 1;
        
        qn1 = p(:,ss) ./ norm(p(:,ss));qn1 = qn1';
        qn2 = [-qn1(2) qn1(1) 0];qn2 = qn2./norm(qn2);
        qn3 = cross(qn1,qn2);qn3 = qn3 ./ norm(qn3);
        
        %
        cw = [qn1;qn2;qn3;0 0 0]';
        Cw = [qn1 0;qn2 0;qn3 0;1 1 1 1];
    
        [M, Cw, Alph] = PrepareData1(p,q,Cw);

        %roubst kernel estimation
        [Km, cidinliers, robustiters] = my_robust_kernel_noise2(M,dims,minerror);
        

        
        %% rely on reporjection error, this has to be run together with a softweight or outlier suppression
%         Km=kernel_noise(M,dims); %Compute kernel M
%         [R, t] = KernelPnP1(cw, Km, dims, 0);
        % reproj_err
%         rep_err = calc_reproj_err(p,q,R,t,eye(3));
%         cidinliers = id(rep_err < 0.2);
        ninliers = numel(cidinliers);
        
        if ninliers > bestcost
            bestcost = ninliers;%sum(rep_err(cinliers));
            inliers = cidinliers;
%             ninliers = numel(rinlier); 
            % Update estimate of N, the number of trials to ensure we pick,
            % with probability p, a data set with no outliers.
            fracinliers =  ninliers/n;
            pNoOutliers = 1 -  fracinliers^minset;
            pNoOutliers = max(eps, pNoOutliers);  % Avoid division by -Inf
            pNoOutliers = min(1-eps, pNoOutliers);% Avoid division by 0.
            maxiter = log(1-pd)/log(pNoOutliers);
        end
        iter = iter + 1;
    end
    
    % refine
    qn1 = [1 0 0];
    qn2 = [-qn1(2) qn1(1) 0];qn2 = qn2./norm(qn2);
    qn3 = cross(qn1,qn2);qn3 = qn3 ./ norm(qn3);
    cw = [qn1;qn2;qn3;0 0 0]';
    Cw = [qn1 0;qn2 0;qn3 0;1 1 1 1];
    [M, Cw, Alph] = PrepareData1(p(:,inliers),q(:,inliers),Cw);
    Km=kernel_noise(M,dims); %Compute kernel M    
    [R, t, err] = KernelPnP1(cw, Km, dims, sol_iter);
    
end

function err = calc_reproj_err(p,q,R,t,K)
    pc = R*p+repmat(t,1,size(p,2));
    pc = pc./pc(3,:);
    pq = K*pc;
    err = pq(1:2,:) - q(1:2,:);%2xn
    err = sqrt(diag(err'*err));
end

function [M,Cw, Alph] = PrepareData1(Pts,impts,Cw)
    U=impts;
    
    %compute alphas (linear combination of the control points to represent the 3d points)
%     Alph=compute_alphas(Xw,Cw');
    Pts = [Pts;ones(1,size(Pts,2))];
    Alph = Cw\Pts;
    Alph = Alph';
    %Compute M
    M=ComputeM1(U(:),Alph);
end

function M = ComputeM1(U,Alph)
    %ATTENTION U must be multiplied by K previously
    M = kron(Alph,[1 0 -1; 0 1 -1]);
    M(:,[[3,6,9,12]]) =  M(:,[3,6,9,12]) .* (U * ones(1,4));
end

function [K, idinliers, i]=my_robust_kernel_noise2(M,dimker, minerror)
    m   = size(M,1);
    prev_cost = Inf;
    maxIter = 30;
    softWeight = ones(round(m * 0.5),1);%softWeight = softWeight./sum(softWeight(:));
    prev_inlier_cnt = -m;
    id1 = 1:2:m-1;
    id2 = 2:2:m;
    for i=1:maxIter
        N(id1,:) = softWeight.*M(id1,:);
        N(id2,:) = softWeight.*M(id2,:);
        [~,S,V] = svd(N);
        id = find(diag(S)>0);
        v = V(:,id(end));
        
        error21    = M(1:2:end,:) * v(:,end);
        error22    = M(2:2:end,:) * v(:,end);
        error2     = sqrt(error21.^2 + error22.^2);
            
        w = minerror./error2;
%         softWeight(id1,1) = softWeight(id1,1).*w;
%         softWeight(id2,1) = softWeight(id2,1).*w;
        softWeight = softWeight.*w;
        %softWeight = softWeight./sum(softWeight(:));
        softWeight(softWeight>1) = 1;
        
        ninliers = sum(error2<minerror);
       
        ccost = sum(error2);
        
        if (abs(ccost-prev_cost) < 1e-6) || (ninliers-prev_inlier_cnt) < 1
            break;
        else
            prev_inlier_cnt = ninliers;
            prev_cost = ccost;
            resv    = v;
        end
    end
    idx = 1:m/2;
    K = resv(:,end-dimker+1:end);   
    idinliers = idx(error2<minerror);
end

function [K, idinliers, i]=my_robust_kernel_noise1(M,dimker, minerror)
    m   = size(M,1);
    id  =round(m/8);
    idx = 1:m;
    prev_sv = Inf;
    pairs = 1; %each correspondence is a couple of equations
    for i=1:30
        N = M(idx,:);
        [~,~,v] = svd(N'*N);
       
        if (pairs)
            error21    = M(1:2:end,:) * v(:,end);
            error22    = M(2:2:end,:) * v(:,end);
            error2     = sqrt(error21.^2 + error22.^2);
            
            [sv, tidx] = sort(error2);        

            med = sv(floor(m/16)); 

        else
            error2    = M * v(:,end);
            [sv, tidx] = sort(error2.^2);
            med = sv(floor(m/2)); 
        end
     
        ninliers = sum(sv<max(med,minerror));

        if (med >= prev_sv)
            break;
        else
            prev_sv = med;
            resv    = v;
            if(pairs)
                residx  = tidx(1:ninliers);
            else
                %always pairs = 1!! :P
               
            end
        end
        
        if(pairs)
            tidx2     = tidx'*2;
            ttidx     = [tidx2-1; tidx2];
            tidx2     = ttidx(:);
            idx       = tidx2(1:2*ninliers);
        else
            idx       = tidx(1:ninliers);
        end
    end
    
    K = resv(:,end-dimker+1:end);   
    idinliers = residx;
end

function [R,T, err] = KernelPnP1(Cw, Km, dims, sol_iter)

    vK = reshape(Km(:,end),3,dims);
    
    %precomputations
    X.P     = Cw;
    X.mP    = mean(X.P,2);
    X.cP    = X.P - X.mP * ones(1,dims);
    X.norm  = norm(X.cP(:));
    X.nP    = X.cP/X.norm;
    
    %procrustes solution for the first kernel vector
    if (mean(vK(3,:)<0))
        vK = -vK;
    end
    [R,b,mc] = myProcrustes1(X,vK);
    
    solV  = b * vK;
    solR  = R;
    solmc = mc;
  
    % procrustes solution using 4 kernel eigenvectors
    if sol_iter
         err = Inf;
         n_iterations=10;%10
         for iter=1:n_iterations
             % projection of previous solution into the null space
             A = R * (- mc +X.P);
             abcd = Km \ A(:);
             newV = reshape(Km * abcd,3,dims);
             
             %eucliedean error
             newerr = norm(R' * newV + mc - X.P);
             
             if ((newerr > err) && (iter>2))
                 break;
             else
                 %procrustes solution
                 [R,b,mc] = myProcrustes1(X,newV);
                 solV = b * newV;
                 
                 solmc = mc;
                 solR = R;
                 err = newerr;
             end
             
         end
    end
       
    R  = solR;
    mV = mean(solV,2);
     
    T = mV - R * X.mP;
end


function [R, b, mc] = myProcrustes1(X,Y)
%X is an structure containing points, points centered in the origin, points
%normalized
%Y are 3D points
    dims = size(Y,2);
    mY = mean(Y,2);
    cY = Y - mY * ones(1,dims);
    ncY = norm(cY(:));
    tcY = cY/ncY;
    
    A = X.nP * tcY';
    [L, D, M] = svd(A);
  
%     R = M * L';
%     
%     if(mY(3)>0 && det(R)<0)
        R = M * diag([1,1,sign(det(M*L'))])* L';
%   end
    
    b = sum(diag(D)) * X.norm/ncY;
    c = X.mP - b*R'*mY;
    mc = c * ones(1,dims);
end
