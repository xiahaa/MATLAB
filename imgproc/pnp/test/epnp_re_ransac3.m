function [R, t] = epnp_re_ransac3(p,q)

    sol_iter = 1; %indicates if the initial solution must be optimized
    dims = 4;     %kernel dimensions
    
    minerror = 0.02;
    
    % 
    cw = [1 0 0;0 1 0;0 0 1;0 0 0]';
    Cw=[1 0 0 0;
        0 1 0 0;
        0 0 1 0;
        1 1 1 1];
    [M, Cw, Alph] = PrepareData1(p,q,Cw);

    % implementation of two point ransac for pnp
    minset = 6;
    maxiter = 3000;
    iter = 0;
    bestcost = 0;
    inliers = [];
    
    % total number
    n = size(p,2);
    pd = 0.9;         % Desired probability of choosing at least one sample
                       % free from outliers (probably should be a parameter)
    while iter < maxiter
        % random pick six
        ss = randperm(n, minset);
        id = [(ss-1).*2+1 ss.*2];
        N = M(id,:);
        [~,~,v] = svd(N);
        error21    = M(1:2:end,:) * v(:,end);
        error22    = M(2:2:end,:) * v(:,end);
        error2     = sqrt(error21.^2 + error22.^2);
        
        inlier = error2 < minerror;
        
        ninliers = sum(inlier);
                
        if ninliers > bestcost
            bestcost = ninliers;
            inliers = inlier;
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
    [M, ~, ~] = PrepareData1(p(:,inliers),q(:,inliers),Cw);
    [~,~,v] = eig(M'*M);
    Km = v(:,dims:-1:1);
    [R, t, ~] = KernelPnP1(cw, Km, dims, sol_iter);
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
