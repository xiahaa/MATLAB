function [R, t] = epnp_re_ransac4(p,q)
% this file trys to solve the pnp problem with outliers by placing the
% ransac in advance.
    if size(q,1) == 2
        q = [q;ones(1,size(q,2))];
        % to normalized vector
        qn = q ./ sqrt(q(1,:).^2+q(2,:).^2+q(3,:).^2);
    end

     % implementation of two point ransac for pnp
    minset = 1;
    maxiter = 1e6;
    iter = 0;
    bestcost = 0;
    inliers = [];
    
    % total number
    n = size(p,2);
    pd = 0.99;         % Desired probability of choosing at least one sample
                       % free from outliers (probably should be a parameter)
    
    inlierthreshold = 0.05;
    while iter < maxiter
        % choosing one point as the control point
        control_point_id = 1;%randperm(n, minset);
        
        control_point_3d = p(:,control_point_id);
        control_point_2d = qn(:,control_point_id);
        
        % generate pairwise constraints
        dummy_id = 1:n;
        dummy_id(control_point_id) = [];
        
        other_points_3d = p(:,dummy_id);
        other_points_2d = qn(:,dummy_id);
        
        v1 = other_points_3d - repmat(control_point_3d,1,n-1);
        v1 = v1 ./ sqrt(v1(1,:).^2+v1(2,:).^2+v1(3,:).^2);
        v1 = v1';
        v2 = cross(other_points_2d, repmat(control_point_2d,1,n-1));
        v2 = v2 ./ sqrt(v2(1,:).^2+v2(2,:).^2+v2(3,:).^2);
        v2 = v2';
        % kernel
%         M = kron(v2', v1');
        M = [v2(:,1).*v1 v2(:,2).*v1 v2(:,3).*v1];
        
        [eigv,eige] = eig(M'*M);
        rvec = eigv(:,1);
        Rbar = rvec([1 2 3;4 5 6;7 8 9]);
        [U,~,V] = svd(Rbar);
        D = V*U';
        if det(D) < 0
            Ropt = V*[1 0 0;0 1 0;0 0 -1]*U';
        else
            Ropt = V*U';
        end
        rvec = vec(Ropt');
        
        error = M*rvec;
        
        inlier = abs(error) < inlierthreshold;
        ninliers = sum(inlier);
                
        if ninliers > bestcost
            bestcost = ninliers;
            inliers = [control_point_id dummy_id(inlier)];
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
    
    sol_iter = 1; %indicates if the initial solution must be optimized
    dims = 4;     %kernel dimensions
    % 
    cw = [1 0 0;0 1 0;0 0 1;0 0 0]';
    Cw=[1 0 0 0;
        0 1 0 0;
        0 0 1 0;
        1 1 1 1];
    % refine
    [M, ~, ~] = PrepareData1(p(:,inliers),q(:,inliers),Cw);
    [~,~,v] = eig(M'*M);
    Km = v(:,dims:-1:1);
    [R, t, ~] = KernelPnP1(cw, Km, dims, sol_iter);
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