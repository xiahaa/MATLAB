function [R,t,aM] = soft_posit(p, q, params, K, varargin)%
%     p = [ -0.5000   -0.5000   -0.5000; 
%          0.5000   -0.5000   -0.5000;
%          0.5000    0.5000   -0.5000;
%         -0.5000    0.5000   -0.5000;
%         -0.5000   -0.5000    0.5000;
%          0.5000   -0.5000    0.5000;
%          0.5000    0.5000    0.5000;
%         -0.5000    0.5000    0.5000]';
% 
% 
%     q = [ 172.3829  -15.4229;
%         174.9147 -183.8248;
%         -28.3942 -147.8052;
%         243.2142  105.4463;
%         252.6934  -72.3310;
%          25.7430  -28.9218;
%          35.9377  149.1948]';
% 
%      params.beta0 = 2.0e-04;
%      params.noise_std = 0;
%      
%      K = [1500 0 0;0 1500 0;0 0 1];
    addpath ./3rdparty/posit/softPosit
    
    M = size(p,2);
    N = size(q,2);
    if size(q,1) == 2
        q = [q;ones(1,N)];
    end
    if size(p,1) == 3
        p = [p;ones(1,M)];
    end
    % normalized
    q = K\q;
    u = q(1,:);v = q(2,:);% 1xN
    
    focal_length = K(1,1);
    
    % R,t
    if nargin == 5
        R = varargin{1};
    else
        R = eye(3);
    end
    if nargin == 6
        t = varargin{2};
    else
        t = [0;0;50];
    end
    
%     R = [ 0.9149    0.1910   -0.3558;
%           -0.2254    0.9726   -0.0577;
%            0.3350    0.1330    0.9328];
%     t = [0; 0; 50];
    
    wk = [R(3,:)/t(3) 1]*p;% 1x4 * 4xM = 1xM
    
    % s = f/Tz = 1/Tz if coordiante has been normalized
    r1T = [R(1,:)/t(3), t(1)/t(3)];%1x4
    r2T = [R(2,:)/t(3), t(2)/t(3)];%1x4
    
    % fetch parameters
    beta_count = 0;
    pose_converged = 0;
    assign_converged = 0;
    found_pose = 0;
    beta = params.beta0;
    
    noise_std = params.noise_std;
    alpha = 9.21*noise_std^2 + 1;        % maximum pixel distance for considering match
    max_delta = sqrt(alpha)/2;          % Max allowed error per world point.
    beta_final = 0.5;                   % Terminate iteration when beta == betaFinal.
    beta_update = 1.05;                 % Update rate on beta.
    epsilon0 = 0.01;                    % Used to initialize assignement matrix.
    
    min_beta_count = 20;
    
    scale = 1/(max(M,N)+1);             % add one column for the slack variable
    
    aM = ones(N+1,M+1) + epsilon0;      % N+1 x M+1
    
    %% precompute
    pp = zeros(M, 16);
    for k = 1:M
        tmp = p(:,k) * p(:,k)';
        pp(k,:) = tmp(:);
    end
    ujpk = kron(u, p);
    vjpk = kron(v, p);
    
    % Deterministic annealing loop.
    while beta < beta_final && ~assign_converged
        % alternating optimization
        
        %% 1st step, assume R,t is known, update the assign matrix
        pu = r1T * p;% 1xM
        pv = r2T * p;
        
        pu = repmat(pu, N, 1);% NxM
        pv = repmat(pv, N, 1);% NxM
        
        ccu = u'*wk;% Nx1 1*M = N*M
        ccv = v'*wk;% Nx1 1*M = N*M, corrected pixel coordinates
        
        dist_mat = focal_length*focal_length*((pu - ccu).^2 + (pv - ccv).^2);% from real distance to pixel metrics
        
        % update assign matrix
        aM(1:N,1:M) = scale*exp(-beta*(dist_mat - alpha));% NxM
        aM(N+1,:) = scale; 
        aM(:,M+1) = scale;
        
        aM = sinkhornImp(aM);    % normalization assign matrix
        %% finish fist step
        num_matches = num_of_matches(aM);
        
        A1 = zeros(4,4);
        % 2 step, use updated assign matrix for pose estimation using POSIT
        sumcol = sum(aM(1:N,1:M), 1);% 1xM
       
        % vectorization
        sumcol = repmat(sumcol', 1, 16);
        A1 = sum(sumcol.*pp, 1);A1 = reshape(A1,4,4);
        if cond(A1) > 1e10
            disp('sumSkSkT is ill-conditioned, termininating search.');
            return;
        end
        
        L = inv(A1);                                    % Inv(L), a 4x4 matrix.

        pose_converged = 0;                              % Initialize for POSIT loop.
        count = 0;

        wkrep = repmat(wk, 1, N);
        amvec = vec(aM(1:N,1:M)')';
        weightedUi = sum(amvec.*wkrep.*ujpk,2);
        weightedVi = sum(amvec.*wkrep.*vjpk,2);
        
        % Compute the pose vectors. M = s(R1,Tx) and N = s(R2,Ty) where the
        % scale factor is s = f/Tz, and f = 1.  These are column 4-vectors.
        r1T= L * weightedUi;               % M
        r2T = L * weightedVi;              % N
        
        % Compute the rotation matrix and translation vector corresponding to
        % the computed pose vectors.
        if 0  % Chang & Tsai calculation of R and T.
            [U, S, V] = svd([r1T(1:3)'; r2T(1:3)']');
            A = U * [1 0; 0 1; 0 0] * V';
            r1 = A(:,1)';
            r2 = A(:,2)';
            r3 = cross(r1,r2);
            tz = 2 / (S(1,1) + S(2,2));
            tx = r1T(4) * tz;
            ty = r2T(4) * tz;
            r3T= [r3 tz];
        else 
            % Standard calculation of R and T.  The rotation matrix may not be
            % orthonormal.  The object must distort when the rotation matrix
            % is not orthonormal.
            r1TSquared = r1T(1)*r1T(1) + r1T(2)*r1T(2)+ r1T(3)*r1T(3);
            r2TSquared = r2T(1)*r2T(1) + r2T(2)*r2T(2)+ r2T(3)*r2T(3);
            tz = sqrt(2/(r1TSquared+r2TSquared));   % Chang & Tsai's suggestion.
            r1N = r1T*tz;                   % Column 4-vectors: (R1,Tx).
            r2N = r2T*tz;                   %                   (R2,Ty).
            r1 = r1N(1:3)';                  % Three rows of the rotation matrix.
            r2 = r2N(1:3)';
            r3 = cross(r1,r2);
            r3T= [r3 tz];                  % Column 4-vector: (R3,Tz).
            tx = r1N(4);
            ty = r2N(4);
        end 
        r1T = [r1, tx]/tz;
        r2T = [r2, ty]/tz;
        
        wk = [r3/tz 1]*p;% 1x4 * 4xM = 1xM
        
        delta = sqrt(sum(sum(aM(1:N,1:M) .* dist_mat))/M);
        pose_converged = delta < max_delta;
        % end % of POSIT loop
        
        % Update the "annealing temperature" and determine if the assignments
        % have converged.
        beta = beta_update * beta;
        beta_count = beta_count + 1;   % Number of deterministic annealing iterations.
        assign_converged = pose_converged & beta_count > min_beta_count;

        % Form the best estimates for the translation vector and rotation matrix.
        t = [tx; ty; tz];
        R = [r1;r2;r3];

        % Has the pose converged?
        found_pose = (delta < max_delta & beta_count > min_beta_count);
    end
end

function num = num_of_matches(M)
    num = 0;
    
    [m,n] = size(M);
    
    for i = 1:m-1
        % for each feature
        [vmax, maxid] = max(M(i,:));
        if (maxid == n) continue; end % slack is the maximum means no match
        if all(vmax > M([1:maxid-1, maxid+1:m],maxid))
            % also the maximum along the column direction, then a correct
            % match
            num = num + 1;
        end
    end
end