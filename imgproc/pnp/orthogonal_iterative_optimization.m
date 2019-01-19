function varargout = orthogonal_iterative_optimization(varargin)
% Implementation of the iterative optimization algorithm proposed in
% Lu C P, Hager G D, Mjolsness E. 
% Fast and globally convergent pose estimation from video images[J]. 
% IEEE Transactions on Pattern Analysis & Machine Intelligence, 
% 2000 (6): 610-622.
%
% Author: xiahaa@space.dtu.dk
    pi = varargin{1};%% 3D points,  3xN
    vi = varargin{2};%% normalized 2D points,  2xN
    N = size(pi,2);
    
    %% precompute some constant variables
    if size(vi,1) == 2
        vih = [vi;ones(1,N)];%% to homogeneous
    else
        vih = vi;%% assume already homogeneous
    end
    
    %% to center
    pmean = mean(pi,2);
    pic = pi - repmat(pmean,1,N);
    
    % projection matrix
    Fi = zeros(3,3,N);
    Vi = zeros(3,3,N);
    C1 = zeros(3,3);
    for i = 1:N
        Vi(:,:,i) = vih(:,i) * vih(:,i)' ./ (vih(:,i)'*vih(:,i));
        C1 = C1 + Vi(:,:,i);
        Fi(:,:,i) = Vi(:,:,i) - eye(3);
    end
    C1 = eye(3) - C1 ./ N;
    C2 = inv(C1)./N;
    
    %% iteration exit conditions
    cond1 = 1e-10;
    cond2 = 1e-10;
    maxiter = 100;
    if nargin == 5
        cond1 = varargin{5};
        cond2 = varargin{5};
    elseif nargin == 6
        cond1 = varargin{5};
        cond2 = varargin{6};
    end
    %% maximum iteration
    if nargin == 7
        maxiter = varargin{7};
    end
    
    error1 = 1e6;
    error2 = 1e6;
    
    if nargin > 2
        R0 = varargin{3};
        t0 = varargin{4};
    else
        %% initialize with weak projection
        [R0] = svd_3d23d(pic, vih);
    end
    
    Rk = R0;
    tk = optimize_t(Fi, pic, Rk, C2);
    f0 = evalf(Rk, pic, tk, Fi);
    iter = 1;
    
    while (error1 > cond1) && (error2 > cond2) && iter < maxiter
        % current hypothesis
        q = Rk * pic + repmat(tk, 1, N);
        qhat = prepareq(Vi, q);
        [Rnew] = svd_3d23d(pic, qhat);
        tnew = optimize_t(Fi, pic, Rnew, C2);%% estimate t
        fnew = evalf(Rnew, pic, tnew, Fi);
        
        Rk = Rnew;
        tk = tnew;
        
        %% new error
        error1 = norm(Rk-R0,'fro');
        error2 = abs(fnew - f0);
        R0 = Rk;
        t0 = tk;
        f0 = fnew;
        iter = iter + 1;
    end
    Ropt = Rk;
    topt = tk;
    objerr = fnew;
    objerr = sqrt(objerr) / N;
    if nargout >= 5
        qproj = Ropt*pic + repmat(topt,1,N);
        qproj = qproj ./ qproj(3,:);
        imgerr = qproj - vi;
        sumerr = sum(diag(imgerr'*imgerr));
        imgerravg = sqrt(sumerr) / N;
        varargout{5} = imgerravg;
    end
    topt = topt - Ropt * pmean;
    
    varargout{1} = Ropt;
    varargout{2} = topt;
    varargout{3} = iter;
    varargout{4} = objerr;
end

function qhat = prepareq(V, q)
    qhat = q;
    for i = 1:size(q,2)
        qhat(:,i) = V(:,:,i) * q(:,i);
    end
end

function f = evalf(R, p, t, Fi)
    q = R*p + repmat(t,1,size(p,2));
    f = 0;
    for i = 1:size(p,2)
        err = -Fi(:,:,i) * q(:,i);
        f = f + err'*err;
    end
end

function topt = optimize_t(Fi, p, R, C2)
    ps = R * p;
    s = zeros(3,1);
    for i = 1:size(Fi,3)
        s = s + Fi(:,:,i) * ps(:,i); 
    end
    topt = C2 * s;
end

function [Ropt] = svd_3d23d(ptsrc, ptdst)
    %% SVD solution for estimating motion from two 3D point sets
    ptsrcmean = mean(ptsrc,2);
    ptdstmean = mean(ptdst,2);

    ptsrcrefine = ptsrc - repmat(ptsrcmean, 1, size(ptsrc,2));
    ptdstrefine = ptdst - repmat(ptdstmean, 1, size(ptsrc,2));

    Y = ptdstrefine';
    X = ptsrcrefine;
    S = X*Y;
    [U,~,V] = svd(S);

    D = V*U';
    if det(D) < 0
        Ropt = V*[1 0 0;0 1 0;0 0 -1]*U';
    else
        Ropt = D;
    end
%     topt = ptdstmean - Ropt * ptsrcmean;
end