function varargout = sol_cvx_sdp(TA,TB,N)
%% solving hand eye calibration using SDP

%% Author: xiahaa@space.dtu.dk
    if N < 2
        error('At least two samples needed for unique solution!');
        varargout{1} = [];
        return;
    end
    
    format long;
        
    Asq = preprocessing_q(TA,TB,N);
    x0q = [1 0 0 0 1 0 0 0]';
    [Ask,x0k] = preprocessing_kron(TA,TB,N);
    
    [T2, time2] = solve_SQP(x0q, Asq);
    
%     [T4, time4] = solv_SQP_k(x0k,Ask, options);

    varargout{1} = T2;
%     varargout{2} = T2;
%     varargout{5} = time1;
%     varargout{6} = time2;
%     varargout{7} = time3;
%     varargout{8} = time4;
end

function [As,x0] = preprocessing_kron(TA,TB,N)
    dim = size(TA,2);
    C = zeros(N*12,13);
    As2 = zeros(13,13);
    for i = 1:N
        T1 = TA(i,:,:);T1 = reshape(T1,dim,dim,1);
        T2 = TB(i,:,:);T2 = reshape(T2,dim,dim,1);
        Ra = T2(1:3,1:3);ta = T2(1:3,4);
        Rb = T1(1:3,1:3);tb = T1(1:3,4);

        id = (i-1)*12;
        C(id+1:id+12,:) = [kron(eye(3),Ra)-kron(Rb',eye(3)) zeros(9,3) zeros(9,1); ...
                           kron(tb', eye(3)) eye(3)-Ra -ta];
       As2 = As2 + C(id+1:id+12,:)'*C(id+1:id+12,:);
    end
%         As = C'*C;
    As = (As2+As2')/2;
    %% vector to optimize
    %     x0 = C\d;
    x0 = [1 0 0 0 1 0 0 0 1 0 0 0 1]';
    
end

function [T, time] = solv_SQP_k(x0,As)
    x = x0;
    AA = 2*As;
    tic
    for k = 1:1000
        [L,S] = conslin2(x);
        soln = quadprog(AA, [], [],[],L,S,[],[],[]);
        xnew = soln(1:13);
        if norm(xnew-x) < 1e-12
            disp(['converge at step ', num2str(k)]);
            break;
        end
        x = xnew;
    end
    time = toc;
    disp(['SQP :',num2str(time)]);
    R12 = [x(1) x(4) x(7); ...
           x(2) x(5) x(8); ...
           x(3) x(6) x(9)];
    t12 = x(10:12);
    R12 = R12*inv(sqrtm(R12'*R12));
    T = [R12 t12;[0 0 0 1]];
end

function As = preprocessing_q(TA,TB,N)
    dim = size(TA,2);
    Nv = 0;
    ids = [];
    thetaas = zeros(N,1);das = zeros(N,1);
    thetabs = zeros(N,1);dbs = zeros(N,1);
    las = zeros(N,3);mas = zeros(N,3);
    lbs = zeros(N,3);mbs = zeros(N,3);
    for i = 1:N
        T1 = TA(i,:,:);T1 = reshape(T1,dim,dim,1);
        T2 = TB(i,:,:);T2 = reshape(T2,dim,dim,1);

        [thetaa, da, la, ma] = screwParams(T2(1:3,1:3),T2(1:3,4));
        [thetab, db, lb, mb] = screwParams(T1(1:3,1:3),T1(1:3,4));

        thetaas(i) = thetaa;das(i) = da;las(i,:) = la';mas(i,:) = ma';
        thetabs(i) = thetab;dbs(i) = db;lbs(i,:) = lb';mbs(i,:) = mb';

    %             if abs(thetaa-thetab) < 0.15 && abs(da-db) < 0.15
            Nv = Nv + 1;
            ids = [ids;i];
    %             end
    end
%         T = zeros(6*Nv, 8);%% formulation 1
    T = zeros(8, 8);
    As = zeros(8,8);
    for i = 1:Nv
        ii = ids(i);
        qas = [cos(thetaas(ii)*0.5);sin(thetaas(ii)*0.5).*las(ii,:)'];
        qdas = [(-das(ii)*0.5)*sin(thetaas(ii)*0.5);sin(thetaas(ii)*0.5).*mas(ii,:)'+(das(ii)*0.5*cos(thetaas(ii)*0.5)).*las(ii,:)'];
        qbs = [cos(thetabs(ii)*0.5);sin(thetabs(ii)*0.5).*lbs(ii,:)'];
        qdbs = [(-dbs(ii)*0.5)*sin(thetabs(ii)*0.5);sin(thetabs(ii)*0.5).*mbs(ii,:)'+(dbs(ii)*0.5*cos(thetabs(ii)*0.5)).*lbs(ii,:)'];
        %% formulation 1
        %             a = qas(2:4);b = qbs(2:4);
        %             ad = qdas(2:4);bd = qdbs(2:4);
        %             T((i-1)*6+1:(i-1)*6+3,:) = [a-b skewm(a+b) zeros(3,1) zeros(3,3)];
        %             T((i-1)*6+4:(i-1)*6+6,:) = [ad-bd skewm(ad+bd) a-b skewm(a+b)];
        %             As = As + T((i-1)*6+1:(i-1)*6+6,:)'*T((i-1)*6+1:(i-1)*6+6,:);
        %% formulation 2
        a = qas(1:4);b = qbs(1:4);
        ad = qdas(1:4);bd = qdbs(1:4);
        T(1:4,:) = [q2m_left(a)-q2m_right(b) zeros(4,4)];
        T(5:8,:) = [q2m_left(ad)-q2m_right(bd) q2m_left(a)-q2m_right(b)];
        As = As + T'*T;
    end
    As = (As+As')/2;
end

%%
% SDP£º
% cost min: xTAAx -> X = xTx Q = AA
function [T, time] = solve_SQP(x0, As)
    AA = As;
    n = size(x0,1);
    K = 3*n;
    A1 = blkdiag(eye(4),zeros(4));
    A2 = zeros(8);A2(1,5)=1;A2(2,6)=1;A2(3,7)=1;A2(4,8)=1;
    tic
    cvx_begin quiet
        variable X(n, n) symmetric
        minimize ( trace(AA*X) )
        subject to
            X == semidefinite(n);
            trace(A1*X) == 1;
            trace(A2*X) == 0;
    cvx_end
    
    lb = cvx_optval;
    X = (X + X')/2; % Force symmetry
    
    % Randomized algorithm
    mu = zeros(n,1); Sigma = X;
    % Using eigenvalue decomposition because chol() complains about
    % near-zero negative eigenvalues
    [V, D] = eig(Sigma);
    A = V*sqrt(max(D, 0));
%     A = chol(Sigma, 'lower');
    
    ub = 1e6; xhat = zeros(n, 1);
    for k = 1:K
        xf = (mulrandn_cached(mu, A));
        xf = make_feasible(xf, 1, A1);
        val = xf'*AA*xf;
%         [x, val] = ONEOPT(x, AA);
        if ub > val
            ub = val; xhat = xf;
        end
    end
    time = toc;
    x = xhat;
    disp(['scf :',num2str(time)]);
    q12 = x(1:4);
    q12d = x(5:8);
    R12 = q2r(q12);
    t12 = 2.*qprod(q12d, conjugateq(q12));
    t12 = t12(2:4);
    T = [R12 t12;[0 0 0 1]];
end

% Subroutine for multivariate normal distribution
%
% Samples from a multivariate normal distribution with mean mu
% and square root matrix of covariance, A (covariance = A*A').
%
% Reference:
% http://en.wikipedia.org/wiki/Multivariate_normal_distribution
function x = mulrandn_cached(mu, A)
    n = size(mu, 1);
    z = randn(n, 1);
    x = mu + A*z;
end

function xf = make_feasible(x, varargin)
    n = varargin{1};
    scales = zeros(n,1);
    for i = 1:n
        A = varargin{i+1};
        scales(i) = 1./sqrt(x'*A*x);
    end
    [~,id] = max(scales);
    xf = x .* scales(id);
end

% 1-opt greedy descent subroutine
%
% Starting from x, greedily chooses and optimizes over a single coordinate 
% to reduce the objective x^T P x + 2 q^T x the most. Stops when no single 
% coordinate change can reduce the objective. Returns the 1-opt point and 
% the function value.
function [x, val] = ONEOPT(x, P)
    g = 2*(P*x);
    v = diag(P);
    iters = 0;
    while true
        iters = iters + 1;
        if v >= abs(g)
            break;
        end
        c = (-g./(2*v));
        diffs = (c.^2).*v + c.*g;
        [~, i] = min(diffs);
        x(i) = x(i) + c(i);
        g = g + 2*c(i)*P(:, i);
    end
    val = x'*P*x;
end

% function [L,S] = conslin2(x)
%     f1 = @(x) (x(1)*x(4)+x(2)*x(5)+x(3)*x(6));
%     f2 = @(x) (x(1)*x(7)+x(2)*x(8)+x(3)*x(9));
%     f3 = @(x) (x(4)*x(7)+x(5)*x(8)+x(6)*x(9));
%     f4 = @(x) (x(1)*x(1)+x(2)*x(2)+x(3)*x(3));
%     f5 = @(x) (x(4)*x(4)+x(5)*x(5)+x(6)*x(6));
%     f6 = @(x) (x(7)*x(7)+x(8)*x(8)+x(9)*x(9));
%     
%     df1 = [x(4) x(5) x(6) x(1) x(2) x(3) 0 0 0 0 0 0 0];
%     df2 = [x(7) x(8) x(9) 0 0 0 x(1) x(2) x(3) 0 0 0 0];
%     df3 = [0 0 0 x(7) x(8) x(9) x(4) x(5) x(6) 0 0 0 0];
%     df4 = [2*x(1) 2*x(2) 2*x(3) 0 0 0 0 0 0 0 0 0 0];
%     df5 = [0 0 0 2*x(4) 2*x(5) 2*x(6) 0 0 0 0 0 0 0];
%     df6 = [0 0 0 0 0 0 2*x(7) 2*x(8) 2*x(9) 0 0 0 0];
%     
%     h = [0 0 0 1 1 1]';
%     L = [df1;df2;df3;df4;df5;df6];
%     b1 = [f1(x)-df1*x;f2(x)-df2*x;f3(x)-df3*x;f4(x)-df4*x;f5(x)-df5*x;f6(x)-df6*x];
%     S = h-b1;
%     L = [L;zeros(1,12) 1];
%     S = [S;1];
% end

%% dual quaternion based semidefinite programming
% X = x*x'=[  x1x1 x1x2 x1x3 x1x4 x1x5 x1x6 x1x7 x1x8;
%             x1x2 x2x2 x2x3 x2x4 x2x5 x2x6 x2x7 x2x8;
%             x1x3 x2x3 x3x3 x3x4 x3x5 x3x6 x3x7 x3x8;
%             x1x4 x2x4 x3x4 x4x4 x4x5 x4x6 x4x7 x4x8;
%             x1x5 x2x5 x3x5 x4x5 x5x5 x5x6 x5x7 x5x8;
%             x1x6 x2x6 x3x6 x4x6 x5x6 x6x6 x6x7 x6x8;
%             x1x7 x2x7 x3x7 x4x7 x5x7 x6x7 x7x7 x7x8;
%             x1x8 x2x8 x3x8 x4x8 x5x8 x6x8 x7x8 x8x8;];
% A1 = blkdiag(eye(4),zeros(4));
% A2 = zeros(8);A2(1,5)=1;A2(2,6)=1;A2(3,7)=1;A2(4,8)=1;
% function [L,S] = conslin5(x)
%     f1 = @(x) (x(1)*x(1)+x(2)*x(2)+x(3)*x(3)+x(4)*x(4));
%     f2 = @(x) (x(1)*x(5)+x(2)*x(6)+x(3)*x(7)+x(4)*x(8));
%     
%     df1 = [2*x(1) 2*x(2) 2*x(3) 2*x(4) 0 0 0 0];
%     df2 = [x(5) x(6) x(7) x(8) x(1) x(2) x(3) x(4)];
% 
%     h = [1 0]';
%     a1 = [df1;df2];
%     L = a1;
%     b1 = [f1(x)-df1*x;f2(x)-df2*x];
%     S = h - b1;
% end