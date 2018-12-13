function varargout = sol_cvx_sdp(TA,TB,N)
%% solving hand eye calibration using SDP
%% Author: xiahaa@space.dtu.dk
    addpath('./solver/sdp');
    if N < 2
        error('At least two samples needed for unique solution!');
        varargout{1} = [];
        return;
    end
    
    format long;
        
    Asq = preprocessing_q(TA,TB,N);
    x0q = [1 0 0 0 1 0 0 0]';
    
    [T, time2] = solve_SDP(x0q, Asq);
    

    varargout{1} = T;
%     varargout{5} = time2;
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

        if abs(thetaa-thetab) < 0.15 && isempty(find(isnan(la),1)) && isempty(find(isnan(ma),1))
            Nv = Nv + 1;
            ids = [ids;i];
        end
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
% SDP??
% cost min: xTAAx -> X = xTx Q = AA
function [T, time] = solve_SDP(x0, As)
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
%         [xf, val] = ONEOPT(xf, AA);
        if ub > val
            ub = val; xhat = xf;
        end
    end
    xstar = xhat;
    
    time = toc;
    x = xstar;
    disp(['SDP + randomization :',num2str(time)]);
    q12 = x(1:4);
    q12d = x(5:8);
    R12 = q2r(q12);
    t12 = 2.*qprod(q12d, conjugateq(q12));
    t12 = t12(2:4);
    T = [R12 t12;[0 0 0 1]];
end





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