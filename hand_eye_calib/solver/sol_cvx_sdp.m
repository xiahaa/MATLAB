function varargout = sol_cvx_sdp(TA,TB,N, varargin)
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
    [Ask,x0k] = preprocessing_kron(TA,TB,N);
        [T, time2] = solve_SDP(x0q, Asq);

%     if nargin == 4
%         va = varargin{1};
% %         vb = varargin{2};
%         [T, time4] = solve_SDP_kron(x0k, Ask, va);
%     else
%         [T, time4] = solve_SDP_kron(x0k, Ask);
%     end
    

    varargout{1} = T;
%     varargout{5} = time2;
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
% 
% f1 = @(x) (x(1)*x(4)+x(2)*x(5)+x(3)*x(6));
%     f2 = @(x) (x(1)*x(7)+x(2)*x(8)+x(3)*x(9));
%     f3 = @(x) (x(4)*x(7)+x(5)*x(8)+x(6)*x(9));
%     f4 = @(x) (x(1)*x(1)+x(2)*x(2)+x(3)*x(3));
%     f5 = @(x) (x(4)*x(4)+x(5)*x(5)+x(6)*x(6));
%     f6 = @(x) (x(7)*x(7)+x(8)*x(8)+x(9)*x(9));
function [T, time] = solve_SDP_kron(x0, As, varargin)
    AA = As;
    n = size(x0,1);
    K = 3*n;
    
    A1 = zeros(n,n);
    A1(1,4) = 1;A1(2,5) = 1;A1(3,6) = 1;
    A2 = zeros(n,n);
    A2(1,7) = 1;A2(2,8) = 1;A2(3,9) = 1;
    A3 = zeros(n,n);
    A2(4,7) = 1;A2(5,8) = 1;A2(6,9) = 1;
    A4 = zeros(n,n);
    A4(1,1) = 1;A4(2,2) = 1;A4(3,3) = 1;
    A5 = zeros(n,n);
    A5(4,4) = 1;A5(5,5) = 1;A5(6,6) = 1;
    A6 = zeros(n,n);
    A6(7,7) = 1;A6(8,8) = 1;A6(9,9) = 1;
    A7 = zeros(n,n);
    A7(n,n) = 1;
    
    ext = 0;
    if nargin == 3
        va = varargin{1};
        Aext = zeros(n,n); Aext(1,13) = va(1);Aext(2,13) = va(2);Aext(3,13) = va(3);
    end
    
    
    tic
    cvx_begin quiet
        variable X(n, n) symmetric
        minimize ( trace(AA*X) )
        subject to
            X == semidefinite(n);
            trace(A1*X) == 0;
            trace(A2*X) == 0;
            trace(A3*X) == 0;
            trace(A4*X) == 1;
            trace(A5*X) == 1;
            trace(A6*X) == 1;
            trace(A7*X) == 1;
            if ext == 1
                trace(Aext*X) > 0;
            end
    cvx_end
    lb = cvx_optval;
    X = (X + X')/2; % Force symmetry
    % Randomized algorithm
    mu = zeros(n,1); Sigma = X;
    % Using eigenvalue decomposition because chol() complains about
    % near-zero negative eigenvalues
    [V, D] = eig(Sigma);
    A = V*sqrt(max(D, 0));    
    ub = 1e6; xhat = zeros(n, 1);
    for k = 1:K
        xf = (mulrandn_cached(mu, A));
        xf = make_feasible(xf, 4, A4, A5, A6, A7);
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
    R12 = [x(1) x(4) x(7); ...
           x(2) x(5) x(8); ...
           x(3) x(6) x(9)];
    t12 = x(10:12);
%     R12 = R12*inv(sqrtm(R12'*R12));
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
    K = 5*n;
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