function [R,t] = pnp_sdr(P, q, K, R, t)
    
    addpath('../../../hand_eye_calib/solver/sdp/');

    n = size(q,2);
    if size(q,1) == 2
        q = [q;ones(1,n)];
    end
    qn = K\q;
%     vi = qn;
    vi = qn ./ sqrt(qn(1,:).^2+qn(2,:).^2+qn(3,:).^2);
    
    %% formulate
    nopt = 9 + 1;% rotation + translation + slack
    Qi = zeros(3,3,n);
    QiCci = zeros(3,nopt,n);
    F1 = zeros(3,3);
    F2 = zeros(3,nopt);
    for i = 1:n
        Vi = vi(:,i)*vi(:,i)'/(vi(:,i)'*vi(:,i));
        Qi(:,:,i) = (eye(3)-Vi)'*(eye(3)-Vi);
        F1 = F1 + Qi(:,:,i);
        
        Cd = [P(:,i)' 0 0 0 0 0 0 -qn(1,i); ...
              0 0 0 P(:,i)' 0 0 0 -qn(2,i); ...
              0 0 0 0 0 0 P(:,i)' -qn(3,i)];
        QiCci(:,:,i) = Qi(:,:,i) * Cd;
        F2 = F2 + QiCci(:,:,i);
    end
    F1inv = inv(F1);
    F3 = -F1inv * F2;
    M = zeros(nopt, nopt);
    for i = 1:n
        Cd = [P(:,i)' 0 0 0 0 0 0 -qn(1,i); ...
              0 0 0 P(:,i)' 0 0 0 -qn(2,i); ...
              0 0 0 0 0 0 P(:,i)' -qn(3,i)];
        M = M + (Cd + F3)'*Qi(:,:,i)*(Cd + F3);
    end
    M = (M+M').*0.5;% force symmetric
    xxx = [R(1,:)';R(2,:)';R(3,:)';1];
    
    ext = [0 0 0 0 0 0 P(:,i)' 0] + F3(3,:);
    [R, x, time] = solve_SDP_kron(nopt, M, ext);
    t = F3*x;
end

function [R, x, time] = solve_SDP_kron(nopt, Q, varargin)
    AA = Q;
    n = nopt;
    
    K = 3*n;
    
    A1 = zeros(n,n);
    A1(1,2) = 1;A1(4,5) = 1;A1(7,8) = 1;
    A2 = zeros(n,n);
    A2(1,3) = 1;A2(4,6) = 1;A2(7,9) = 1;
    A3 = zeros(n,n);
    A2(2,3) = 1;A2(5,6) = 1;A2(8,9) = 1;
    A4 = zeros(n,n);
    A4(1,1) = 1;A4(4,4) = 1;A4(7,7) = 1;
    A5 = zeros(n,n);
    A5(2,2) = 1;A5(5,5) = 1;A5(8,8) = 1;
    A6 = zeros(n,n);
    A6(3,3) = 1;A6(6,6) = 1;A6(9,9) = 1;
%     A7 = zeros(n,n);
%     A7(n,n) = 1;
%     format long;
    ext = 0;
    if nargin == 3
        va = varargin{1};
        Aext = zeros(n,n); 
        %% r31*P(1)+r32*P(2)+r33*P(3)+t(3) > 0 chirality constraint
%         Aext(7,end) = va(1);Aext(8,end) = va(2);Aext(9,end) = va(3);Aext(12,13) = -1;
        Aext(:,end) = va';
        ext = 1;
    end
    
    tic
    cvx_begin quiet
        variable X(n, n) symmetric
        minimize ( trace(AA*X) )
        subject to
            X == semidefinite(n);
%             X >= 0
            trace(A1*X) == 0;
            trace(A2*X) == 0;
            trace(A3*X) == 0;
            trace(A4*X) == 1;
            trace(A5*X) == 1;
            trace(A6*X) == 1;
%             trace(A7*X) == 1;
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
        xf = make_feasible(xf, 3, A4, A5, A6);
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
    R12 = [x(1) x(2) x(3); ...
           x(4) x(5) x(6); ...
           x(7) x(8) x(9)];
%     t12 = x(10:12);
%     R12 = R12*inv(sqrtm(R12'*R12));
    R = R12;
%     t = t12;
end