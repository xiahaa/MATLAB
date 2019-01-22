function [R,t] = pnp_sdr(P, q, K, R, t)
    
    addpath('../../../hand_eye_calib/solver/sdp/');

    n = size(q,2);
    if size(q,1) == 2
        q = [q;ones(1,n)];
    end
    qn = K\q;
    vi = qn;
%     vi = qn ./ sqrt(qn(1,:).^2+qn(2,:).^2+qn(3,:).^2);
    
    %% formulate
    nopt = 9 + 3;% rotation + translation + slack
    Q = zeros(nopt,nopt);
    for i = 1:n
        Vi = vi(:,i)*vi(:,i)'/(vi(:,i)'*vi(:,i));
        IVVI = (eye(3)-Vi)'*(eye(3)-Vi);
        A = [P(1,i) 0 0 P(2,i) 0 0 P(3,i) 0 0 1 0 0; ...
             0 P(1,i) 0 0 P(2,i) 0 0 P(3,i) 0 0 1 0; ...
             0 0 P(1,i) 0 0 P(2,i) 0 0 P(3,i) 0 0 1];
        Q = Q + A'*IVVI*A;
    end
    Q = (Q+Q').*0.5;% force symmetric
%     xxx = [R(1,:)';R(2,:)';R(3,:)';t];
    
    [R, t, time] = solve_SDP_kron(nopt, Q, P(:,1));
end

function [R, t, time] = solve_SDP_kron(nopt, Q, varargin)
    AA = Q;
    n = nopt;
    
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
%     format long;
    ext = 0;
    if nargin == 3
        va = varargin{1};
        Aext = zeros(n,n); 
        %% r31*P(1)+r32*P(2)+r33*P(3)+t(3) > 0 chirality constraint
        Aext(7,13) = va(1);Aext(8,13) = va(2);Aext(9,13) = va(3);Aext(12,13) = -1;
    end
    
    tic
    cvx_begin 
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
%             X(n,n) == 1;
            if ext == 1
%                 trace(Aext*X) > 0;
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
    R12 = [x(1) x(4) x(7); ...
           x(2) x(5) x(8); ...
           x(3) x(6) x(9)];
    t12 = x(10:12);
%     R12 = R12*inv(sqrtm(R12'*R12));
    R = R12;
    t = t12;
end