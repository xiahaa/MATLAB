function dT = sol_manifold_opt_SE3(TA, TB, N)
%% input should be Nx4x4 3D matrices.
    if iscell(TA) == false
        %% 1st, prepare data, transform to cell array
        for j = 1:N
            Ts = TA(j,:,:);Ts = reshape(Ts,4,4);
            B{j} = Ts;
            Ts = TB(j,:,:);Ts = reshape(Ts,4,4);
            A{j} = Ts;
        end
    else
        A = TA;
        B = TB;
    end
    T0 = [eye(3) [0;0;0];0 0 0 1];
%     T1 = sol_horaud(TA, TB, N);%sol_dual_quaternion
%     if ~isempty(T1)
%         T0 = T1;
%     end
    dT = se3optimization(A, B, N, T0);
end

function dTopt = se3optimization(A, B, N, T0)
    maxIter = 200;
    % Solve for pose using our algorithm
    T = T0;
    xo = zeros(6,1);
    
    %% some precompute values
    se31 = zeros(6,N);
    e1 = zeros(6,N);
    curlya = zeros(6, 6, N);
    for k=1:N
        se31(:,k) = tran2vec(B{k});
        e1(:,k) = tran2vec((A{k}));
        curlya(:,:,k) = curlyhat(e1(:,k)); 
    end
    
    for i=1:maxIter      % Gauss-Newton iterations
        LHS = zeros(6);
        RHS = zeros(6,1);
        
        AdT = tranAd(T);
        se32 = AdT * se31;
        es = -e1 + se32;
        for k=1:N
            [G] = eJSE3fast(se32(:,k));
            G2 = (eye(6)- 0.5.*curlya(:,:,k))*G;
            e2 = es(:,k) - 0.5.*(curlya(:,:,k)*se32(:,k));
            LHS = LHS + G2'*G2;
            RHS = RHS + G2'*e2;
        end
%         tic
        xi = -LHS \ RHS;
%         toc
%         tic
%         R = chol(LHS);
%         xi = -R\(R'\RHS);
%         t1 = toc;
%         disp(['chol ',num2str(t1)]);
        T = vec2tran( xi ) * T;
        if norm(xi-xo) < 1e-12
            break;
        end
        xo = xi;
    end
    dTopt = T;
    disp(['iter: ', num2str(i)]);
end

function [G] = eJSE3fast(se32)
%     se32 = AdT*se31;
    invJ = vec2jacInv( se32 );
    That = vec2tran(se32);
    %% jacobian
%     G = invJ * (eye(6) - tranAd(That));
    G = invJ - invJ * tranAd(That);
    %% residual
%     e = e1 + AdT * se31;
end