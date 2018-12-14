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
    T1 = sol_park_martin(TA, TB, N);%sol_dual_quaternion
    if ~isempty(T1)
        T0 = T1;
    end
    dT = se3optimization(A, B, N, T0);
end

function dTopt = se3optimization(A, B, N, T0)
    maxIter = 200;
    % Solve for pose using our algorithm
    T = T0;
    xo = zeros(6,1);
    for i=1:maxIter      % Gauss-Newton iterations
        LHS = zeros(6);
        RHS = zeros(6,1);

        for k=1:N
            [G,e] = eJSE3(A{k},B{k},T);
            LHS = LHS + G'*G;
            RHS = RHS + G'*e;
        end
        xi = -LHS \ RHS;
        T = vec2tran( xi ) * T;
        if norm(xi-xo) < 1e-15
            break;
        end
        xo = xi;
    end
    dTopt = T;
    disp(['iter: ', num2str(i)]);
end

function [G,e] = eJSE3(A,B,T)
    se31 = tran2vec(B);
    AdT = tranAd(T);
    se32 = AdT*se31;
    invJ = vec2jacInv( se32 );
    That = vec2tran(se32);
    %% jacobian
    G = invJ * (eye(6) - tranAd(That));
    %% residual
    e1 = tran2vec(inv(A));
    e2 = tran2vec(B);
    e = e1 + AdT * e2;
end