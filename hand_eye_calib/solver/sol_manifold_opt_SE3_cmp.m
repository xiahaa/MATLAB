function varargout = sol_manifold_opt_SE3_cmp(TA, TB, N)
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
    convSolver = {@sol_andreff, ...                               %% LIE
%                   @, ...                                   %% KR
%                   @HandEye_DQ, ...                                    %% DQ
                    };
    num = size(convSolver,2)+1;
    dT = cell(num);
    
    dT{1} = se3optimization(A, B, N, T0);
    varargout{1} = dT{1};
    for i = 1:size(convSolver,2)
        handle_sol = convSolver{i};
        T0 = handle_sol(TA, TB, N);
        dT{i+1} = se3optimization(A, B, N, T0);
        varargout{i+1} = dT{i+1};
    end
end

function dTopt = se3optimization(A, B, N, T0)
    maxIter = 1000;    
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
        
        if norm(xi-xo) < 1e-12
            break;
        end
        xo = xi;
    end
    dTopt = T;
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