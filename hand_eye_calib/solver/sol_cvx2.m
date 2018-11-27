function varargout = sol_cvx2(TA,TB,N)
%% solving hand eye calibration using SCFS

%% Author: xiahaa@space.dtu.dk
    dim = size(TA,2);
    if N < 2
        error('At least two samples needed for unique solution!');
        varargout{1} = [];
        return;
    end
    if dim < 4
        error('Only work for T!');
        varargout{1} = [];
        return;
    end
    
    format long;
    
    C = zeros(N*12,12);
    d = zeros(N*12,1);
    for i = 1:N
        T1 = TA(i,:,:);T1 = reshape(T1,dim,dim,1);
        T2 = TB(i,:,:);T2 = reshape(T2,dim,dim,1);
        Ra = T2(1:3,1:3);ta = T2(1:3,4);
        Rb = T1(1:3,1:3);tb = T1(1:3,4);
        
        id = (i-1)*12;
        C(id+1:id+12,:) = [kron(eye(3),Ra)-kron(Rb',eye(3)) zeros(9,3); ...
                           kron(tb', eye(3)) eye(3)-Ra];
        d(id+1:id+12) = [zeros(9,1);ta];
    end
    
    As = C'*C;
    bs = d'*C;
        
    %% vector to optimize
%     x0 = C\d;
    x0 = [1 0 0 0 1 0 0 0 1 0 0 0]';
    x = x0;
    for k = 1:1000
        [L,S] = conslin(x);
        disp(S);
        soln = quadprog(2*As,(-2*bs)', [],[],L, S);
        xnew = soln;
        if norm(xnew-x) < 1e-12
            disp(['converge at step ', num2str(k)]);
            break;
        end
        x = xnew;
    end
    
    R12 = [x(1) x(4) x(7); ...
           x(2) x(5) x(8); ...
           x(3) x(6) x(9)];
    t12 = x(10:12);
    
    R12 = R12*inv(sqrtm(R12'*R12));
    
    if dim == 4
        varargout{1} = [R12 t12;[0 0 0 1]];
    else
        varargout{1} = R12;
    end
end

function [L,S] = conslin(x)
    f1 = @(x) (x(1)*x(4)+x(2)*x(5)+x(3)*x(6));
    f2 = @(x) (x(1)*x(7)+x(2)*x(8)+x(3)*x(9));
    f3 = @(x) (x(4)*x(7)+x(5)*x(8)+x(6)*x(9));
    f4 = @(x) (x(1)*x(1)+x(2)*x(2)+x(3)*x(3));
    f5 = @(x) (x(4)*x(4)+x(5)*x(5)+x(6)*x(6));
    f6 = @(x) (x(7)*x(7)+x(8)*x(8)+x(9)*x(9));
    
    df1 = [x(4) x(5) x(6) x(1) x(2) x(3) 0 0 0 0 0 0];
    df2 = [x(7) x(8) x(9) 0 0 0 x(1) x(2) x(3) 0 0 0];
    df3 = [0 0 0 x(7) x(8) x(9) x(4) x(5) x(6) 0 0 0];
    df4 = [2*x(1) 2*x(2) 2*x(3) 0 0 0 0 0 0 0 0 0];
    df5 = [0 0 0 2*x(4) 2*x(5) 2*x(6) 0 0 0 0 0 0];
    df6 = [0 0 0 0 0 0 2*x(7) 2*x(8) 2*x(9) 0 0 0];
    
    b = [0 0 0 1 1 1]';
    
    L1 = [df1;df2;df3;df4;df5;df6];
    fstar = [f1(x);f2(x);f3(x);f4(x);f5(x);f6(x)];
    fhatstar = L1 * x;
    S1 = b + fhatstar - fstar;
    L = [L1;];%-1.*L1];
    S = [S1;];%-1.*S1];
end