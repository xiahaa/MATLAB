function [X] = sol_andreff(TA,TB,N)
% Solves the problem AX=XB
% using the formulation of
%
% On-line Hand-Eye Calibration.
% N. Andreff, R. Horaud, B. Espiau 
%
% Mili Shah
% July 2014
    dim = size(TA,2);
    n = N;
    A = zeros(12*n,12);
    b = zeros(12*n,1);
    for i = 1:n
        T1 = TA(i,:,:);
        T2 = TB(i,:,:);
        AA = reshape(T2,dim,dim,1);BB = reshape(T1,dim,dim,1);
        Ra = AA(1:3,1:3);
        Rb = BB(1:3,1:3);
        ta = AA(1:3,4);
        tb = BB(1:3,4);
        
        A(12*i-11:12*i-3,1:9) = eye(9) - kron(Rb,Ra);
        A(12*i-2:12*i,:) = [kron(tb',eye(3)) eye(3)-Ra];
        b(12*i-2:12*i) = ta;
    end
    x = A\b;

    X = reshape(x(1:9),3,3)';
    X = sign(det(X))/abs(det(X))^(1/3)*X;

    [u, ~, v] = svd(X); X = u*v'; if det(X)<0, X = u*diag([1 1 -1])*v'; end
    X = [X' x(10:12);[0 0 0 1]];
end
