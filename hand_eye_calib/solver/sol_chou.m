function [X] = sol_chou(TA,TB,N)
% Calculates the least squares solution of
% AX = XB
% From
% Finding the Position and Orientation of a
% Sensor on a Robot Manipulator Using Quaternions.
% Chou and Kamel
%
% Mili Shah
% July 2014
%
% Uses q2rot.m and rot2q.m
    n = N;
    AA = zeros(4*n,4);
    dim = size(TA,2);
    %Calculate best rotation R
    for i = 1:n
        T1 = TA(i,:,:);
        T2 = TB(i,:,:);
        A = reshape(T2,dim,dim,1);B = reshape(T1,dim,dim,1);
        
        A1 = (A(1:3,1:3));
        B1 = (B(1:3,1:3));
        a = rot2quat(A1);
        b = rot2quat(B1);
        AA(4*i-3:4*i,:) = ...
            [a(1) -a(2:4)';a(2:4) a(1)*eye(3)+skewm(a(2:4))]...
            -...
            [b(1) -b(2:4)';b(2:4) b(1)*eye(3)-skewm(b(2:4))];
    end
    [~,~,v]=svd(AA); v = v(:,4); v = [v(1);v(2:4)];
    R = q2r(v);
    %Calculate best translation t
    C = zeros(3*n,3);
    d = zeros(3*n,1);
    I = eye(3);
    for i = 1:n
        T1 = TA(i,:,:);
        T2 = TB(i,:,:);
        A = reshape(T2,dim,dim,1);B = reshape(T1,dim,dim,1);
        C(3*i-2:3*i,:) = I - A(1:3,1:3);
        d(3*i-2:3*i,:) = A(1:3,4)-R*B(1:3,4);
    end
    t = C\d;
    %Put everything together to form X
    X = [R t;0 0 0 1];
end

