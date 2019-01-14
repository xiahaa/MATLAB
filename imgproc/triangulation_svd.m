function [ Q ] = triangulation_svd( q1, P1, q2, P2 )
%% triangulation 1 point
% Author: xiahaa@sapce.dtu.dk
    A = zeros(4,3);
    b = zeros(4,1);
    Q = zeros(3,size(q1,2));
    for i = 1:size(q1,2)
        u1 = q1(1,i);
        v1 = q1(2,i);

        u2 = q2(1,i);
        v2 = q2(2,i);

        A(1,:) = [P1(1,1)-u1*P1(3,1) P1(1,2)-u1*P1(3,2) P1(1,3)-u1*P1(3,3)];
        A(2,:) = [P1(2,1)-v1*P1(3,1) P1(2,2)-v1*P1(3,2) P1(2,3)-v1*P1(3,3)];
        A(3,:) = [P2(1,1)-u2*P2(3,1) P2(1,2)-u2*P2(3,2) P2(1,3)-u2*P2(3,3)];
        A(4,:) = [P2(2,1)-v2*P2(3,1) P2(2,2)-v2*P2(3,2) P2(2,3)-v2*P2(3,3)];

        b(1,1) = [u1*P1(3,4)-P1(1,4)];
        b(2,1) = [v1*P1(3,4)-P1(2,4)];
        b(3,1) = [u2*P2(3,4)-P2(1,4)];
        b(4,1) = [v2*P2(3,4)-P2(2,4)];

        Q(:,i) = A\b;
    end
end

