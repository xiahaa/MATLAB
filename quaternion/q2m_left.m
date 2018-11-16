function Q = q2m_left(q)
%% quaternion to left matrix
%
%% Author: xiahaa@space.dtu.dk
    Q = eye(4).*q(1) + blkdiag(0,[0 -q(4) q(3);q(4) 0 -q(2);-q(3) q(2) 0]);
    Q(1,2:4) = -q(2:4);
    Q(2:4,1) = q(2:4);
end