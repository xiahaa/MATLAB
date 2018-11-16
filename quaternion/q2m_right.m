function W = q2m_right(q)
%% quaternion to left matrix
%
%% Author: xiahaa@space.dtu.dk
    W = eye(4).*q(1) + blkdiag(0,[0 q(4) -q(3);-q(4) 0 q(2);q(3) -q(2) 0]);
    W(1,2:4) = -q(2:4);
    W(2:4,1) = q(2:4);
end