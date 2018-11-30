function q = dqconj(q)
    %% dual quaternion conjugate
    q(1:3) = -q(1:3);
    q(5:7) = -q(5:7);

end