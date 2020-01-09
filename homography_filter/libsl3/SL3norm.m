function n = SL3norm(H1,H2)
% Frobenius norm of SL3
    n = sqrt(trace(H1'*H2));
end