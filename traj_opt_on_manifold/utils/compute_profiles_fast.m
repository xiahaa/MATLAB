function [speed, acc] = compute_profiles_fast(dtau,X)
    Nd = size(X,3);
    speed = zeros(1, Nd);
    for k = 1 : Nd-1
        fw = logSO3(X(:, :, k)'*X(:, :, k+1));
        v = fw/dtau;
        speed(k) = sqrt(2)*norm(v);
    end
    % Backward difference for last point.
    speed(end) = speed(end-1);
    % Acceleration is NaN at first and last point. For all the others,
    % using a symmetric difference formula.
    acc = NaN(1, Nd);
    for k = 2 : Nd-1
        fw = logSO3(X(:, :, k)'*X(:, :, k+1));
        bw = logSO3(X(:, :, k)'*X(:, :, k-1));
        a = ( fw + bw ) / ( dtau^2 );
        acc(k) = sqrt(2)*norm(a);
    end

end
