function [p, q, R, t] = quasi_sigular_case(npt)
    % generate 3d coordinates in camera space
    q = [xrand(1,npt,[1 2]); xrand(1,npt,[1 2]); xrand(1,npt,[4 8])];
    t = mean(Xc,2);
    R = rodrigues(randn(3,1));
    p = inv(R)*(q-repmat(t,1,npt));
end

