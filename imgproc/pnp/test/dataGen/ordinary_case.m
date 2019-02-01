function [p, q, R, t] = ordinary_case(npt)
    % generate 3d coordinates in camera space
    q = [xrand(1,npt,[-2 2]); xrand(1,npt,[-2 2]); xrand(1,npt,[4 8])];
    t = mean(q,2);
    R= rodrigues(randn(3,1));
    p = inv(R)*(q-repmat(t,1,npt));
end


