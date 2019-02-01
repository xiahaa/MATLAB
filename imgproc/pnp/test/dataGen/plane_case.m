function [p, q, R, t] = plane_case(npt)
    % generate 3d coordinates (all on a plane) in camera space
    p = [xrand(2,npt,[-2 2]); zeros(1,npt)];
    R = rodrigues(randn(3,1));
    t = [rand-0.5;rand-0.5;rand*8+4];
    q = R*p+repmat(t,1,npt);
end

