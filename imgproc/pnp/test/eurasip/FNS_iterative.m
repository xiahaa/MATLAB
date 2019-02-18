function M = FNS_iterative(SO3s,M_1st)
% My implementation of iterative optimizing the mean using manfiold
% optimization technique. It works directly with the definition of the mean
% and applies optimization directly on SO3.
    M = M_1st;
    iter = 1;
    maxIter = 50;
    while iter < maxIter
        Minv = (M)';
        LHS = zeros(3,3);
        RHS = zeros(3,1);
        for k=1:size(SO3s,3)
            v1 = rot2vec(Minv*SO3s(:,:,k));
            RHS = RHS + v1;
            LHS = LHS + vec2jacInv(v1);
        end
        xhat = LHS\RHS;
        M = M*vec2rot(xhat);
        if norm(xhat) < 1e-3
            break;
        end
        iter = iter + 1;
    end
end