% 1-opt greedy descent subroutine
%
% Starting from x, greedily chooses and optimizes over a single coordinate 
% to reduce the objective x^T P x + 2 q^T x the most. Stops when no single 
% coordinate change can reduce the objective. Returns the 1-opt point and 
% the function value.
function [x, val] = ONEOPT(x, P)
    g = 2*(P*x);
    v = diag(P);
    iters = 0;
    while true
        iters = iters + 1;
        if v >= abs(g)
            break;
        end
        c = (-g./(2*v));
        diffs = (c.^2).*v + c.*g;
        [~, i] = min(diffs);
        x(i) = x(i) + c(i);
        g = g + 2*c(i)*P(:, i);
    end
    val = x'*P*x;
end