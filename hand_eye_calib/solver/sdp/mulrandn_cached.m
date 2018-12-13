% Subroutine for multivariate normal distribution
%
% Samples from a multivariate normal distribution with mean mu
% and square root matrix of covariance, A (covariance = A*A').
%
% Reference:
% http://en.wikipedia.org/wiki/Multivariate_normal_distribution
function x = mulrandn_cached(mu, A)
    n = size(mu, 1);
    z = randn(n, 1);
    x = mu + A*z;
end

