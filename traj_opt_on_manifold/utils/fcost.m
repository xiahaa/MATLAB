function y = cost(xi,v,tau,lambda,miu,varargin)
    if nargin > 5
        alpha = varargin{1};
    else
        alpha = 1;
    end
    % cost term 1, data cost
    cost1 = sum(vecnorm(xi,2).^2.*2);

    % cost term 2, first order smooth cost, integrate with trapezoidal
    % rule, consistent with Boumal's paper. TODO change in paper.
    N = size(v,2)+1;
    wv = [1 ones(1,N-2)];
    cost2 = sum(vecnorm(v,2).^2.*(2/tau).*wv);

    % cost term 3, second order smooth cost, integrate with trapezoidal
    % rule
    a = zeros(3,N-2);
    for i = 2:N-1
        a(:,i-1)=v(:,i)-v(:,i-1);
    end
    cost3 = sum(vecnorm(a,2).^2.*(2/tau^3));

    y = alpha.*cost1 * 0.5 + cost2 * 0.5 * lambda + cost3 * 0.5 * miu;
end
