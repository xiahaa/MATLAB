function [theta, d, l, m] = screwParams(R,t)
    l = rot2vec(R);
    theta = norm(l);
    l = l./(theta+eps);
    d = t'*l;
    txl = skewm(t)*l;
    m = 0.5.*(txl+skewm(l)*txl.*cot(theta*0.5));
end