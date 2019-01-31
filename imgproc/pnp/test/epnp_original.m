function [R,t] = epnp_original(p,q,K)
    n = size(p,2);
    if size(q,1) == 2
        q = [q;ones(1,n)];
    end
    qn = K\q;
    
    %% random sampling
    addpath ../../../../3rdparty/code3/epnp/
    
    [R,t]= efficient_pnp(p',q',K);
end