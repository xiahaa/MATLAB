function [R,t] = pnp_robust(p,q,K)
%% implementation of PnP algorithm proposed in
%% A Robust O(n) Solution to the Perspective-n-Point Problem
%% Shiqi Li, Chi Xu and Ming Xie, IEEE PAMI, 2012.
    n = size(p,2);
    if size(q,1) == 2
        q = [q;ones(1,n)];
    end
    qn = K\q;
    
    %% random sampling
    addpath C:/Users/xiahaa/Documents/MATLAB/imgproc/pnp/3rdparty/OPnP_Toolbox_Original/rpnp1.0/code3/func
    
    [R,t]= RPnP(p,qn(1:2,1:n));
    
end