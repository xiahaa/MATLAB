function [R,t] = pnp_gOp(P, q, K)
% Interface to 
% Globally Optimal O(n) Solution to the PnP Problem for General Camera Models
% Gerald Schweighofer and Axel Pinz, BMVC.

    addpath C:\Users\xiahaa\Documents\MATLAB\imgproc\pnp\3rdparty/GlobalOptimalPose/gOp
    addpath C:\Users\xiahaa\Documents\MATLAB\imgproc\pnp\3rdparty/GlobalOptimalPose/gOp/util

    opt.acc = 1e-8; 
    opt.methode ='3D';%'3D';planar

    % generate a test model 
    n = size(P,2);
    c = rand(3,n)*5;
    c=c*0;     %% for simple perspective camera 

    if size(q,1) == 2
        q = [q;ones(1,n)];
    end
    qn = K\q;
    
    %% generate a measurement 
    v = qn ./ sqrt(qn(1,:).^2+qn(2,:).^2+qn(3,:).^2);

    %% call gOp
    [R,t] = gOp([v;c],P,opt);
    
end