function M = mean_1st_order(SO3s)
% first order mean.
    addpath C:\Users\xiahaa\Documents\MATLAB\MatrixLieGroup\barfoot_tro14\
    M = zeros(3,3);
    for i = 1:size(SO3s,3)
        M = M + SO3s(:,:,i);
    end
    M = M ./ size(SO3s,3);
    %% orthogonalization
%     [U, ~, V] = svd(M); 
%     M = U*V';
%     if det(M)<0, M = U*diag([1 1 -1])*V'; end 
     M = orthog(M);
end