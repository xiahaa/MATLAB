function Sigma = cov_SE3(MX, X, order)
    %% compute the covariance of a batch of matrices
    
    n = size(X,2) / 4;
    Sigma = zeros(6,6);
    
    if order == 2
        %% do second order approximation
        for i = 1:n
            X_i = X(:,4*(i-1)+1:4*i);
            P = (MX\X_i - eye(4)) - 0.5*((MX\X_i - eye(4)))^2;
            Sigma = Sigma + se3_vec(P)*se3_vec(P)';
        end
    else
        %% do first order approximation
        for i = 1:n
            X_i = X(:,4*(i-1)+1:4*i);
%             P = (MX\X_i - eye(4));
%             Sigma = Sigma + se3_vec(P)*se3_vec(P)';
              P = MX\X_i;
              Sigma = Sigma + se3_vec(logm(P))*se3_vec(logm(P))';
        end
    end
    Sigma = Sigma ./ n;
end