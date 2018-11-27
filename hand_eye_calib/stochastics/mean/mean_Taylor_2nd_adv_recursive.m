function MX = mean_Taylor_2nd_adv_recursive( X, noise, num )
%% This function calculates the 2nd order approximation of the mean of a
% bunch of matrices based on the Taylor expansion of the matrix logarithm
% and the definition of mean of a probability density function.

% Input: X is a cell of 4 by 4*n rigid transformation matrices
% Output: M_T1 is the mean of the 1st order approximation of Taylor
% expansion

% Output: MX is the 2nd order approximation of the Taylor expansion
% coder.extrinsic('mean_Taylor_1st_mex');

    %%
    n = size(X,2)/4;
    M_t1 = zeros(4);
    M_t1 = mean_Taylor_1st( X );

    % Perturb M_t1 along each direction defined by the Lie algebra element of
    % SE(3)
    costfun_min = 10^-15;%10^-12;

    E = zeros(4,4,6);
    E(:,:,1) = [0  0 0 0; 0 0 -1 0;  0 1 0 0; 0 0 0 0];
    E(:,:,2) = [0  0 1 0; 0 0  0 0; -1 0 0 0; 0 0 0 0];
    E(:,:,3) = [0 -1 0 0; 1 0  0 0;  0 0 0 0; 0 0 0 0];
    E(:,:,4) = [0  0 0 1; 0 0  0 0;  0 0 0 0; 0 0 0 0];
    E(:,:,5) = [0  0 0 0; 0 0  0 1;  0 0 0 0; 0 0 0 0];
    E(:,:,6) = [0  0 0 0; 0 0  0 0;  0 0 0 1; 0 0 0 0];

    MX = M_t1;

    count = 0;


    % fprintf('---------------------------------------------- \n')
    % fprintf('Search gradient descent direction from %f to 6 \n', m)
    % fprintf('---------------------------------------------- \n')

    while(1)
        count = count + 1;

        MXsum = zeros(4);
        for k = 1:n
            X_k = X(:,(k-1)*4+1:k*4);
            MXsum = MXsum + X_k/MX*X_k;
        end

        MXcost = 2*M_t1 - 1/2/n*MXsum - 3/2*MX;
        
        if noise == 0
            CA = norm(MXcost(1:3,4))^2;
        elseif noise == 1
            CA = norm(MXcost)^2;
        else
            CA = -1;
        end
        
        fprintf('Number of iterations is %f \n', count);
        fprintf('Final Error CA for this loop is %e \n', CA );
        
        if abs(CA) < costfun_min
            break;
        end
        %% equ 43b
        b = - [ MXcost(:,1);MXcost(:,2);MXcost(:,3);MXcost(:,4)];

        J_1 = zeros(16,16);
        %% equ 43a
        for k = 1:n
            X_k = X(:,(k-1)*4+1:k*4);
            C = MX\X_k;
            J_k = kron(transpose(C),X_k);
            J_1 = J_1 + J_k;
        end
        J = 1/2/n * J_1 -3/2 * kron(eye(4),MX);

        X_vector = J\b;
        X_gradient = [X_vector(1:4) X_vector(5:8) X_vector(9:12) X_vector(13:16)];
        MX_translation = MX * (eye(4) + X_gradient);

        if noise == 0
            MX(1:3,4) = MX_translation(1:3,4);
        elseif noise == 1
            MX = MX_translation;
        end
    end
    MX = orthog(MX);
    fprintf('loop \n');
end


