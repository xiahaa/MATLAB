function MX = mean_Taylor_2nd_SO3( X, num )
    %% Calculate the 2nd order approximation of the mean of a bunch of matrices 
    % based on the Taylor expansion of the matrix logarithm and the definition 
    % of mean of a probability density function.
    %
    % Input: 
    %       X : a cell of 4 x 4 x n SE3 matrices
    % Output: 
    %       MX : the 2nd order approximation of the Taylor expansion
    %
    % Note : M_T1 is the mean of the 1st order approximation of Taylor
    % expansion
    %
    %
    % coder.extrinsic('mean_Taylor_1st_mex');
    %%
    M_t1 = mean_Taylor_1st( X );
    n = size(X,3);
    % Perturb M_t1 along each direction defined by the Lie algebra element of
    % SE(3)
    CA = inf;
    diff = 10^-12;
    eps  =  0.001;
    
    E = zeros(3,3,3);
    E(:,:,1) = [0  0 0; 0 0 -1;  0 1 0];
    E(:,:,2) = [0  0 1; 0 0  0; -1 0 0];
    E(:,:,3) = [0 -1 0; 1 0  0;  0 0 0];

    MX = M_t1;
    count = 0;
 
    step = 1; % Change the step lengths if descent direction is not found

    while( abs(CA) > diff )
        count = count + 1;
        n_step = 0;
        % Approximation of rotation is good so far so we only further optimize
        % on translation
        for j = 1:3
            MX1 = MX*expm( step*eps*E(:,:,j));
            MX2 = MX*expm(-step*eps*E(:,:,j));

            MX1sum = zeros(3);
            MX2sum = zeros(3);

            for k = 1:n
                X_k = X(:,:,k);
                MX1sum = MX1sum + X_k/MX1*X_k;
                MX2sum = MX2sum + X_k/MX2*X_k;
            end

            %% equation 39 in his paper
            MX1cost = 2*M_t1 - 1/2/n*MX1sum - 3/2*MX1;
            MX2cost = 2*M_t1 - 1/2/n*MX2sum - 3/2*MX2;

            CA1 = norm(MX1cost)^2;
            CA2 = norm(MX2cost)^2;

            if CA1 <= CA2 && CA1 < CA
                CA = CA1;
                MX = MX1;
                step = 1;
                % fprintf('Found a descent direction along %d \n', j)
                % disp(MX1cost)
                % disp(CA)
            elseif CA2 < CA1 && CA2 < CA
                CA = CA2;
                MX = MX2;
                step = 1;
                % fprintf('Found a descent direction along -%d \n', j)
                % disp(MX2cost)
                % disp(CA)
            else
                % fprintf('Not a descent direction along +/-%d \n', j)
                n_step = n_step + 1;
                if n_step == 3
                   step = step + 1;
                end
            end
        end
        if count == num
            break;
        end
    end
    % fprintf('\n')
    fprintf('Number of iterations is %f', count)
end