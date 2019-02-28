
function [M_1, M_hat] = mean_Taylor_1st1( X ) 
    % This function calculates the 1st order approximation of the mean of a
    % bunch of matrices based on the Taylor expansion of the matrix logarithm
    % and the definition of mean of a probability density function.

    % Input: X is a 4 by 4*n rigid transformation matrices
    % Output: M_T1 is the mean of the 1st order approximation of Taylor
    % expansion

    % Change of this m file doesn't automatically change the executable generated by 
    % mean_Taylor_2nd.m
    n =  size(X, 3);

    M_hat = zeros(3);
    M_1 = zeros(3);

    for i = 1:n
        M_hat = M_hat + X(1:3,1:3,i);
    end

    M_hat = 1/n*M_hat;  % Note that M_hat doesn't belong to SE(3)

    M_1 = orthog(M_hat);
end
