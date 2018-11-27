function A = ivecm(a)
    %% reverse vectorize a matrix, only works for square matrix
    n = numel(a);
    ns = sqrt(n);
    A = reshape(a,ns,ns);
end