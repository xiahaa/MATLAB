function vecA = vecm(A)
    %% vectorize a matrix
    vecA = reshape(A,size(A,1)*size(A,2),1);
end