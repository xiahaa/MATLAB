function error=tran_error(X1, X2)
    %% return the translation error, inputs should be 2 4x4 transformation matrix
    error=norm(X1(1:3,4)-X2(1:3,4))/norm(X1(1:3,4));
end