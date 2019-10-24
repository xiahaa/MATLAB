function error=f_rot_error(X1,X2)
%     C1 = X1(1:3,1:3);
%     C1 = C1 * inv(sqrtm(C1'*C1));
%     C2 = X2(1:3,1:3);
%     C2 = C2 * inv(sqrtm(C2'*C2));
%     X1(1:3,1:3) = C1;
%     X2(1:3,1:3) = C2;

    M = (inv(X1)*X2);
    C1 = M(1:3,1:3);
    C1 = C1 * inv(sqrtm(C1'*C1));
    M(1:3,1:3) = C1;
%     error = norm(trace(M'*M));
    error = norm(tran2vec(M));
end