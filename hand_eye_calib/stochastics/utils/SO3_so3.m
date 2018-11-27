function g=SO3_so3(X)
    % exp on skew(so3) or log on SO3
    if (X(3,3)~=0) %If input is SO3, log
        g=logm(X);
    else %If input is so3, exponentiate
        g=expm(X);    
    end
end