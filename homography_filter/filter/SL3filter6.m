function Hest = SL3filter6(Hest, CC, CCmeas, dt)
    Kk = eye(3)*1;
    maxiter = 1000;
    for iter = 1:maxiter
        invHest = inv(Hest);
        Corr = zeros(3,3);
        for i = 1:size(CCmeas,3)
            ek = invHest'*CCmeas(:,:,i)*invHest;
            Corr = Corr + ek*(ek - CC(:,:,i))*Kk+ek*Kk*(ek-CC(:,:,i));
        end
        Delta = tosl3(Corr);
        Delta = invHest * Delta * Hest;
        %%%%%%%%%%%%%%%%%% Armijo line search %%%%%%%%%%%%%%%%
        alpha = 0.05;
        sigma = 0.25;
        beta = 0.75;
        tk = 0;
        fprintf('start Armijo ... ');
        m = 0;
        temp_fun_cost = inf;
        tangent_inner = trace(Delta*Delta');
        old_fun_cost = cost(Hest, CC, CCmeas, Kk);
        while (old_fun_cost - temp_fun_cost) <= tk*sigma*tangent_inner
            tk = beta^m * alpha;
            tempH = Hest * expm(Delta.*tk);
            temp_fun_cost = cost(tempH, CC, CCmeas, Kk);
            m = m+1;
            if m > 100
                break;
            end
        end
        fprintf('end Armijo at %d\n',m);
        Hnew = Hest * expm(Delta.*tk);
        if norm(Hest-Hnew,'fro') < 1e-6
            break;
        end
        Hest = Hnew;
    end
end

function c = cost(H, CC, CCm, Kk)
    c = 0;
    invH = inv(H);
    for i = 1:size(CCm,3)
        ek = invH'*CCm(:,:,i)*invH - CC(:,:,i);
        c = c + trace(ek*Kk*ek');
    end
    c = c * 0.5;
end