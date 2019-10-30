function Rreg = coarse_to_fine_seq_sol(Rdata, Rreg, indices, tau, lambda, miu, N)
    % start from 30
    Ns = 10;
    N0 = Ns;
    Rreg = reshape(Rreg,3,3,[]);
    options = optimoptions('quadprog','MaxIterations',100,'OptimalityTolerance',1e-5,'StepTolerance',1e-5,'Display','off');
    while N0 <= N
        Ncur = round(linspace(1,N,N0));
        for i = 1:length(indices)
            if isempty(find(Ncur == indices(i),1))
                Ncur = [Ncur indices(i)];
            end
        end
        Ncur = sort(Ncur,'ascend');
        indicescur = indices;
        for i = 1:length(indices)
            indicescur(i) = find(Ncur==indices(i),1);
        end

        iter = 1;
        maxiter = 200;
        oldcost = inf;
        tol1 = 1e-6;
%         tr = 1;
        N2 = length(Ncur);
        Rregcur = Rreg(:,:,Ncur);
        Rregcur = reshape(Rregcur,3,[]);
        while iter < maxiter
            newcost = -1e6;
            for j = 1:N2
                id = j;
                xi = data_term_error(Rdata,Rregcur,indicescur,id);
                v = numerical_diff_v(Rregcur,id);
                dxi = seq_sol(xi, v, indicescur, tau, lambda, miu, N2, id, Rregcur, options);
%                 if norm(dxi) > tr
%                     dxi = dxi ./ norm(dxi) .* tr;
%                 end
                Rregcur(:,id*3-2:id*3) = Rregcur(:,id*3-2:id*3) * expSO3(dxi);
                if norm(dxi) > newcost
                    newcost = norm(dxi);
                end
            end
            if abs(newcost - oldcost) < tol1
                break;
            else
                oldcost = newcost;
            end
            iter = iter + 1;
        end
        Rregcur = reshape(Rregcur,3,3,[]);
%         [speed0, acc0] = compute_profiles_fast(tau,Rregcur);
%         figure(1);
%         plot(1:N2,speed0,1:N2,acc0);
        Rreg(:,:,Ncur) = Rregcur;
        if N0 >= N
            break;
        else
            N0 = min(N0+Ns,N);
        end

    end
    Rreg = reshape(Rreg,3,[]);
end
