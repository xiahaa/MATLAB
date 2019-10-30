function Rreg = seg2seg_seq_sol(Rdata, Rreg, indices, tau, lambda, miu, N)
    Rreg = reshape(Rreg,3,3,[]);
    options = optimoptions('quadprog','MaxIterations',100,'OptimalityTolerance',1e-5,'StepTolerance',1e-5,'Display','off');
    cost2 = inf;
    for k = 1:10
        for i = 2:length(indices)-1
            Ncur = indices(i-1):indices(i+1);
            indicescur = [];
            for ii = 1:length(indices)
                indicescur = [indicescur find(Ncur==indices(ii),1)];
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
                stidx = 1;
                if i > 2
                    stidx = N2 - (indices(i+1)-indices(i)) + 1;
                end
                for j = stidx:N2
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
            [speed0, acc0] = compute_profiles_fast(tau,Rregcur);
            figure(1);
            plot(1:N2,speed0,1:N2,acc0);
            Rreg(:,:,Ncur) = Rregcur;
        end
        xi = data_term_error(Rdata,Rreg,indices);
        v = numerical_diff_v(Rreg);
        cost1 = cost(xi,v,tau,lambda,miu);
        if abs(cost1 - cost2) < 1e-3
            break;
        else
            cost2 = cost1;
        end
    end
    Rreg = reshape(Rreg,3,[]);
end
