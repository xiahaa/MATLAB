function [LHS, RHS] = seq_sol2(Rdata,Rreg,indices,tau,lambda,miu,N,varargin)
    % another way of approximating by expanding around R.
    if isempty(varargin) 
        so3reg = zeros(3,N);
        for i = 1:N
            so3reg(:,i)=logSO3(Rreg(:,i*3-2:i*3));
        end
        so3data = zeros(3,round(size(Rdata,2)/3));
        for i = 1:length(indices)
            so3data(:,i)=logSO3(Rdata(:,i*3-2:i*3));
        end

        lhs = zeros(3,3,N);
        rhs = zeros(3,N);

        % fist deal with term 1
        for i = 1:length(indices)
            Jr = rightJinv(so3reg(:,indices(i)));% * Rreg(:,indices(i)*3-2:indices(i)*3);
            A = (Jr-0.5.*hat(so3data(:,i))*Jr);
            b = so3reg(:,indices(i))-so3data(:,i)-0.5.*hat(so3data(:,i))*so3reg(:,indices(i));
            lhs(:,:,indices(i)) = lhs(:,:,indices(i)) + A'*A;
            rhs(:,indices(i)) = rhs(:,indices(i)) + A'*b;
        end

        c1 = lambda / tau;
        if lambda == 0
        end

        % third term
        c2 = miu / (tau^3);

        % end points
        [AtA,Atb] = calcAb(so3reg(:,1),so3reg(:,2),so3reg(:,3),1);
        lhs(:,:,1) = lhs(:,:,1) + AtA.*c2;
        rhs(:,1) = rhs(:,1) + Atb.*c2;

        [AtA,Atb] = calcAb(so3reg(:,end-2),so3reg(:,end-1),so3reg(:,end),3);
        lhs(:,:,end) = lhs(:,:,end) + AtA.*c2;
        rhs(:,end) = rhs(:,end) + Atb.*c2;

        % 2, two times
        [AtA,Atb] = calcAb(so3reg(:,1),so3reg(:,2),so3reg(:,3),2);
        lhs(:,:,2) = lhs(:,:,2) + AtA.*c2;
        rhs(:,2) = rhs(:,2) + Atb.*c2;

        [AtA,Atb] = calcAb(so3reg(:,2),so3reg(:,3),so3reg(:,4),1);
        lhs(:,:,2) = lhs(:,:,2) + AtA.*c2;
        rhs(:,2) = rhs(:,2) + Atb.*c2;

        % end - 1, two times
        [AtA,Atb] = calcAb(so3reg(:,end-2),so3reg(:,end-1),so3reg(:,end),2);
        lhs(:,:,end-1) = lhs(:,:,end-1) + AtA.*c2;
        rhs(:,end-1) = rhs(:,end-1) + Atb.*c2;

        [AtA,Atb] = calcAb(so3reg(:,end-3),so3reg(:,end-2),so3reg(:,end-1),3);
        lhs(:,:,end-1) = lhs(:,:,end-1) + AtA.*c2;
        rhs(:,end-1) = rhs(:,end-1) + Atb.*c2;

        % 3 times
        for i = 3:N-2
            [AtA,Atb] = calcAb(so3reg(:,i-1),so3reg(:,i),so3reg(:,i+1),2);
            lhs(:,:,i) = lhs(:,:,i) + AtA.*c2;
            rhs(:,i) = rhs(:,i) + Atb.*c2;

            [AtA,Atb] = calcAb(so3reg(:,i),so3reg(:,i+1),so3reg(:,i+2),1);
            lhs(:,:,i) = lhs(:,:,i) + AtA.*c2;
            rhs(:,i) = rhs(:,i) + Atb.*c2;

            [AtA,Atb] = calcAb(so3reg(:,i-2),so3reg(:,i-1),so3reg(:,i),3);
            lhs(:,:,i) = lhs(:,:,i) + AtA.*c2;
            rhs(:,i) = rhs(:,i) + Atb.*c2;
        end

        if c1 == 0 && c2 == 0
            ii = 1:N;
            ii(indices) = [];
            for i = 1:length(ii)
                lhs(:,:,ii(i)) = eye(3);
            end
        end

        LHS = spblkdiag(lhs);
        RHS = rhs(:);
    else
        LHS = zeros(3,3);
        RHS = zeros(3,1);
        
        id = varargin{1};
        so3c = logSO3(Rreg(:,id*3-2:id*3));
        
        mid = find(indices==id,1);
        if ~isempty(mid)
            so3data=logSO3(Rdata(:,mid*3-2:mid*3));
            Jr = rightJinv(so3c);% * Rreg(:,indices(i)*3-2:indices(i)*3);
            A = (Jr-0.5.*hat(so3data)*Jr);
            b = so3c-so3data-0.5.*hat(so3data)*so3c;
            LHS = LHS + A'*A;
            RHS = RHS + A'*b;
        end
        
        c1 = lambda / tau;
        if lambda == 0
        end

        % third term
        c2 = miu / (tau^3);
        
        if id == 1
            so3cc = logSO3(Rreg(:,4:6));
            so3ccc = logSO3(Rreg(:,7:9));
            % end points
            [AtA,Atb] = calcAb(so3c,so3cc,so3ccc,1);
            LHS = LHS + AtA.*c2;
            RHS = RHS + Atb.*c2;
        elseif id == N
            so3cp = logSO3(Rreg(:,end-5:end-3));
            so3cpp = logSO3(Rreg(:,end-8:end-6));
            [AtA,Atb] = calcAb(so3cpp,so3cp,so3c,3);
            LHS = LHS + AtA.*c2;
            RHS = RHS + Atb.*c2;
        elseif id == 2
            so3cc = logSO3(Rreg(:,7:9));
            so3ccc = logSO3(Rreg(:,10:12));
            so3cp = logSO3(Rreg(:,1:3));
            % 2, two times
            [AtA,Atb] = calcAb(so3cp,so3c,so3cc,2);
            LHS = LHS + AtA.*c2;
            RHS = RHS + Atb.*c2;

            [AtA,Atb] = calcAb(so3c,so3cc,so3ccc,1);
            LHS = LHS + AtA.*c2;
            RHS = RHS + Atb.*c2;
        elseif id == N-1
            so3cp = logSO3(Rreg(:,end-8:end-6));
            so3cpp = logSO3(Rreg(:,end-11:end-9));
            so3cc = logSO3(Rreg(:,end-2:end));
            [AtA,Atb] = calcAb(so3cp,so3c,so3cc,2);
            LHS = LHS + AtA.*c2;
            RHS = RHS + Atb.*c2;
            
            [AtA,Atb] = calcAb(so3cpp,so3cp,so3c,3);
            LHS = LHS + AtA.*c2;
            RHS = RHS + Atb.*c2;
        else
            so3cc = logSO3(Rreg(:,(id+1)*3-2:(id+1)*3));
            so3ccc = logSO3(Rreg(:,(id+2)*3-2:(id+2)*3));
            so3cp = logSO3(Rreg(:,(id-1)*3-2:(id-1)*3));
            so3cpp = logSO3(Rreg(:,(id-2)*3-2:(id-2)*3));
            
            [AtA,Atb] = calcAb(so3cp,so3c,so3cc,2);
            LHS = LHS + AtA.*c2;
            RHS = RHS + Atb.*c2;
            
            [AtA,Atb] = calcAb(so3cpp,so3cp,so3c,3);
            LHS = LHS + AtA.*c2;
            RHS = RHS + Atb.*c2;
            
            [AtA,Atb] = calcAb(so3c,so3cc,so3ccc,1);
            LHS = LHS + AtA.*c2;
            RHS = RHS + Atb.*c2;
        end

        if c1 == 0 && c2 == 0
            if isempty(mid)
                LHS = eye(3);
            end
        end
    end
end

function [AtA,Atb] = calcAb(a,b,c,type)
    if type == 1
        bc = logSO3(expSO3(b)'*expSO3(c));
        Jrinv = rightJinv(a);
        s = hat(b);
        A = Jrinv - 0.5.*s*Jrinv;
        b = bc+a-b-0.5.*s*a;
    elseif type == 2
        Jrinv = rightJinv(b);
        sa = hat(a);
        sc = hat(c);
        A = -Jrinv + 0.5.*sa*Jrinv - Jrinv + 0.5.*sc*Jrinv;
        b = c+a-b-b+0.5.*sa*b+0.5.*sc*b;
    else
        ba = logSO3(expSO3(b)'*expSO3(a));
        Jrinv = rightJinv(c);
        s = hat(b);
        A = Jrinv - 0.5.*s*Jrinv;
        b = ba+c-b-0.5.*s*c;
    end
    AtA = A'*A;
    Atb = A'*b;
end

