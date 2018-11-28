function varargout = sol_cvx2(TA,TB,N)
%% solving hand eye calibration using SCFS

%% Author: xiahaa@space.dtu.dk
    dim = size(TA,2);
    if N < 2
        error('At least two samples needed for unique solution!');
        varargout{1} = [];
        return;
    end
    if dim < 4
        error('Only work for T!');
        varargout{1} = [];
        return;
    end
    
    format short;
    
    options = optimoptions('quadprog','Display','off');
    useSCF = 1;
    useQ = 0;
    
    if useQ == 1
        Nv = 0;
        ids = [];
        thetaas = zeros(N,1);das = zeros(N,1);
        thetabs = zeros(N,1);dbs = zeros(N,1);
        las = zeros(N,3);mas = zeros(N,3);
        lbs = zeros(N,3);mbs = zeros(N,3);
        for i = 1:N
            T1 = TA(i,:,:);T1 = reshape(T1,dim,dim,1);
            T2 = TB(i,:,:);T2 = reshape(T2,dim,dim,1);

            [thetaa, da, la, ma] = screwParams(T2(1:3,1:3),T2(1:3,4));
            [thetab, db, lb, mb] = screwParams(T1(1:3,1:3),T1(1:3,4));

            thetaas(i) = thetaa;das(i) = da;las(i,:) = la';mas(i,:) = ma';
            thetabs(i) = thetab;dbs(i) = db;lbs(i,:) = lb';mbs(i,:) = mb';

            if abs(thetaa-thetab) < 0.15 && abs(da-db) < 0.15
                Nv = Nv + 1;
                ids = [ids;i];
            end
        end
        T = zeros(6*Nv, 8);
        As = zeros(8,8);
        for i = 1:Nv
            ii = ids(i);
            qas = [cos(thetaas(ii)*0.5);sin(thetaas(ii)*0.5).*las(ii,:)'];
            qdas = [(-das(ii)*0.5)*sin(thetaas(ii)*0.5);sin(thetaas(ii)*0.5).*mas(ii,:)'+(das(ii)*0.5*cos(thetaas(ii)*0.5)).*las(ii,:)'];
            qbs = [cos(thetabs(ii)*0.5);sin(thetabs(ii)*0.5).*lbs(ii,:)'];
            qdbs = [(-dbs(ii)*0.5)*sin(thetabs(ii)*0.5);sin(thetabs(ii)*0.5).*mbs(ii,:)'+(dbs(ii)*0.5*cos(thetabs(ii)*0.5)).*lbs(ii,:)'];
            a = qas(2:4);b = qbs(2:4);
            ad = qdas(2:4);bd = qdbs(2:4);

            T((i-1)*6+1:(i-1)*6+3,:) = [a-b skewm(a+b) zeros(3,1) zeros(3,3)];
            T((i-1)*6+4:(i-1)*6+6,:) = [ad-bd skewm(ad+bd) a-b skewm(a+b)];
            
            As = As + T((i-1)*6+1:(i-1)*6+6,:)'*T((i-1)*6+1:(i-1)*6+6,:);
        end
        x0 = [1 0 0 0 0 0 0 0]';
        x = x0;
        if useSCF == 1
            AA = 2*As;
            tic
            for k = 1:1000
                [L,S] = conslin(x);
                soln = quadprog(AA, ft, L, S,[],[],[],[],[],options);
                xnew = soln(1:12);
                if norm(xnew-x) < 1e-12
                    disp(['converge at step ', num2str(k)]);
                    break;
                end
                x = xnew;
            end
            time = toc;
            disp(['scf :',num2str(time)]);
        else
            
        end
    else
        C = zeros(N*12,13);
%         d = zeros(N*12,1);
        As2 = zeros(13,13);
        for i = 1:N
            T1 = TA(i,:,:);T1 = reshape(T1,dim,dim,1);
            T2 = TB(i,:,:);T2 = reshape(T2,dim,dim,1);
            Ra = T2(1:3,1:3);ta = T2(1:3,4);
            Rb = T1(1:3,1:3);tb = T1(1:3,4);

            id = (i-1)*12;
            C(id+1:id+12,:) = [kron(eye(3),Ra)-kron(Rb',eye(3)) zeros(9,3) zeros(9,1); ...
                               kron(tb', eye(3)) eye(3)-Ra -ta];
           As2 = As2 + C(id+1:id+12,:)'*C(id+1:id+12,:);
        end
%         As = C'*C;
        As = (As2+As2')/2;
        %% vector to optimize
    %     x0 = C\d;
        x0 = [1 0 0 0 1 0 0 0 1 0 0 0 1]';
        x = x0;
        
        if useSCF == 1
            regualrization = 0;
            AA = 2*blkdiag(As,regualrization);
            tic
            for k = 1:1000
%                 [L,S, Aeq, beq] = conslin(x);
                [L,S] = conslin(x);
                soln = quadprog(AA, [], L, S, [], [], [], [], [],options);
                xnew = soln(1:13);
                if norm(xnew-x) < 1e-12
                    disp(['converge at step ', num2str(k)]);
                    break;
                end
                x = xnew;
            end
            time = toc;
            disp(['scf :',num2str(time)]);

            tic
            for k = 1:1000
                [L,S, Aeq, beq] = conslin3(x);
%                 [L,S] = conslin(x);
                soln = quadprog(AA, [], L, S, Aeq, beq, [], [], [],options);
                xnew = soln(1:13);
                if norm(xnew-x) < 1e-12
                    disp(['converge at step ', num2str(k)]);
                    break;
                end
                x = xnew;
            end
            time = toc;
            disp(['scf eq:',num2str(time)]);

        else
            AA = 2*As;
            tic
            for k = 1:1000
                [L,S] = conslin2(x);
                soln = quadprog(AA, [], [],[],L,S,[],[],[],options);
                xnew = soln(1:13);
                if norm(xnew-x) < 1e-12
                    disp(['converge at step ', num2str(k)]);
                    break;
                end
                x = xnew;
            end
            time = toc;
            disp(['SQP :',num2str(time)]);
        end
    end

    R12 = [x(1) x(4) x(7); ...
           x(2) x(5) x(8); ...
           x(3) x(6) x(9)];
    t12 = x(10:12);
    
    R12 = R12*inv(sqrtm(R12'*R12));
    
    if dim == 4
        varargout{1} = [R12 t12;[0 0 0 1]];
    else
        varargout{1} = R12;
    end
end


function [L,S] = conslin2(x)
    f1 = @(x) (x(1)*x(4)+x(2)*x(5)+x(3)*x(6));
    f2 = @(x) (x(1)*x(7)+x(2)*x(8)+x(3)*x(9));
    f3 = @(x) (x(4)*x(7)+x(5)*x(8)+x(6)*x(9));
    f4 = @(x) (x(1)*x(1)+x(2)*x(2)+x(3)*x(3));
    f5 = @(x) (x(4)*x(4)+x(5)*x(5)+x(6)*x(6));
    f6 = @(x) (x(7)*x(7)+x(8)*x(8)+x(9)*x(9));
    
    df1 = [x(4) x(5) x(6) x(1) x(2) x(3) 0 0 0 0 0 0 0];
    df2 = [x(7) x(8) x(9) 0 0 0 x(1) x(2) x(3) 0 0 0 0];
    df3 = [0 0 0 x(7) x(8) x(9) x(4) x(5) x(6) 0 0 0 0];
    df4 = [2*x(1) 2*x(2) 2*x(3) 0 0 0 0 0 0 0 0 0 0];
    df5 = [0 0 0 2*x(4) 2*x(5) 2*x(6) 0 0 0 0 0 0 0];
    df6 = [0 0 0 0 0 0 2*x(7) 2*x(8) 2*x(9) 0 0 0 0];
    
    h = [0 0 0 1 1 1]';
    L = [df1;df2;df3;df4;df5;df6];
    b1 = [f1(x)-df1*x;f2(x)-df2*x;f3(x)-df3*x;f4(x)-df4*x;f5(x)-df5*x;f6(x)-df6*x];
    S = h-b1;
    L = [L;zeros(1,12) 1];
    S = [S;1];
end

function [L,S] = conslin(x)
    f1 = @(x) (x(1)*x(4)+x(2)*x(5)+x(3)*x(6));
    f2 = @(x) (x(1)*x(7)+x(2)*x(8)+x(3)*x(9));
    f3 = @(x) (x(4)*x(7)+x(5)*x(8)+x(6)*x(9));
    f4 = @(x) (x(1)*x(1)+x(2)*x(2)+x(3)*x(3));
    f5 = @(x) (x(4)*x(4)+x(5)*x(5)+x(6)*x(6));
    f6 = @(x) (x(7)*x(7)+x(8)*x(8)+x(9)*x(9));
    f7 = @(x) (x(13));

    df1 = [x(4) x(5) x(6) x(1) x(2) x(3) 0 0 0 0 0 0 0];
    df2 = [x(7) x(8) x(9) 0 0 0 x(1) x(2) x(3) 0 0 0 0];
    df3 = [0 0 0 x(7) x(8) x(9) x(4) x(5) x(6) 0 0 0 0];
    df4 = [2*x(1) 2*x(2) 2*x(3) 0 0 0 0 0 0 0 0 0 0];
    df5 = [0 0 0 2*x(4) 2*x(5) 2*x(6) 0 0 0 0 0 0 0];
    df6 = [0 0 0 0 0 0 2*x(7) 2*x(8) 2*x(9) 0 0 0 0];
    df7 = [0 0 0 0 0 0 0 0 0 0 0 0 1];

    h = [0 0 0 1 1 1]';
    a1 = [df1;df2;df3;df4;df5;df6];
    L = [[-a1 -h];[-a1 h]];
    b1 = [f1(x)-df1*x;f2(x)-df2*x;f3(x)-df3*x;f4(x)-df4*x;f5(x)-df5*x;f6(x)-df6*x];
    S = [b1;b1];
%     A = [df7 0];
%     b = 1;
    L = [L;[zeros(1,12) 1 0];[zeros(1,12) -1 0]];
    S = [S;1;-1];
end

function [L,S,A,b] = conslin3(x)
    f1 = @(x) (x(1)*x(4)+x(2)*x(5)+x(3)*x(6));
    f2 = @(x) (x(1)*x(7)+x(2)*x(8)+x(3)*x(9));
    f3 = @(x) (x(4)*x(7)+x(5)*x(8)+x(6)*x(9));
    f4 = @(x) (x(1)*x(1)+x(2)*x(2)+x(3)*x(3));
    f5 = @(x) (x(4)*x(4)+x(5)*x(5)+x(6)*x(6));
    f6 = @(x) (x(7)*x(7)+x(8)*x(8)+x(9)*x(9));

    df1 = [x(4) x(5) x(6) x(1) x(2) x(3) 0 0 0 0 0 0 0];
    df2 = [x(7) x(8) x(9) 0 0 0 x(1) x(2) x(3) 0 0 0 0];
    df3 = [0 0 0 x(7) x(8) x(9) x(4) x(5) x(6) 0 0 0 0];
    df4 = [2*x(1) 2*x(2) 2*x(3) 0 0 0 0 0 0 0 0 0 0];
    df5 = [0 0 0 2*x(4) 2*x(5) 2*x(6) 0 0 0 0 0 0 0];
    df6 = [0 0 0 0 0 0 2*x(7) 2*x(8) 2*x(9) 0 0 0 0];
    df7 = [0 0 0 0 0 0 0 0 0 0 0 0 1];

    h = [0 0 0 1 1 1]';
    a1 = [df1;df2;df3;df4;df5;df6];
    L = [[-a1 -h];[-a1 h]];
    b1 = [f1(x)-df1*x;f2(x)-df2*x;f3(x)-df3*x;f4(x)-df4*x;f5(x)-df5*x;f6(x)-df6*x];
    S = [b1;b1];
    A = [df7 0];
    b = 1;
end