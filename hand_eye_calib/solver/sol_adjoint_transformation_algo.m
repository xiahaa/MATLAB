function varargout = sol_adjoint_transformation_algo(TA,TB,N)
%% Adjoint Transformation Algorithm for Hand?Eye Calibration with Applications in Robotic Assisted Surgery
% KRITTIN PACHTRACHAI, FRANCISCO VASCONCELOS, FRANCOIS CHADEBECQ, MAX ALLAN,
% STEPHEN HAILES, VIJAY PAWAR, and DANAIL STOYANOV
% Annals of Biomedical Engineering, Vol. 46, No. 10, October 2018
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
    
    converge_thres = 10;
    counter_rotation = 0;
    counter_translation = 0;
    max_ite = 50;
    v_a = zeros(3, Nv);
    v_b = zeros(3, Nv);
    om_a = zeros(3, Nv);
    om_b = zeros(3, Nv);
    K = zeros(8*Nv, 4);

    X = eye(4);
    
    LHS = zeros(3*Nv, 3);
    RHS = zeros(3*Nv, 1);
    
    for j = 1:Nv
        ii = ids(j);
        T1 = TA(ii,:,:);T2 = TB(ii,:,:);
        A = reshape(T2,dim,dim,1);
        B = reshape(T1,dim,dim,1);
        
%         a = logm(A(:, :));
%         b = logm(B(:, :));
%         v_a(:, j) = a(1:3 ,4);
%         v_b(:, j) = b(1:3, 4);
%         om_a(:, j) = rot2vec(A(1:3, 1:3));
%         om_b(:, j) = rot2vec(B(1:3, 1:3));

        a = tran2vec(A);
        b = tran2vec(B);
        v_a(:, j) = a(1:3);
        v_b(:, j) = b(1:3);
        om_a(:, j) = a(4:6);
        om_b(:, j) = b(4:6);
        
        qa = rot2quat(A(1:3,1:3));%% need
        qb = rot2quat(B(1:3,1:3));%% need
%         qta = [0; A(1:3,4)];
%         qtb = [0; B(1:3,4)];
%         ap = qdual(qta,a);%% need
%         bp = qdual(qtb,b);%% need
        
        x = qa - qb;
        y = qa + qb;

        K(8*j - 7:8*j - 4, 1:4) = [x(1), -(x(2:4))'; x(2:4), skewm(y(2:4)) + x(1)*eye(3)];
    end
    
    counter = 0;
    eps = 1e-6;
    while (((counter_rotation < converge_thres) || (counter_translation < converge_thres)) && (counter < max_ite))
        R_init = X(1:3, 1:3);
        t_init = X(1:3, 4);
        for j = 1:Nv
            vec_buf = v_a(:, j) - cross(t_init, om_a(:, j));
            K(8*j - 3:8*j, 1) = [vec_buf - v_b(:, j); 0];
            K(8*j - 3:8*j, 2:4) = [skewm(vec_buf + v_b(:, j)); (-vec_buf + v_b(:, j))'];
        end
        
        [~, ~, v_basis] = svd(K);
        v_basis = v_basis(:, 4);

        qR = [v_basis(1);v_basis(2:4)];
        
        X(1:3, 1:3) = q2r(qR);

        for j = 1:Nv
            LHS(3*j - 2:3*j, :) = skewm(om_a(:, j));
            RHS(3*j - 2:3*j, 1) = X(1:3, 1:3)*v_b(:, j) - v_a(:, j);
        end
        X = [X(1:3, 1:3) LHS\RHS;0 0 0 1];
        
        diff = X/[R_init, t_init;0 0 0 1];
        
        if (norm(rot2vec(diff(1:3, 1:3))) < eps)
            counter_rotation = counter_rotation + 1;
        else
            counter_rotation = 0;
        end
        
        if (norm(diff(1:3, 4)) < eps)
            counter_translation = counter_translation + 1;
        else
            counter_translation = 0;
        end
        
        counter = counter + 1;
%         ccc = rodrigues(X(1:3, 1:3));
%         C(1:3, counter) = ccc;
%         C(4:6, counter) = X(1:3, 4);
    end
    
    %Determine eye to hand tranformation
    R12 = X(1:3,1:3);
    t12 = X(1:3,4);
   
    if dim == 4
        varargout{1} = [R12 t12;[0 0 0 1]];
    else
        varargout{1} = R12;
    end
    
end