function varargout = sol_improved_dual_quaternion(TA,TB,N)
%% Robust Hand-Eye Calibration for Computer Aided Medical Endoscopy
% Abed Malti and Joao P. Barreto
% 2010 IEEE International Conference on Robotics and Automation
% May 3-8, 2010, Anchorage, Alaska, USA
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
        
        isvalid = (ma(1) == ma(1)) & (ma(2) == ma(2)) & (ma(3) == ma(3)) ...
                 & (mb(1) == mb(1)) & (mb(2) == mb(2)) & (mb(3) == mb(3));
        
        if abs(thetaa-thetab) < 1e-3 && abs(da-db) < 1e-3 && isvalid
            Nv = Nv + 1;
            ids = [ids;i];
        end
    end
    
    %Construct L and Lp matrix
    L = zeros(4*Nv, 4);
    Lp = zeros(4*Nv, 4);
    
    for i = 1:Nv
        ii = ids(i);
        %% we go either 1 or 2
        %% 1, directly from R,t to dual quaternion, 
%         T1 = TA(ii,:,:);T1 = reshape(T1,dim,dim,1);
%         T2 = TB(ii,:,:);T2 = reshape(T2,dim,dim,1);
%         qa = rot2quat(T1(1:3,1:3));%% need
%         qb = rot2quat(T2(1:3,1:3));%% need
%         qta = [0; T1(1:3,4)];
%         qtb = [0; T2(1:3,4)];
%         qda = qdual(qta,qa);%% need
%         qdb = qdual(qtb,qb);%% need
%         a = qa(2:4);b = qb(2:4);
%         ad = qda(2:4);bd = qdb(2:4);
        
        %% 2, from screw parameters to dual quaternion
        qas = [cos(thetaas(ii)*0.5);sin(thetaas(ii)*0.5).*las(ii,:)'];
        qdas = [(-das(ii)*0.5)*sin(thetaas(ii)*0.5);sin(thetaas(ii)*0.5).*mas(ii,:)'+(das(ii)*0.5*cos(thetaas(ii)*0.5)).*las(ii,:)'];
        qbs = [cos(thetabs(ii)*0.5);sin(thetabs(ii)*0.5).*lbs(ii,:)'];
        qdbs = [(-dbs(ii)*0.5)*sin(thetabs(ii)*0.5);sin(thetabs(ii)*0.5).*mbs(ii,:)'+(dbs(ii)*0.5*cos(thetabs(ii)*0.5)).*lbs(ii,:)'];
        a = qas(1:4);b = qbs(1:4);
        ad = qdas(1:4);bd = qdbs(1:4);
        
        x = a - b;
        y = a + b;
        z = ad - bd;
        w = ad + bd;
        
        L(4*i-3:4*i, :) = [x(1), -x(2:4)';x(2:4), skewm(y(2:4))+x(1)*eye(3)];
        Lp(4*i-3:4*i, :) = [z(1), -z(2:4)';z(2:4), skewm(w(2:4))+z(1)*eye(3)];
    end
    
    %Determine real part and dual part of matrix X
    [~, ~, V] = svd(L);
    q = V(:, 4);
    qprime = L\(-Lp*q);

    q = [q(1);q(2:4)];
    qprime = [qprime(1);qprime(2:4)];
    
    %Determine eye to hand tranformation
    R12 = q2r(q);
    t = 2*qprod(qprime,conjugateq(q));
    t12 = t(2:4);
   
    if dim == 4
        varargout{1} = [R12 t12;[0 0 0 1]];
    else
        varargout{1} = R12;
    end
    
end
