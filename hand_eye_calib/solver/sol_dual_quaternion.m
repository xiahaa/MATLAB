function varargout = sol_dual_quaternion(TA,TB,N)
%% implementation of hand eye calibration proposed by:
% Konstantinos Daniilidis 
% Hand-Eye Calibration Using Dual Quaternions
% Int. Journ. Robotics Res, 18: 286-298, 1999 

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
        a = qas(2:4);b = qbs(2:4);
        ad = qdas(2:4);bd = qdbs(2:4);
        
        T((i-1)*6+1:(i-1)*6+3,:) = [a-b skewm(a+b) zeros(3,1) zeros(3,3)];
        T((i-1)*6+4:(i-1)*6+6,:) = [ad-bd skewm(ad+bd) a-b skewm(a+b)];
        
%         q1 = rot2vec(T1(1:3,1:3));
%         q2 = rot2vec(T2(1:3,1:3));
%         q1 = q1./norm(q1);
%         q2 = q2./norm(q2);
%         vi = q1;pi = T1(1:3,4);
%         vj = q2;pj = T2(1:3,4);
    end
    %% 2
    [U,S,V] = svd(T);
    s7 = S(7,7);
    s8 = S(8,8);
    if abs(s7) > 0.3 || abs(s8) > 0.3
        error('illness problem due to heavy noises!');
        varargout{1} = [];
        return;
    end
    v7 = V(:,7);
    v8 = V(:,8);
    
    %% 3
    u1 = v7(1:4);v1 = v7(5:8);u2 = v8(1:4);v2 = v8(5:8);
    p2 = u1'*v1; p1 = u1'*v2+u2'*v1; p0 = u2'*v2;
    den = p1*p1 - 4*p0*p2;
    if den < 0
        error('illness problem due to heavy noises!');
        varargout{1} = [];
        return;
    end
    s1 = (-p1 + sqrt(den))/(2*p2);%% 
    s2 = (-p1 - sqrt(den))/(2*p2);
    
    den1 = s1*s1*(u1'*u1)+2*s1*(u1'*u2)+u2'*u2;
    den2 = s2*s2*(u1'*u1)+2*s2*(u1'*u2)+u2'*u2;
    if den1 > den2
        lambda2 = sqrt(1/den1);
        lambda1 = lambda2 * s1;
    else
        lambda2 = sqrt(1/den2);
        lambda1 = lambda2 * s2;
    end
    sol = lambda1.*v7 + lambda2.*v8;
    
    if sol(1)<0
        sol = -sol;
    end
    
    q12 = sol(1:4);
    q12 = q12 ./ norm(q12);%% in case of not unit quaternion
    q12d = sol(5:8);
    
    R12 = q2R(q12);
    t12 = 2.*qprod(q12d, conjugateq(q12));
    t12 = t12(2:4);
   
    if dim == 4
        varargout{1} = [R12 t12;[0 0 0 1]];
    else
        varargout{1} = R12;
    end
end

