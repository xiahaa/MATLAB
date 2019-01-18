function [Rres,tres] = posit_algo(P, q, K)
% Implementation of the POSIT algorithm described in 
% Iterative Pose Estimation Using Coplanar Feature Points,
% DENIS OBERKAMPF, DANIEL F. DEMENTHON, AND LARRY S. DAVIS
% COMPUTER VISION AND IMAGE UNDERSTANDING
% Vol. 63, No. 3, May, pp. 495?511, 1996
%
% Author: xiahaa@space.dtu.dk

    %% 
    n = size(q,2);
    if size(q,1) == 2
        q = [q;ones(1,n)];% to homogeneous
    end
    %% normalize
    qn = K\q;

    %% precompute
    M0Mi = P(:,2:end) - repmat(P(:,1),1,n-1);
    xix0 = qn(1,2:end) - repmat(qn(1,1),1,n-1);
    yiy0 = qn(2,2:end) - repmat(qn(2,1),1,n-1);
    xi = qn(1,2:end);
    yi = qn(2,2:end);
    
    [U,S,V] = svd(M0Mi');
%         I0 = A\b1;
%         J0 = A\b2;
%     AA = V*inv(S'*S)*S'*U';
%     AA = pinv(M0Mi');
    validS = abs(diag(S)) > 1e-6;
    Sinv = zeros(3,n-1);
    Sinv(1:2,1:2) = diag(1./diag(S(validS,validS)));
    AA = V*Sinv*U';

    
    eps_old = zeros(2,n-1);
    eps_new = zeros(2,n-1);
    
    Ropt = zeros(3,3,2);
    topt = zeros(3,1,2);
    
    current_trace = 1;
    
    iter = 1;
    maxiter = 100;
    
    oldE = 1e6;
    
    while iter < maxiter
        Isucc = zeros(3,4);
        Jsucc = zeros(3,4);
        flag = zeros(1,4);
        Rs = zeros(3,3,4);
        ts = zeros(3,1,4);
        succcnt = 0;
        idsucc = [-1 -1 -1 -1];
        for i = 1:current_trace
            [I, J] = solveIJ(AA, V, xix0, yiy0, eps_new(i,:), xi, yi);
            for j = 1:size(I,2)
                [R,t] = decompose(I(:,j),J(:,j),qn(1:2,1));
                infront = check_res(R,t,M0Mi);
                
                if infront == 1
                    flag(1,j+(i-1)*2) = 1;
                    Isucc(:,j+(i-1)*2) = I(:,j);
                    Jsucc(:,j+(i-1)*2) = J(:,j);
                    Rs(:,:,j+(i-1)*2) = R;
                    ts(:,:,j+(i-1)*2) = t;                    
                    succcnt = succcnt + 1;
                    idsucc(succcnt) = j+(i-1)*2;
                end
            end
        end
        
        if current_trace == 1
            if succcnt == 0
                idb1 = -1; idb2 = -1;
            else
                idb1 = idsucc(1);
                idb2 = idsucc(2);
            end
        elseif current_trace == 2
            if succcnt == 0
                idb1 = -1; idb2 = -1;
            elseif succcnt <= 2
                idb1 = idsucc(1);
                idb2 = idsucc(2);
            else
                %% filter
                idb1 = -1;
                %% branch 1
                if flag(1,1) == 1 && flag(1,2) == 1
                    %% select only 1
                    err1 = calc_err(Rs(:,:,1),ts(:,:,1),P,q);
                    err2 = calc_err(Rs(:,:,2),ts(:,:,2),P,q);
                    if err1 > err2
                        idb1 = 2;
                    else
                        idb1 = 1;
                    end
                elseif flag(1,1) == 1
                    idb1 = 1;
                elseif flag(1,2) == 1
                    idb1 = 2;
                end
                idb2 = -1;
                %% branch 2
                if flag(1,3) == 1 && flag(1,4) == 1
                    %% select only 1
                    err3 = calc_err(Rs(:,:,3),ts(:,:,3),P,q);
                    err4 = calc_err(Rs(:,:,4),ts(:,:,4),P,q);
                    if err3 > err4
                        idb2 = 4;
                    else
                        idb2 = 3;
                    end
                elseif flag(1,3) == 1
                    idb2 = 3;
                elseif flag(1,4) == 1
                    idb2 = 4;
                end
            end
        end
        
        if idb2 == -1 && idb1 == -1
            break;
        elseif idb2 ~= -1 && idb1 ~= -1
            eps_old = eps_new;
            current_trace = 2;
            
            invZ01 = norm(Isucc(:,idb1));
            k = Rs(3,:,idb1)';
            eps_new(1,:) = M0Mi'*k.*invZ01;
            
            invZ02 = norm(Isucc(:,idb2));
            k = Rs(3,:,idb2)';
            eps_new(2,:) = M0Mi'*k.*invZ02;
            
            Ropt(:,:,1) = Rs(:,:,idb1);
            Ropt(:,:,2) = Rs(:,:,idb2);
            topt(:,:,1) = ts(:,:,idb1);
            topt(:,:,2) = ts(:,:,idb2);
            
            err1 = calc_err(Rs(:,:,idb1),ts(:,:,idb1),P,q);
            err2 = calc_err(Rs(:,:,idb2),ts(:,:,idb2),P,q);
            E = min(err1,err2);
            
        elseif idb2 == -1 && idb1 ~= -1
            eps_old = eps_new;
            current_trace = 1;
            invZ01 = norm(Isucc(:,idb1));
            k = Rs(3,:,idb1)';
            eps_new(1,:) = M0Mi'*k.*invZ01;
            Ropt(:,:,1) = Rs(:,:,idb1);
            topt(:,:,1) = ts(:,:,idb1);
            err1 = calc_err(Rs(:,:,idb1),ts(:,:,idb1),P,q);
            E = err1;
        else
            eps_old = eps_new;
            current_trace = 1;
            invZ02 = norm(Isucc(:,idb2));
            k = Rs(3,:,idb2)';
            eps_new(1,:) = M0Mi'*k.*invZ02;
            Ropt(:,:,1) = Rs(:,:,idb2);
            topt(:,:,1) = ts(:,:,idb2);
            err1 = calc_err(Rs(:,:,idb2),ts(:,:,idb2),P,q);
            E = err1;
        end
        
        if E > oldE
            break;
        else
            oldE = E;
        end
        
        err = max(norm(eps_new(1,:)-eps_old(1,:)),norm(eps_new(2,:)-eps_old(2,:)));
        if err < 1e-6
            break;
        end
        iter = iter + 1;
    end
    
    if current_trace == 1
        Rres = Ropt(:,:,1);tres = topt(:,:,1);
    else
        Rres = Ropt;tres = topt;
    end
    
end

function err = calc_err(R,t,P,q)
    pp = R*P + repmat(t,1,size(P,2));
    err = 0;
    for i = 1:size(P,2)
        pp1 = pp(:,i)./norm(pp(:,i));
        pp1 = pp1./pp1(3);
        err = err + norm(pp1-q(:,i));
    end
end

function infront = check_res(R,t,P)
    Z = R(3,:)*P + repmat(t(3),1,size(P,2));
    id = find(Z<0,1);
    if isempty(id)
        infront = 1;
    else
        infront = 0;
    end
end

function [R,t] = decompose(I,J,q0)
    Z0 = 1./norm(I);
    X0 = q0(1)*Z0;
    Y0 = q0(2)*Z0;
    r1t = (I)'./norm(I);
    r2t = (J)'./norm(J);
    r3t = cross(r1t',r2t')';
    r3t = r3t./norm(r3t);
    R = [r1t;r2t;r3t];
    t = [X0;Y0;Z0];
end

function [I, J] = solveIJ(AA, V, xix0, yiy0, eps, xi, yi)
%     A = M0Mi';% nx3
    
    b1 = xix0 + xi.*eps;
    b2 = yiy0 + yi.*eps;
    b1 = b1';
    b2 = b2';
    
    %% check svd
    if rank(AA) == 3
        I = AA*b1;
        J = AA*b2;
    else
%         [U,S,V] = svd(A);
% %         I0 = A\b1;
% %         J0 = A\b2;
%         CA = V*inv(S'*S)*S'*U';
        I0 = AA*b1;
        J0 = AA*b2;
        
        u = V(:,end);
        u = u ./ norm(u);
        %% solve for lambda, miu
        I02 = norm(I0)^2;
        J02 = norm(J0)^2;
        I0J0 = I0'*J0;
        diff = -I02 + J02;
        
        if abs(diff) < 1e-8
            THETA = -sign(I0J0)*pi/2;
            R = abs(I0J0);
        else
            R = sqrt(diff^2+4*I0J0^2);
            if diff > 0
                THETA = atan2(-2*I0J0, diff);
            else
                THETA = atan2(-2*I0J0, diff) + pi;
            end
        end
        rho = sqrt(R);
        theta = THETA / 2;
        
        rhost = rho*sin(theta);
        rhoct = rho*cos(theta);
        
        I(:,1) = I0 + rhoct*u;
        I(:,2) = I0 - rhoct*u;
        
        J(:,1) = J0 + rhost*u;
        J(:,2) = J0 - rhost*u;
    end
end