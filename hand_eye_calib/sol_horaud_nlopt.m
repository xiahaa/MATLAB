function varargout = sol_horaud_nlopt(TA,TB,N)
%% implementation of hand eye calibration proposed by:
% R. Horaud and F. Dornaika. 
% Hand-eye calibration. 
% International Journal of Robotics Research, 14(3):195?210, 1995.

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
    
    format long;
    
    l = 2*1e6;
    A = zeros(4,4);
    B = zeros(4,4);
    C = zeros(3,3);
    d = zeros(1,3);
    e = zeros(1,3);
    for i = 1:N
        T1 = TA(i,:,:);T1 = reshape(T1,dim,dim,1);
        T2 = TB(i,:,:);T2 = reshape(T2,dim,dim,1);
        q1 = rot2vec(T1(1:3,1:3));
        q2 = rot2vec(T2(1:3,1:3));
        
        vi = q1;pi = T1(1:3,4);
        vj = q2;pj = T2(1:3,4);
        %% A
        A1 = formA([0;q2],[0;q1]);
        A = A + A1;
        %% B
        B1 = (pi'*pi+pj'*pj).*eye(4) - ...
              q2m_right([0;pi])'*q2m_left([0;pj]) - ...
              q2m_left([0;pj])'*q2m_right([0;pi]);
        B = B+B1;
        %% C
        K = T2(1:3,1:3);
        C = C + (K'*K-K-K'+eye(3));
        %% d and e
        d = d + 2.*pj'*(K-eye(3));
        e = e - 2.*pi'*(T1(1:3,1:3)-eye(3));
    end
    
    [U,S,V] = svd(A);
    q12 = V(:,end);
    R12 = q2R(q12);
    A1 = zeros(3*N,3);
    b1 = zeros(3*N,1);
    for i = 1:N
        T1 = TA(i,:,:);T1 = reshape(T1,dim,dim,1);
        T2 = TB(i,:,:);T2 = reshape(T2,dim,dim,1);
        A1((i-1)*3+1:i*3,:) = -T2(1:3,1:3)+eye(3);
        b1((i-1)*3+1:i*3) = -R12*T1(1:3,4)+T2(1:3,4);
    end
    t12 = A1\b1;
    
    x0 = [q12;t12];
%     options = optimoptions('fmincon','Algorithm','Trust-region-reflective'); % run Trust-region-reflective algorithm
    x = fmincon(@(x) fcost(x,A,B,C,d,e,l),x0);
    q12 = x(1:4);
    q12 = q12./norm(q12);
    R12 = q2R(q12);
    if dim == 4
        t12 = x(5:7);
        varargout{1} = [R12 t12;[0 0 0 1]];
    else
        varargout{1} = R12;
    end
end

function cost = fcost(x,A,B,C,d,e,l)
    q = x(1:4);
    t = x(5:7);
    cost = q'*(A+B)*q+t'*C*t+d*t+e*(q2R(q))'*t+l*(1-q'*q)^2;
end

function A = formA(q1,q2)
    A1 = q2m_left(q1);
    A2 = q2m_right(q2);
    A = (A1-A2)'*(A1-A2);
end