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
    
    format short;
    
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
        q1 = q1./norm(q1);
        q2 = q2./norm(q2);
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
    
%     [f1,g1] = fcost(x0,A,B,C,d,e,l);
%     f2 = zeros(3,1);
%     q = x0(1:4);
%     t = x0(5:7);
%     for i = 1:N
%         T1 = TA(i,:,:);T1 = reshape(T1,dim,dim,1);
%         T2 = TB(i,:,:);T2 = reshape(T2,dim,dim,1);
%         vi = rot2vec(T1(1:3,1:3));
%         vj = rot2vec(T2(1:3,1:3));
%         pi = T1(1:3,4);
%         pj = T2(1:3,4);
%         vi = vi./norm(vi);
%         vj = vj./norm(vj);
%   
%         ff1 = norm([0;vj]-prod_qpqconj(q,[0;vi]))^2;
%         A1 = formA([0;vj],[0;vi]);
%         ff11 = q'*A1*q;
%         t1 =  prod_qpqconj(q,[0;pi]);
%         K = T2(1:3,1:3);
%         ff2 = norm(t1(2:4) - (K-eye(3))*t - pj)^2;
%         f2(1) = f2(1) + ff1;
%         f2(2) = f2(2) + ff2;
%     end
%     f2(3) = l*(1-q'*q)^2;
%     f2
    
    options = optimoptions('fminunc','Algorithm','trust-region', ...
        'SpecifyObjectiveGradient',true,'MaxIterations',1000,'OptimalityTolerance',1e-6);
    x = fminunc(@(x) fcost(x,A,B,C,d,e,l),x0,options);
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

function [f,g] = fcost(x,A,B,C,d,e,l)
    q = x(1:4);
    t = x(5:7);
    f1 = q'*(A+B)*q+t'*C*t+d*t+e*(q2R(q))'*t+l*(1-q'*q)^2;
       
%     disp(strcat('1:',num2str(q'*(A)*q)));
%     disp(strcat('2:',num2str(q'*B*q)));
%     disp(strcat('3:',num2str(t'*C*t)));
%     disp(strcat('4:',num2str(d*t+e*(q2R(q))'*t)));
%     disp(strcat('5:',num2str(l*(1-q'*q)^2)));
    f = 0.5*f1'*f1;

    if nargout > 1
        %% provide gradient
        g = f1.*[2.*(A+B)*q + grad1(e,q,t) + l.*2.*(1-q'*q)*(-2.*q); ...
             2.*C*t + d' + (q2R(q))*e'];
    end
end

function g1 = grad1(e,q,t)
    t1 = t(1);t2 = t(2);t3 = t(3);
    e1 = e(1);e2 = e(2);e3 = e(3);
    q1 = q(1);q2 = q(2);q3 = q(3);q4 = q(4);
    
    e1q1 = e1*q1;e1q2 = e1*q2;e1q3 = e1*q3;e1q4 = e1*q4;
    e2q1 = e2*q1;e2q2 = e2*q2;e2q3 = e2*q3;e2q4 = e2*q4; 
    e3q1 = e3*q1;e3q2 = e3*q2;e3q3 = e3*q3;e3q4 = e3*q4;
    
%     t1*(2*e1*q1 - 2*e2*q4 + 2*e3*q3) + t2*(2*e2*q1 + 2*e1*q4 - 2*e3*q2) + t3*(2*e2*q2 - 2*e1*q3 + 2*e3*q1)
%     t1*(2*e1*q2 + 2*e2*q3 + 2*e3*q4) - t2*(2*e2*q2 - 2*e1*q3 + 2*e3*q1) + t3*(2*e2*q1 + 2*e1*q4 - 2*e3*q2)
%     t1*(2*e2*q2 - 2*e1*q3 + 2*e3*q1) + t2*(2*e1*q2 + 2*e2*q3 + 2*e3*q4) - t3*(2*e1*q1 - 2*e2*q4 + 2*e3*q3)
%     t2*(2*e1*q1 - 2*e2*q4 + 2*e3*q3) - t1*(2*e2*q1 + 2*e1*q4 - 2*e3*q2) + t3*(2*e1*q2 + 2*e2*q3 + 2*e3*q4)
%  
%     
    g1 = zeros(4,1);
    g1(1) = t1*(e1q1 - e2q4 + e3q3)*2 + t2*(e2q1 + e1q4 - e3q2)*2 + t3*(e2q2 - e1q3 + e3q1)*2;
    g1(2) = t1*(e1q2 + e2q3 + e3q4)*2 - t2*(e2q2 - e1q3 + e3q1)*2 + t3*(e2q1 + e1q4 - e3q2)*2;
    g1(3) = t1*(e2q2 - e1q3 + e3q1)*2 + t2*(e1q2 + e2q3 + e3q4)*2 - t3*(e1q1 - e2q4 + e3q3)*2;
    g1(4) = t2*(e1q1 - e2q4 + e3q3)*2 - t1*(e2q1 + e1q4 - e3q2)*2 + t3*(e1q2 + e2q3 + e3q4)*2;
 
end

function A = formA(q1,q2)
    A1 = q2m_left(q1);
    A2 = q2m_right(q2);
    A = (A1-A2)'*(A1-A2);
end