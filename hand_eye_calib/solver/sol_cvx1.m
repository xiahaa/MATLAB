function varargout = sol_cvx1(TA,TB,N)
%% implementation of hand eye calibration proposed by:
% Zhao, Zijian. "Hand-eye calibration using convex optimization." Robotics and Automation (ICRA), 2011 IEEE International Conference on. IEEE, 2011.

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
    
    C = zeros(N*12,12);
    d = zeros(N*12,1);
    for i = 1:N
        T1 = TA(i,:,:);T1 = reshape(T1,dim,dim,1);
        T2 = TB(i,:,:);T2 = reshape(T2,dim,dim,1);
        Ra = T2(1:3,1:3);ta = T2(1:3,4);
        Rb = T1(1:3,1:3);tb = T1(1:3,4);
        
        id = (i-1)*12;
        C(id+1:id+12,:) = [eye(9)-kron(Ra,Rb) zeros(9,3); ...
                           kron(eye(3),tb') eye(3)-Ra];
        d(id+1:id+12) = [zeros(9,1);ta];
    end
    
    %% formulation proposed SOCP
    cvx_begin quiet
        variable t
        variable x(12)
        minimize (t)
        subject to:
            for i = 1:N
                id = (i-1)*12;
                C1 = C(id+1:id+12,:);
                d1 = d(id+1:id+12);
                norm(C1*x-d1,2) <= t
            end
    cvx_end
    
    %% QP
%     As = zeros(12,12);
%     bs = zeros(1,12);
%     for i = 1:N
%         id = (i-1)*12;
%         
%         As = As + C(id+1:id+12,:)'*C(id+1:id+12,:);
%         bs = d(id+1:id+12)'*C(id+1:id+12,:);
%     end
%     A1 = zeros(12,12);A1(1,2)=1;A1(4,5)=1;A1(7,8)=1;
%     A2 = zeros(12,12);A2(1,3)=1;A2(4,6)=1;A2(7,9)=1;
%     A3 = zeros(12,12);A3(2,3)=1;A3(5,6)=1;A3(8,9)=1;
%     A4 = zeros(12,12);A4(1,1)=1;A4(4,4)=1;A4(7,7)=1;
%     A5 = zeros(12,12);A5(2,2)=1;A5(5,5)=1;A5(8,8)=1;
%     A6 = zeros(12,12);A5(3,3)=1;A5(6,6)=1;A5(9,9)=1;
%     cvx_begin quiet
%         variable x(12)
%         minimize (x'*As*x-2*bs*x)
%         subject to:
%             x'*A1*x == 0
%             x'*A2*x == 0
%             x'*A3*x == 0
%             x'*A4*x == 1
%             x'*A5*x == 1
%             x'*A6*x == 1
%     cvx_end
    
    R12 = [x(1) x(2) x(3); ...
           x(4) x(5) x(6); ...
           x(7) x(8) x(9)];
    t12 = x(10:12);
    
    R12 = R12*inv(sqrtm(R12'*R12));
    
    if dim == 4
        varargout{1} = [R12 t12;[0 0 0 1]];
    else
        varargout{1} = R12;
    end
end