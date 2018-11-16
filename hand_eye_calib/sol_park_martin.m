function varargout = sol_park_martin(TA,TB,N)
%% implementation of hand eye calibration proposed by:
%   F. Park and B. Martin 
%   Robot Sensor Calibration: Solving AX = XB on the Euclidean Group, 
%   IEEE Transactions on Robotics and Automation, (10) 1994, p. 717?721

%% Author: xiahaa@space.dtu.dk
    if N < 2
        error('At least two samples needed for unique solution!');
        varargout{1} = [];
        return;
    end
    dim = size(TA,2);
    M = zeros(3,3);
    
    format long;
    
    if N == 2
        T1 = TA(1,:,:);T1 = reshape(T1,dim,dim,1);
        T2 = TB(1,:,:);T2 = reshape(T2,dim,dim,1);
        T3 = TA(2,:,:);T3 = reshape(T3,dim,dim,1);
        T4 = TB(2,:,:);T4 = reshape(T4,dim,dim,1);
        so31 = rot2vec(T1(1:3,1:3));
        so32 = rot2vec(T2(1:3,1:3));
        so33 = rot2vec(T3(1:3,1:3));
        so34 = rot2vec(T4(1:3,1:3));
        R12 = [so31,so33,cross(so31,so33)]*inv([so32,so34,cross(so32,so34)]);
    else
        for i = 1:N
            T1 = TA(i,:,:);T1 = reshape(T1,dim,dim,1);
            T2 = TB(i,:,:);T2 = reshape(T2,dim,dim,1);
            %% log
            so31 = rot2vec(T1(1:3,1:3));
            so32 = rot2vec(T2(1:3,1:3));
            M = M + so31*so32';
        end
        R12 = inv(sqrtm(M'*M))*M';
    end
    
    if dim == 4
        A = zeros(3*N,3);
        b = zeros(3*N,1);
        for i = 1:N
            T1 = TA(i,:,:);T1 = reshape(T1,dim,dim,1);
            T2 = TB(i,:,:);T2 = reshape(T2,dim,dim,1);
            
            A((i-1)*3+1:i*3,:) = -T2(1:3,1:3)+eye(3);
            b((i-1)*3+1:i*3) = -R12*T1(1:3,4)+T2(1:3,4);
        end
        t12 = A\b;
        varargout{1} = [R12 t12;[0 0 0 1]];
    else
        varargout{1} = R12;
    end
end
