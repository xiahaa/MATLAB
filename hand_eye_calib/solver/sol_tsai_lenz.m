function varargout = sol_tsai_lenz(TA,TB,N)
%% implementation of hand eye calibration proposed by:
%   Shiu, Y.C., Ahmad, S. 
%   A New Technique for Fully Autonomous and Efficient 3D Robotics Hand/Eye Calibration
%   IEEE Transactions onRobotics and Automation, 5 (3) 1989 p.345-358

%% Author: xiahaa@space.dtu.dk
    if N < 2
        error('At least two samples needed for unique solution!');
        varargout{1} = [];
        return;
    end
    dim = size(TA,2);
    A = zeros(3*N,3);
    b = zeros(3*N,1);
    for i = 1:N
        T1 = TA(i,:,:);T1 = reshape(T1,dim,dim,1);
        T2 = TB(i,:,:);T2 = reshape(T2,dim,dim,1);
        pr1 = R2pr(T1(1:3,1:3));
        pr2 = R2pr(T2(1:3,1:3));
        A((i-1)*3+1:i*3,:) = skewm(pr1+pr2);
        b((i-1)*3+1:i*3) = pr1-pr2;
    end
    y = A\b;
    p12 = 2.*y./sqrt(1+y'*y);
    R12 = pr2R(p12);
    
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

function pr = R2pr(R)
    format long
    pr = rot2vec(R);
    theta = norm(pr);
    pr = pr ./ theta .* (2*sin(theta*0.5));
%     R'*pr2R(pr)-eye(3)
end

function R = pr2R(pr)
    format long
    pr2 = pr'*pr;
    R = (1-pr2*0.5)*eye(3) + 0.5.*(pr*pr'+sqrt(4-pr2).*skewm(pr));
end