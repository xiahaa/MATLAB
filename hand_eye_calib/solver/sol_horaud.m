function varargout = sol_horaud(TA,TB,N)
%% implementation of hand eye calibration proposed by:
% R. Horaud and F. Dornaika.
% Hand-eye calibration.
% International Journal of Robotics Research, 14(3):195?210, 1995.

%% Author: xiahaa@space.dtu.dk
    if N < 2
        error('At least two samples needed for unique solution!');
        varargout{1} = [];
        return;
    end
    dim = size(TA,2);

    format long;
    At = zeros(4,4);
    for i = 1:N
        T1 = TA(i,:,:);T1 = reshape(T1,dim,dim,1);
        T2 = TB(i,:,:);T2 = reshape(T2,dim,dim,1);
        q1 = rot2vec(T1(1:3,1:3));
        q2 = rot2vec(T2(1:3,1:3));

        if (norm(q1) < 1e-3 || norm(q2) < 1e-3) %mall motion
            continue;
        end

        q1 = q1./norm(q1);
        q2 = q2./norm(q2);
        A = formA([0;q2],[0;q1]);
        At = At + A;
    end
    [U,S,V] = svd(At);
    q12 = V(:,end);
    R12 = q2R(q12);
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

function A = formA(q1,q2)
    A1 = q2m_left(q1);
    A2 = q2m_right(q2);
    A = (A1-A2)'*(A1-A2);
end