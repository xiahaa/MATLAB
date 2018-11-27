function varargout = sol_shiu_admad(TA,TB,N)
%% implementation of hand eye calibration proposed by:
%   Shiu, Y.C., Ahmad, S. 
%   Calibration of wrist-mounted robotic sensors by solving homogeneous transform equations of the form AX=XB 
%   IEEE Transactions on Robotics and Automation, 5 (1) 1989, p.16-29

%% Author: xiahaa@space.dtu.dk
    if N < 2
        error('At least two samples needed for unique solution!');
        varargout{1} = [];
        return;
    end
    
    dim = size(TA,2);
    
    TA1 = TA(1,:,:);TA1 = reshape(TA1,dim,dim,1);
    TA2 = TA(2,:,:);TA2 = reshape(TA2,dim,dim,1);
    TB1 = TB(1,:,:);TB1 = reshape(TB1,dim,dim,1);
    TB2 = TB(2,:,:);TB2 = reshape(TB2,dim,dim,1);

    
    Ra1 = TA1(1:3,1:3);
    Ra2 = TA2(1:3,1:3);
    Rb1 = TB1(1:3,1:3);
    Rb2 = TB2(1:3,1:3);
    
    ka1 = rot2vec(Ra1);ka2 = rot2vec(Ra2);
    kb1 = rot2vec(Rb1);kb2 = rot2vec(Rb2);

    ka1 = ka1 ./ norm(ka1);kb1 = kb1 ./ norm(kb1);
    ka2 = ka2 ./ norm(ka2);kb2 = kb2 ./ norm(kb2);

    v1 = cross(kb1,ka1);
    w1 = atan2(norm(v1),dot(kb1,ka1));
    v1 = v1./norm(v1);
    R1 = vec2rot(w1.*v1);
    
    v2 = cross(kb2,ka2);
    w2 = atan2(norm(v2),dot(kb2,ka2));
    v2 = v2./norm(v2);
    R2 = vec2rot(w2.*v2);
    
    %% (42-44)
    n1 = Ra1(:,1);o1 = Ra1(:,2);a1 = Ra1(:,3);
    n2 = Ra2(:,1);o2 = Ra2(:,2);a2 = Ra2(:,3);
    
    nx1 = R1(:,1);ox1 = R1(:,2);ax1 = R1(:,3);
    nx2 = R2(:,1);ox2 = R2(:,2);ax2 = R2(:,3);
    kx1 = rot2vec(R1);kx2 = rot2vec(R2);
    kx1 = kx1 ./ norm(kx1);kx2 = kx2 ./ norm(kx2);
    
    %% temporal values
    dot_n1_ka1 = dot(n1,ka1);cross_n1_ka1 = cross(n1,ka1);
    dot_o1_ka1 = dot(o1,ka1);cross_o1_ka1 = cross(o1,ka1);
    dot_a1_ka1 = dot(a1,ka1);cross_a1_ka1 = cross(a1,ka1);

    dot_n2_ka2 = dot(n2,ka2);cross_n2_ka2 = cross(-n2,ka2);
    dot_o2_ka2 = dot(o2,ka2);cross_o2_ka2 = cross(-o2,ka2);
    dot_a2_ka2 = dot(a2,ka2);cross_a2_ka2 = cross(-a2,ka2);
    
    A = [-nx1(1)+kx1(1)*dot_n1_ka1 cross_n1_ka1(1) nx2(1)-kx2(1)*dot_n2_ka2 cross_n2_ka2(1); ...
         -ox1(1)+kx1(1)*dot_o1_ka1 cross_o1_ka1(1) ox2(1)-kx2(1)*dot_o2_ka2 cross_o2_ka2(1); ...
         -ax1(1)+kx1(1)*dot_a1_ka1 cross_a1_ka1(1) ax2(1)-kx2(1)*dot_a2_ka2 cross_a2_ka2(1); ...
         -nx1(2)+kx1(2)*dot_n1_ka1 cross_n1_ka1(2) nx2(2)-kx2(2)*dot_n2_ka2 cross_n2_ka2(2); ...
         -ox1(2)+kx1(2)*dot_o1_ka1 cross_o1_ka1(2) ox2(2)-kx2(2)*dot_o2_ka2 cross_o2_ka2(2); ...
         -ax1(2)+kx1(2)*dot_a1_ka1 cross_a1_ka1(2) ax2(2)-kx2(2)*dot_a2_ka2 cross_a2_ka2(2); ...
         -nx1(3)+kx1(3)*dot_n1_ka1 cross_n1_ka1(3) nx2(3)-kx2(3)*dot_n2_ka2 cross_n2_ka2(3); ...
         -ox1(3)+kx1(3)*dot_o1_ka1 cross_o1_ka1(3) ox2(3)-kx2(3)*dot_o2_ka2 cross_o2_ka2(3); ...
         -ax1(3)+kx1(3)*dot_a1_ka1 cross_a1_ka1(3) ax2(3)-kx2(3)*dot_a2_ka2 cross_a2_ka2(3)];

    b = [-kx2(1)*dot_n2_ka2+kx1(1)*dot_n1_ka1; ...
         -kx2(1)*dot_o2_ka2+kx1(1)*dot_o1_ka1; ...
         -kx2(1)*dot_a2_ka2+kx1(1)*dot_a1_ka1; ...
         -kx2(2)*dot_n2_ka2+kx1(2)*dot_n1_ka1; ...
         -kx2(2)*dot_o2_ka2+kx1(2)*dot_o1_ka1; ...
         -kx2(2)*dot_a2_ka2+kx1(2)*dot_a1_ka1; ...
         -kx2(3)*dot_n2_ka2+kx1(3)*dot_n1_ka1; ...
         -kx2(3)*dot_o2_ka2+kx1(3)*dot_o1_ka1; ...
         -kx2(3)*dot_a2_ka2+kx1(3)*dot_a1_ka1];
     
    y = A\b;
    batahat1 = atan2(y(2),y(1));
    batahat2 = atan2(y(4),y(3));
    Rr1 = vec2rot(ka1.*batahat1)*R1;
    C = Rr1(1:3,1:3);
    Rr1(1:3,1:3) = C * inv(sqrtm(C'*C));
    if dim == 3        
        %% solve for R
        varargout{1} = Rr1;
    elseif dim == 4
        A1 = [Ra1 - eye(3);Ra2 - eye(3)];
        b1 = [Rr1*TB1(1:3,4)-TA1(1:3,4);Rr1*TB2(1:3,4)-TA2(1:3,4)];
        tr1 = A1\b1;
        varargout{1} = [Rr1 tr1;[0 0 0 1]];
    end
end