clc;
clear all;
close all;

%% simulation of homography decomposition
addpath('../../../MatrixLieGroup');
addpath('../../../quaternion');
addpath('../../../beautiful_plot');
addpath('../');
T1 = fakeRT();

N = 50;
p = rand([3,N]) * 5 - 2.5;
% p(3,:) = 5;%p(3,:);
% p(1,:) = p(1,:);
% p(2,:) = p(2,:);
s = 1;
q = s.*(T1(1:3,1:3)*p + repmat(T1(1:3,4),1,N));
q = q(:,1:N) + rand([3,N]).*0.0;

font_size = 14;
blue_color = [0 116 186]/255;
orange_color = [223 80 35]/255;
plot3(q(1,:),q(2,:),q(3,:),'o','MarkerSize',5,'MarkerEdgeColor',orange_color,'MarkerFaceColor',orange_color)
hold on;
plot3(p(1,:),p(2,:),p(3,:),'o','MarkerSize',5,'MarkerEdgeColor',blue_color,'MarkerFaceColor',blue_color)
xlabel('$x$','Interpreter','latex')
ylabel('$y$','Interpreter','latex')
zlabel('$z$','Interpreter','latex')

T1



% 
% p2(1,:) = rand([1,N]) * 10 - 2.5;
% p2(2,:) = rand([1,N]) * 7 - 2.5;
% p2(3,:) = rand([1,N]) * 3 - 2.5;
% q2 = T1(1:3,1:3)*p2 + repmat(T1(1:3,4),1,N);
% q2 = s.*q2(:,1:N);
% %% svd
% % [~,~,V1]=svd(p');
% % [~,~,V2]=svd(p2');
% % [~,~,V3]=svd(q');
% % [~,~,V4]=svd(q2');
% V1 = p - repmat(mean(p,2),1,size(p,2));
% V2 = p2 - repmat(mean(p2,2),1,size(p2,2));
% V3 = q - repmat(mean(q,2),1,size(q,2));
% V4 = q2 - repmat(mean(q2,2),1,size(q2,2));
% V1 = V1*V1'/(1-size(p,2));
% V2 = V2*V2'/(1-size(p2,2));
% V3 = V3*V3'/(1-size(q,2));
% V4 = V4*V4'/(1-size(q2,2));
% % [~,~,V1]=svd(p');
% % [~,~,V2]=svd(p2');
% % [~,~,V3]=svd(q');
% % [~,~,V4]=svd(q2');
% 
% SO3_p(:,:,1) = V1;SE3_p(1,:,:) = [V1 [0;0;0];[0 0 0 1]];
% SO3_p(:,:,2) = V2;SE3_p(2,:,:) = [V2 [0;0;0];[0 0 0 1]];
% SO3_q(:,:,1) = V3;SE3_q(1,:,:) = [V3 [0;0;0];[0 0 0 1]];
% SO3_q(:,:,2) = V4;SE3_q(2,:,:) = [V4 [0;0;0];[0 0 0 1]];
% 
% A = zeros(9*2,9);
% for i = 1:2
%     Ra = SO3_q(1:3,1:3,i);
%     Rb = SO3_p(1:3,1:3,i);
% 
%     A(9*i-8:9*i,1:9) = eye(9) - kron(Ra,Rb);
% end
% [~,Sx,Vx] = svd(A);
% X = reshape(Vx(1:9),3,3)';
% X = sign(det(X))/abs(det(X))^(1/3)*X;
% [u, ~, v] = svd(X); X = u*v'; if det(X)<0, X = u*diag([1 1 -1])*v'; end
% X

    addpath C:\Users\xiahaa\Documents\MATLAB\hand_eye_calib\solver\
    addpath C:\Users\xiahaa\Documents\MATLAB\hand_eye_calib\stochastics
    addpath C:\Users\xiahaa\Documents\MATLAB\hand_eye_calib\stochastics\mean
    addpath C:\Users\xiahaa\Documents\MATLAB\hand_eye_calib\stochastics\utils
    addpath C:\Users\xiahaa\Documents\MATLAB\MatrixLieGroup\barfoot_tro14\
%     [a1,a2,a3]  = size(SE3_q);
%     A_noise_mex = reshape(SE3_q, a1, a2*a3);
%     [a1,a2,a3]  = size(SE3_p);
%     B_mex = reshape(SE3_p, a1, a2*a3);
tic
    opt = 2;
%     [X_solved, MA, MB, SigA, SigB] = batchSolveNew(SE3_q, SE3_p, opt);
%     dT = sol_andreff(SE3_p, SE3_q, size(SE3_q,1));

    nq = size(q,2);
    np = size(p,2);
    
    %% step 1: forme SO(3) for qn and P
    set_q = nchoosek(1:nq,2);
    set_p = nchoosek(1:np,2);
    
    SO3_q = formeSO3(q, set_q);
    SO3_p = formeSO3(p, set_p);
    
    
    %%% try directly work with lie algebra, doesn't work %%%
%     b1 = zeros(3,1);
%     b2 = zeros(3,1);
%     A = zeros(3,3);
%     for i = 1:size(SO3_q,3)
%         so31 = rot2vec(SO3_p(:,:,i));
%         so32 = rot2vec(SO3_q(:,:,i));
%         A = 0.5.*skewm(so31);
%         b1 = b1 + so31;
%         b2 = b2 + so32;
%     end
%     b1 = b1./N;
%     b2 = b2./N;
%     A = A./N;
%     A = eye(3) - A;
%     b = b2 - b1;
%     x = A\b;
%     R = vec2rot(x)
    %%% try directly work with lie algebra, doesn't work %%%
    
    
%     SO3_q = rankingSO3(SO3_q);
%     SO3_p = rankingSO3(SO3_p);
    
    %% step 2: forme relative SO(3) for qn and P
%     raw_SO3_q = forme_relative_SO3(SO3_q);
%     raw_SO3_p = forme_relative_SO3(SO3_p);


    %% planar motion
%     ph = [p;ones(1,np)];
%     qh = [q;ones(1,nq)];
%     
%     A1 = ph';
%     A2 = qh';
%     [~,~,V1] = svd(A1);
%     [~,~,V2] = svd(A2);
%     plane1 = V1(:,end);
%     plane2 = V2(:,end);
%     
%     n1 = plane1(1:3);n1 = n1./norm(n1);
%     n2 = plane2(1:3);n2 = n2./norm(n2);
%     
%     % ok
%     n3 = cross(n1,n2);
%     R1 = eye(3) + skewm(n3) + (1-n1'*n2)/(n3'*n3)*(skewm(n3))^2;
    R1 = eye(3);

    SO3_p1 = SO3_p;
%     for i = 1:size(SO3_p,3)
%         SO31 = SO3_p(:,:,i);
% %         SO32 = SO3_q(:,:,i);
%         SO3_p1(:,:,i) = R1 * SO31;
%     end
    n_search = int16(2*10^2);
    j = 0;
%     threshold = 5/57.3;
%     old_inlier_p = -1;
    old_var = -1;
    while j < 100
        Mq = mean_Taylor_1st( SO3_q );
        SigA = zeros(3,3);
        for i = 1:size(SO3_p1,3)
            X_i = SO3_p1(:,:,i);
            P = Mp1\X_i;
            SigA = SigA + so3_vec(logm(P))*so3_vec(logm(P))';
        end
        stdA = diag(sqrtm(SigA));
        for i = 1:size(SO3_p1,3)
            
        end
        
        Mp1 = mean_Taylor_1st( SO3_p1 );
        SigB = zeros(3,3);
        for i = 1:size(SO3_q,3)
            X_i = SO3_q(:,:,i);
            P = Mq\X_i;
            SigB = SigB + so3_vec(logm(P))*so3_vec(logm(P))';
        end
        
        p1 = 1/((2*pi)^3*det(SigA));
        p2 = 1/((2*pi)^3*det(SigB));
        
        if p1 > 0.8 && p2 > 0.8
            break;
        end
        
        
        
        j = j + 1;
    end
        R2 = Mq*inv(Mp1);
        R = R2*R1;
        
        
        
        
%         err = [];
%         for i = 1:size(SO3_p1,3)
%             err(i) = norm(logm(SO3_q(:,:,i)'*R*SO3_p1(:,:,i)),'fro');
%         end
%         new_var = var(err);

%         std = sqrt(var(err));
%         inliers = err < threshold;
%         SO3_p1 = SO3_p1(:,:,inliers);
%         SO3_q = SO3_q(:,:,inliers);
%         new_inlier_p = numel(find(inliers))/numel(err);
%         if (new_inlier_p - old_inlier_p) < 0.01 || size(SO3_p1,3) < 5
%             break;
%         end
%         old_inlier_p = new_inlier_p;

    R = R2*R1
toc
    norm(logm(R'*T1(1:3,1:3)),'fro')*57.3

%     SigA = zeros(3,3);
%     for i = 1:size(SO3_p1,3)
%         X_i = SO3_p1(:,:,i);
%         P = Mp1\X_i;
%         SigA = SigA + so3_vec(logm(P))*so3_vec(logm(P))';
%     end
%     SigA = zeros(3,3);
%     for i = 1:size(SO3_p1,3)
%         X_i = SO3_p1(:,:,i);
%         P = Mp1\X_i;
%         SigA = SigA + so3_vec(logm(P))*so3_vec(logm(P))';
%     end
    
    
function T = fakeRT()
    euler(1) = (rand(1)*pi/2 - pi/4);
    euler(2) = (rand(1)*pi/2 - pi/4);
    euler(3) = (rand(1)*2*pi - pi);
    R1 = euler2rot(euler(1),euler(2),euler(3));
    t1 = rand([1,3]) * 2 + 3;
    t1 = t1';
    t1(1) = t1(1);
    t1(2) = t1(2);
    t1(3) = 1;%t1(3);
    T = [R1 t1;[0 0 0 1]];
end

function goodSO3s = rankingSO3(SO3s)
    angles = zeros(1,size(SO3s,3));
    for i = 1:size(SO3s,3)
        angles(i) = abs(norm(rot2vec(SO3s(:,:,i))));
    end
    [angless,id] = sort(angles);
    num = 100;
    if size(SO3s,3) < 100
        num = size(SO3s,3);
    end
    goodSO3s = SO3s(:,:,id(1:num));
end

function SO3s = formeSO3(p, setp)
    n_so3 = size(setp,1);
    SO3s = zeros(3,3,n_so3);
    center = mean(p,2);
    for i=1:n_so3
        ii = setp(i,1);
        jj = setp(i,2);
%         kk = setp(i,3);
        v1 = center-p(:,ii); 
        v2 = center-p(:,jj);
        v1 = v1./norm(v1);v2 = v2./norm(v2);
        v3 = cross(v1,v2);v3 = v3./norm(v3);
        v2 = cross(v3,v1);v2 = v2./norm(v2);
        SO3s(:,:,i) = [v1 v2 v3];
    end
end

function relative_SO3s = forme_relative_SO3(SO3s)
    nSO3s = size(SO3s,3);
    set_r_SO3s = nchoosek(1:nSO3s, 2);
    n_r_SO3s = size(set_r_SO3s,1);
    relative_SO3s = zeros(3,3,n_r_SO3s);
    for i = 1:1:n_r_SO3s
        ii = set_r_SO3s(i,1);
        jj = set_r_SO3s(i,2);
        relative_SO3s(:,:,i) = SO3s(:,:,ii)*SO3s(:,:,jj)';
    end
end

function T = orthog(M_hat)
% This function takes in a 4 by 4 matrix which has the form of
% [3 by 4; 0 0 0 n] and return a rigig transformation matrix by
% orthogonalizing the upper left 3 by 3 matrix into a rotation matrix
    opt = 2;   % opt = 1: Orthogonalize only the rotation part of M_hat 
               % opt = 2: SVD the whole M_hat to obtain a transformation matrix
    if opt == 1

        R_hat = M_hat(1:3,1:3); % Note that R_hat doesn't belong to SO(3)

        % Project onto SO(3)
        Rx = R_hat*(R_hat'*R_hat)^(-1/2);
        Rx = sign(det(Rx))/abs(det(Rx))^(1/3)*Rx;

        T = Rx;

    elseif opt == 2

        [U, ~, V] = svd(M_hat(1:3,1:3));

        % Note that the order of U and V matters
        Rx = U*V';

        if det(Rx) < 0
            Rx(:,3) = - Rx(:,3);  
        end

        if abs( det(Rx) - 1) > 10^-8
            fprintf('A non rotation matrix is returned');
        end

        T = Rx;
    end
end

function [M_1, M_hat] = mean_Taylor_1st( X ) 
    % This function calculates the 1st order approximation of the mean of a
    % bunch of matrices based on the Taylor expansion of the matrix logarithm
    % and the definition of mean of a probability density function.

    % Input: X is a 4 by 4*n rigid transformation matrices
    % Output: M_T1 is the mean of the 1st order approximation of Taylor
    % expansion

    % Change of this m file doesn't automatically change the executable generated by 
    % mean_Taylor_2nd.m
    n =  size(X, 3);

    M_hat = zeros(3);
    M_1 = zeros(3);

    for i = 1:n
        M_hat = M_hat + X(1:3,1:3,i);
    end

    M_hat = 1/n*M_hat;  % Note that M_hat doesn't belong to SE(3)

    M_1 = orthog(M_hat);
end

function MX = mean_Taylor_2nd_SO3( X, num )
    %% Calculate the 2nd order approximation of the mean of a bunch of matrices 
    % based on the Taylor expansion of the matrix logarithm and the definition 
    % of mean of a probability density function.
    %
    % Input: 
    %       X : a cell of 4 x 4 x n SE3 matrices
    % Output: 
    %       MX : the 2nd order approximation of the Taylor expansion
    %
    % Note : M_T1 is the mean of the 1st order approximation of Taylor
    % expansion
    %
    %
    % coder.extrinsic('mean_Taylor_1st_mex');
    %%
    M_t1 = mean_Taylor_1st( X );
    n = size(X,3);
    % Perturb M_t1 along each direction defined by the Lie algebra element of
    % SE(3)
    CA = inf;
    diff = 10^-12;
    eps  =  0.001;
    
    E = zeros(3,3,3);
    E(:,:,1) = [0  0 0; 0 0 -1;  0 1 0];
    E(:,:,2) = [0  0 1; 0 0  0; -1 0 0];
    E(:,:,3) = [0 -1 0; 1 0  0;  0 0 0];

    MX = M_t1;
    count = 0;
 
    step = 1; % Change the step lengths if descent direction is not found

    while( abs(CA) > diff )
        count = count + 1;
        n_step = 0;
        % Approximation of rotation is good so far so we only further optimize
        % on translation
        for j = 1:3
            MX1 = MX*expm( step*eps*E(:,:,j));
            MX2 = MX*expm(-step*eps*E(:,:,j));

            MX1sum = zeros(3);
            MX2sum = zeros(3);

            for k = 1:n
                X_k = X(:,:,k);
                MX1sum = MX1sum + X_k/MX1*X_k;
                MX2sum = MX2sum + X_k/MX2*X_k;
            end

            %% equation 39 in his paper
            MX1cost = 2*M_t1 - 1/2/n*MX1sum - 3/2*MX1;
            MX2cost = 2*M_t1 - 1/2/n*MX2sum - 3/2*MX2;

            CA1 = norm(MX1cost)^2;
            CA2 = norm(MX2cost)^2;

            if CA1 <= CA2 && CA1 < CA
                CA = CA1;
                MX = MX1;
                step = 1;
                % fprintf('Found a descent direction along %d \n', j)
                % disp(MX1cost)
                % disp(CA)
            elseif CA2 < CA1 && CA2 < CA
                CA = CA2;
                MX = MX2;
                step = 1;
                % fprintf('Found a descent direction along -%d \n', j)
                % disp(MX2cost)
                % disp(CA)
            else
                % fprintf('Not a descent direction along +/-%d \n', j)
                n_step = n_step + 1;
                if n_step == 3
                   step = step + 1;
                end
            end
        end
        if count == num
            break;
        end
    end
    % fprintf('\n')
    fprintf('Number of iterations is %f', count)
    % fprintf('\n')
    % fprintf('\n')
end


