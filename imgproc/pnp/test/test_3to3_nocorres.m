clc;
clear all;
close all;

%% simulation of homography decomposition
addpath('../../../MatrixLieGroup');
addpath('../../../quaternion');
addpath('../../../beautiful_plot');
addpath('../');
T1 = fakeRT();
tform = affine3d(T1');
% ptCloud = pcread('teapot.ply');
% ptCloud = pcread('C:\Users\xiahaa\3rdparty\opencv-3.1.0\samples\cpp\tutorial_code\viz\bunny.ply');
% gridStep = 0.02;
% ptCloud = pcdownsample(ptCloud,'gridAverage',gridStep);
% ptCloudTformed = pctransform(ptCloud,tform);
% figure
% pcshow(ptCloud); hold on;
% pcshow(ptCloudTformed); 
% tic
% tform1 = pcregistericp(ptCloudTformed,ptCloud,'Extrapolate',true);
% toc
% tform2 = invert(tform1);
% disp(tform2.T);
% p = ptCloud.Location;
% q = ptCloudTformed.Location;
% p = p';
% q = q';

N = 50;
p = rand([3,N]) * 5 - 2.5;
% p(3,:) = 5;%p(3,:);
p(1,:) = p(1,:);
p(2,:) = p(2,:);
s = 1;
q = s.*(T1(1:3,1:3)*p + repmat(T1(1:3,4),1,N));
q = q(:,1:N) + rand([3,N]).*0.1;

font_size = 14;
blue_color = [0 116 186]/255;
orange_color = [223 80 35]/255;
plot3(q(1,:),q(2,:),q(3,:),'o','MarkerSize',5,'MarkerEdgeColor',orange_color,'MarkerFaceColor',orange_color)
hold on;
plot3(p(1,:),p(2,:),p(3,:),'o','MarkerSize',5,'MarkerEdgeColor',blue_color,'MarkerFaceColor',blue_color)
xlabel('$x$','Interpreter','latex')
ylabel('$y$','Interpreter','latex')
zlabel('$z$','Interpreter','latex')

disp(T1)

tic    
    nq = size(q,2);
    np = size(p,2);
    %% step 1: forme SO(3) for qn and P
    set_q = nchoosek(1:nq,2);
    set_p = nchoosek(1:np,2);
    
    SO3_q = formeSO3(q, set_q);
    SO3_p = formeSO3(p, set_p);
   
%     SO3_q = rankingSO3(SO3_q);
%     SO3_p = rankingSO3(SO3_p);
    
    R1 = eye(3);

    SO3_p1 = SO3_p;
    
    Mq1 = mean_1st_order(SO3_q);
    Mq2 = mean_iterative(SO3_q,Mq1);
    Mq3 = FNS_iterative(SO3_q,Mq1);

    Mp1 = mean_1st_order(SO3_p1);
    Mp2 = mean_iterative(SO3_p1,Mp1);
    Mp3 = FNS_iterative(SO3_p1,Mp1);
    
    R3 = Mq3*inv(Mp3)
    R2 = Mq2*inv(Mp2)
    R = Mq1*inv(Mp1)
    
toc
%     norm(logm(R2'*T1(1:3,1:3)),'fro')*57.3

    
function M = mean_1st_order(SO3s)
    addpath C:\Users\xiahaa\Documents\MATLAB\MatrixLieGroup\barfoot_tro14\
    M = zeros(3,3);
    for i = 1:size(SO3s,3)
        M = M + SO3s(:,:,i);
    end
    M = M ./ size(SO3s,3);
    %% orthogonalization
    [U, ~, V] = svd(M); 
    M = U*V';
    if det(M)<0, M = U*diag([1 1 -1])*V'; end 
end

function M = mean_iterative(SO3s,M_1st)
    
    M = M_1st;
    iter = 1;
    maxIter = 50;
    err = 0;
    old_err = 1e6;
%     
%     while iter < maxIter
%         e = zeros(3,1);
%         err = 0;
%         for i = 1:size(SO3s,3)
%             Rbar = M'*SO3s(:,:,i);
%             phibar = rot2vec(Rbar);
%             err = err + norm(phibar);
%             e = e + M*phibar;
%         end
%         disp([old_err, err])
%         if abs(err - old_err) < 1e-6 
%             break;
%         end
%         xhat = e./size(SO3s,3);
%         M = vec2rot(xhat)*M;
%         old_err = err;
%         iter = iter + 1;
%         
%     end
%     disp(iter)
%     

%     Solve for pose using our algorithm
    Asum = M .* -2;
    for i=1:maxIter      % Gauss-Newton iterations
        J = zeros(9,9);
        b = zeros(3,3);
        for k=1:size(SO3s,3)
            a1 = kron(SO3s(:,:,k)', SO3s(:,:,k)*M');
            J = J + a1;
            b = b + SO3s(:,:,k)*M'*SO3s(:,:,k);
        end
        b = b ./ size(SO3s,3) ./ 2;
        bb = vec(Asum+b+1.5.*M);
        J = J./size(SO3s,3) ./ 2 - 1.5.*kron(M',eye(3));
        o = J\bb;
        Rx = [o(1:3) o(4:6) o(7:9)];Rx = (Rx-Rx').*0.5;
        
        M = expm(Rx)*M;
        
        if norm(bb) < 1e-20
            break;
        end
    end

%     Asum = M .* 2;
%     for i=1:maxIter      % Gauss-Newton iterations
%         A1 = zeros(3,3);
%         A2 = zeros(9,9);
%         for k=1:size(SO3s,3)
%             AMA = SO3s(:,:,k)*M'*SO3s(:,:,k);
%             A1 = A1 + AMA;
%             A2 = A2 + kron(AMA',SO3s(:,:,k));
%         end
%         A3 = A1./(size(SO3s,3)*2);
%         b = vec(Asum*2 - A3 - 1.5.*M);
%         J = kron(M',1.5.*eye(3)) - A2./(size(SO3s,3)*2);
%         o = J\b;
%         Rx = [o(1:3) o(4:6) o(7:9)];Rx = (Rx-Rx').*0.5;
%         M = M*expm(Rx);
%         if norm(b) < 1e-20
%             break;
%         end
%     end
    
    %% orthogonalization
    [U, ~, V] = svd(M); 
    M = U*V';
    if det(M)<0, M = U*diag([1 1 -1])*V'; end 
end

function M = FNS_iterative(SO3s,M_1st)
    M = M_1st;
    iter = 1;
    maxIter = 50;
    
    while iter < maxIter
        Minv = (M)';
        LHS = zeros(3,3);
        RHS = zeros(3,1);
        for k=1:size(SO3s,3)
            v1 = rot2vec(Minv*SO3s(:,:,k));
            RHS = RHS + v1;
            LHS = LHS + vec2jacInv(v1);
        end
        xhat = LHS\RHS;
        M = M*vec2rot(xhat);
        if norm(xhat) < 1e-3
            break;
        end
        iter = iter + 1;
    end
end
    
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


