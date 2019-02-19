clc;
close all;

%% a data generator to test pose estimation for visual odometry
%% Author: xiahaa@space.dtu.dk

%% addpath if your code is in a nested directory
addpath ./3rdparty
addpath ./dataGen
addpath ./solver
addpath('../../../MatrixLieGroup');
addpath('../../../quaternion');
addpath('../../../beautiful_plot');
addpath('../');
addpath('./eurasip/');

prefix = 'C:/Users/xiahaa/Documents/MATLAB/imgproc/pnp/test/eurasip/data/';
% name = {'ordinary'};% 'plane'   % quasi-singular
suffix = '.mat';

%% experimental configurations
noise_level = 0.01:0.01:0.05;
number_of_points = 100;%:10:100;%:20:50;
iterations = 50;

random_order_percent = 0.5;

%% function handler to call you method
name= {'PROPOSED','ICP'};
f= {@pose_estimation_correspondence_free,@pcregistericp};%@SVD33,@pose_estimation_correspondence_free,@pcregistericp};

cmap = lines(size(f,2));
cmapalpha = [cmap];
color = cmapalpha;
markerfacecolor = cmapalpha;

B = zeros(iterations,1);%% dummy variable
A = zeros(length(number_of_points), numel(noise_level));
method_list= struct('name', name, 'f', f, 'r', B, 't', B, ...
                    'mean_r', A, 'mean_t', A, ...
                    'med_r', A, 'med_t', A, 'std_r', A, 'std_t', A, ...
                    'marker', 'o', 'color', color, 'markerfacecolor', markerfacecolor);

totaltime = length(noise_level)*length(number_of_points)*iterations;
currtime = 0;

usepcl = [0 1];

% f = waitbar(0,'Solving...','Name','Pose Estimation Simulation',...
%     'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');
    gen = @datagen;
    for i = 1:length(number_of_points)
        num = number_of_points(i);
        for j= 1:numel(noise_level)
            nl= noise_level(j);
            fprintf('noise_level = %.3f: ',nl);
            
            Rs = zeros(3,3,iterations);
            ts = zeros(3,1,iterations);
            
            for k = 1:iterations
                disp(k)
                [P, Q, R, t] = gen(num);
                
                Rs(:,:,k) = R;
                ts(:,:,k) = t;
                
                num_scrambles = round(random_order_percent * num);
                source = randperm(num, num_scrambles);
                Ps = P;
                if numel(source)>0
                    raw = P(:,source(1));
                    Ps(:,source(1:end-1)) = P(:,source(2:end));
                    Ps(:,source(end)) = raw;
                end

                placeholder1 = [];
                placeholder2 = [];
                placeholder3 = [];

%                 P = P;
                Q = Q + randn(3,num).*nl;
                placeholder1 = Ps;
                placeholder2 = Q;
                
                fid = fopen(strcat(prefix,'data_',num2str(num),'_',num2str(nl), '_', num2str(k), 'txt'), 'w');
                fprintf(fid,'%6.10f  %6.10f %6.10f\n',P');
                fclose(fid);
                fid = fopen(strcat(prefix,'model_',num2str(num),'_',num2str(nl), '_', num2str(k), 'txt'), 'w');
                fprintf(fid,'%6.10f  %6.10f %6.10f\n',Q');
                fclose(fid);

%                 fontsize = 15;
%                 cmap = lines(2);
%                 cmapalpha = [cmap 0.8*ones(size(cmap,1),1)];
%                 figure;plot3(P(1,:), P(2,:), P(3,:), 'o', 'MarkerSize', 15, 'Color',cmapalpha(1,:));hold on;grid on;
%                 plot3(Q(1,:), Q(2,:), Q(3,:), '*', 'MarkerSize', 15, 'Color',cmapalpha(2,:));
% %                 title('Point Cloud','Interpreter','latex','FontSize',fontsize);
%                 xlabel('$x: (m)$','Interpreter','latex','FontSize',fontsize);
%                 ylabel('$y: (m)$','Interpreter','latex','FontSize',fontsize);
%                 zlabel('$z: (m)$','Interpreter','latex','FontSize',fontsize);
                % pose estimation
                for kk = 1:length(method_list)
                    if usepcl(kk) == 1
                        ptCloud1 = pointCloud(P');  
                        ptCloud2 = pointCloud(Q');  
                        tform = pcregistericp(ptCloud2,ptCloud1,'Extrapolate',true);
                        tform2 = invert(tform);
                        R1 = tform2.T(1:3,1:3);
                        t1 = tform2.T(4,1:3)';
                    else
                        [R1,t1, R2, t2]= method_list(kk).f(placeholder1, placeholder2, placeholder3);%% if you want to do normalization, do it internally.
                    end
         
                    err = cal_pose_err([R1 t1],[R t]);
%                     err1 = cal_pose_err([R2 t2],[R t]);
                    method_list(kk).r(k)= err(1);
                    method_list(kk).t(k)= err(2);
                    
%                     PQ = R1'*Q - R1'*repmat(t1,1,size(Q,2));
%                     
%                     figure;plot3(P(1,:), P(2,:), P(3,:), 'o', 'MarkerSize', 15, 'Color',cmapalpha(1,:));hold on;grid on;
%                     plot3(PQ(1,:), PQ(2,:), PQ(3,:), '*', 'MarkerSize', 15, 'Color',cmapalpha(2,:));
%     %                 title('Point Cloud','Interpreter','latex','FontSize',fontsize);
%                     xlabel('$x: (m)$','Interpreter','latex','FontSize',fontsize);
%                     ylabel('$y: (m)$','Interpreter','latex','FontSize',fontsize);
%                     zlabel('$z: (m)$','Interpreter','latex','FontSize',fontsize);
                    
                end
%                 currtime = currtime + 1;
%                 waitbar(currtime / totaltime, f);
            end
            save(strcat(prefix,'gt_',num2str(num),'_',num2str(nl),suffix),'Rs','ts');
            save(strcat(prefix,'noise_cmp_',num2str(num),'_',num2str(nl),suffix),'method_list');
            fprintf('\n');
        end
    end
    % save result if you want
    
%% do some plot if you want




function [p,q,R,t] = datagen(N)
    T1 = fakeRT();
    p = rand([3,N]) * 5 - 2.5;
    % p(3,:) = 5;%p(3,:);
    p(1,:) = p(1,:);
    p(2,:) = p(2,:);
    s = 1;
    q = s.*(T1(1:3,1:3)*p + repmat(T1(1:3,4),1,N));
    q = q(:,1:N) + rand([3,N]).*0.0;
    R = T1(1:3,1:3);
    t = T1(1:3,4);
end

function T = fakeRT()
    euler(1) = (rand(1)*pi - pi/2);
    euler(2) = (rand(1)*pi - pi/2);
    euler(3) = (rand(1)*2*pi - pi);
    R1 = euler2rot(euler(1),euler(2),euler(3));
    t1 = rand([1,3]) * 2 + 3;
    t1 = t1';
    t1(1) = t1(1);
    t1(2) = t1(2);
    t1(3) = t1(3);
    T = [R1 t1;[0 0 0 1]];
end