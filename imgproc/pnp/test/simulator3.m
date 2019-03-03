clc;
close all;

%% a data generator to test pose estimation for visual odometry
%% Author: xiahaa@space.dtu.dk

%% addpath if your code is in a nested directory
addpath ./3rdparty
addpath ./dataGen
addpath ./solver
addpath(genpath('../../../MatrixLieGroup'));
addpath('../../../quaternion');
addpath('../../../beautiful_plot');
addpath('../');
addpath('./eurasip/');

% prefix = 'C:/Users/xiahaa/Documents/MATLAB/imgproc/pnp/test/eurasip/data/';
% name = {'ordinary'};% 'plane'   % quasi-singular
% suffix = '.mat';
% %     stds = [0 0.1 0.25 0.5 0.75 1];
% stds = [0.0 0.01 0.02 0.03 0.04 0.05];

%% experimental configurations
noise_level = 0:0.001:0.001;%:0.02:0.1;%:0.1:0.2;%0.5:0.5:5;
number_of_points = 20;%:20:50;
iterations = 20;

outlier_percent = 0.1;
random_order_percent = 0.0;

%% function handler to call you method
name= {'SVD'}%,'PROPOSED','ICP'};
f= {@pose_estimation_by_voting};%@SVD33,@pose_estimation_correspondence_free,@pcregistericp};

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
test_case = {@ordinary_case};

totaltime = length(noise_level)*length(number_of_points)*iterations;
currtime = 0;

usepcl = [0 0 1];

% f = waitbar(0,'Solving...','Name','Pose Estimation Simulation',...
%     'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');

    gen = @datagen;
    for i = 1:length(number_of_points)
        num = number_of_points(i);
        for j= 1:numel(noise_level)
            nl= noise_level(j);
            fprintf('noise_level = %.3f: ',nl);
            index_fail = [];
            for k = 1:iterations
                [P, Q, R, t] = gen(num);

%            P = [-1.9496   -1.4197    2.1843   -0.7022    1.6351; ...
%                 -1.1039   -2.3297   -1.1896   -2.3658   -1.2051; ...
%                  1.3382   -0.3172    0.3487    0.0021   -2.2706];
% 
% 
%             Q = [ ...
%                 3.9088    5.2289    7.1965    5.7747    6.9922; ...
%                 1.9878    1.8314    4.6834    2.1515    5.0434; ...
%                 5.3174    3.6570    5.3085    4.1236    2.6644;];
                
                %% outlier
                num_outliers = round(outlier_percent * num);
                source = randperm(num, num_outliers);
%                 [Po, Qo, ~, ~] = gen(num_outliers);
                Po = rand(3,num_outliers);
                Qo = rand(3,num_outliers);
                Q(:,source) = Qo;
%                 P(:,source) = Po;
                R
                num_scrambles = round(random_order_percent * num);
                source = randperm(num, num_scrambles);
%                 source_s = source + 1;
%                 source_s(source_s>num) = 1;
%                 tmp = P(:,source_s);
%                 P(:,source_s) = P(:,source);
%                 P(:,source) = tmp;
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

                % pose estimation
                for kk = 1:length(method_list)
                    if usepcl(kk) == 1
                        ptCloud1 = pointCloud(P');  
                        ptCloud2 = pointCloud(Q');  
                        tform = pcregistericp(ptCloud2,ptCloud1,'Extrapolate',true);
                        tform2 = invert(tform);
                        R1 = tform2.T(1:3,1:3);
                        t1 = tform2.T(4,1:3);
                    else
                        [R1,t1,R2,t2]= method_list(kk).f(placeholder1, placeholder2, placeholder3);%% if you want to do normalization, do it internally.
                    end
                    %In case of no solution
                    if size(t1,2) ~= 1
                        index_fail = [index_fail, k];
                        break;
                    end
                    err = cal_pose_err([R1 t1],[R t]);
                    err1 = cal_pose_err([R2 t2],[R t]);
                    method_list(kk).r(k)= err(1);
                    method_list(kk).t(k)= err(2);
                end
                currtime = currtime + 1;
%                 waitbar(currtime / totaltime, f);
            end
            for k = 1:length(method_list)
                method_list(k).r(index_fail) = [];
                method_list(k).t(index_fail) = [];

                method_list(k).mean_r(i,j)= mean(method_list(k).r);
                method_list(k).mean_t(i,j)= mean(method_list(k).t);
                method_list(k).med_r(i,j)= median(method_list(k).r);
                method_list(k).med_t(i,j)= median(method_list(k).t);
                method_list(k).std_r(i,j)= std(method_list(k).r);
                method_list(k).std_t(i,j)= std(method_list(k).t);
            end
            
            fprintf('\n');
        end
    end
    % save result if you want
    % save('xxxx','xxxx');

%     delete(f);

%% do some plot if you want
yrange= [0 20];w= 350; h= 350;
for i = 1:size(test_case, 2)
    figure('color','w','position',[w*i,100,w,h]);
    figtitle = strcat('Mean Rotation Error: ', func2str(test_case{i}));
    xdrawgraph(noise_level,yrange,method_list,'mean_r', figtitle, ...
        'Gaussian Noise std','Rotation Error (degrees)');
    figure('color','w','position',[w*i,100,w,h]);
    figtitle = strcat('Median Rotation Error: ', func2str(test_case{i}));
    xdrawgraph(noise_level,yrange,method_list,'med_r', figtitle, ...
        'Gaussian Noise std','Rotation Error (degrees)');
    figure('color','w','position',[w*i,100,w,h]);
    figtitle = strcat('Mean Translation Error: ', func2str(test_case{i}));
    xdrawgraph(noise_level,yrange,method_list,'mean_t', figtitle, ...
        'Gaussian Noise std','Rotation Error (degrees)');
    figure('color','w','position',[w*i,100,w,h]);
    figtitle = strcat('Median Translation Error: ', func2str(test_case{i}));
    xdrawgraph(noise_level,yrange,method_list,'med_t', figtitle, ...
        'Gaussian Noise std','Rotation Error (degrees)');
end

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
    euler(1) = (rand(1)*pi/2 - pi/4);
    euler(2) = (rand(1)*pi/2 - pi/4);
    euler(3) = (rand(1)*2*pi - pi);
    R1 = euler2rot(euler(1),euler(2),euler(3));
    t1 = rand([1,3]) * 2 + 3;
    t1 = t1';
    t1(1) = t1(1);
    t1(2) = t1(2);
    t1(3) = t1(3);
    T = [R1 t1;[0 0 0 1]];
end