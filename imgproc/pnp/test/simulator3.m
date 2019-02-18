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

%% experimental configurations
noise_level = 0.0:0.1:0.2;%0.5:0.5:5;
number_of_points = 10:20:50;
iterations = 20;

outlier_percent = 0.0;
random_order_percent = 0.3;

%% function handler to call you method
name= {'SVD'}%,'PROPOSED','ICP'};
f= {@pose_estimation_correspondence_free};%@SVD33,@pose_estimation_correspondence_free,@pcregistericp};

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

usepcl = [0 0 1];

% f = waitbar(0,'Solving...','Name','Pose Estimation Simulation',...
%     'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');

%     gen = test_case{i};
    for i = 1:length(number_of_points)
        num = number_of_points(i);
        for j= 1:numel(noise_level)
            nl= noise_level(j);
            fprintf('noise_level = %.1f: ',nl);
            index_fail = [];
            for k = 1:iterations
                [P, Q, R, t] = gen(num);

                %% outlier
                num_outliers = round(outlier_percent * num);
                source = randperm(num, num_outliers);
                [Po, Qo, ~, ~] = gen(num_outliers);
                Q(:,source) = Qo;
                P(:,source) = Po;

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
                        [R1,t1, R2, t2]= method_list(kk).f(placeholder1, placeholder2, placeholder3);%% if you want to do normalization, do it internally.
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


