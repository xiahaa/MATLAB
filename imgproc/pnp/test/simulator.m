clc;
close all;

%% a data generator to test pose estimation for visual odometry
%% Author: xiahaa@space.dtu.dk

%% addpath if your code is in a nested directory
addpath ./3rdparty
addpath ./dataGen
addpath ./solver

%% experimental configurations
noise_level = 0.0:0.1:0.2;%0.5:0.5:5;
number_of_points = 50;
iterations = 100;
is_2d = false;

% camera's parameters
width= 640;
height= 480;
f= 800;
K = [f 0 width*0.5;0 f height*0.5;0 0 1];
outlier_percent = 0.3;

%% function handler to call you method
name= {'SVD'};
f= {@SVD33};
test_case = {@ordinary_case};

cmap = lines(size(f,2));
cmapalpha = [cmap];
color = cmapalpha;
markerfacecolor = cmapalpha;

B = zeros(iterations,1);%% dummy variable
A = zeros(size(test_case, 2), numel(noise_level));
method_list= struct('name', name, 'f', f, 'r', B, 't', B, ...
                    'mean_r', A, 'mean_t', A, ...
                    'med_r', A, 'med_t', A, 'std_r', A, 'std_t', A, ...
                    'marker', 'o', 'color', color, 'markerfacecolor', markerfacecolor);

totaltime = length(noise_level)*size(test_case, 2)*iterations;
currtime = 0;

f = waitbar(0,'Solving...','Name','Pose Estimation Simulation',...
    'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');

for i = 1:size(test_case, 2)
    gen = test_case{i};
    for j= 1:numel(noise_level)
        nl= noise_level(j);
        fprintf('noise_level = %.1f: ',nl);
        index_fail = [];
        for k = 1:iterations
            [P, Q, R, t] = gen(number_of_points);
            
            %% outlier
            num_outliers = round(outlier_percent * number_of_points);
            source = randperm(number_of_points, num_outliers);
            [Po, Qo, ~, ~] = gen(num_outliers);
            Q(:,source) = Qo;
            P(:,source) = Po;
            
            placeholder1 = [];
            placeholder2 = [];
            placeholder3 = [];
            if is_2d == true
                % projection
                q = K * (Q./Q(3,:));
                qn = q + randn(3,number_of_points).*nl;
                placeholder1 = P;
                placeholder2 = qn;
                placeholder3 = K;
            else
                P = P;
                Q = Q + randn(3,number_of_points).*nl;
                placeholder1 = P;
                placeholder2 = Q;
            end
            
            % pose estimation
            for kk = 1:length(method_list)
                [R1,t1]= method_list(kk).f(placeholder1, placeholder2, placeholder3);%% if you want to do normalization, do it internally.
                %In case of no solution
                if size(t1,2) ~= 1
                    index_fail = [index_fail, k];
                    break;
                end
                err = cal_pose_err([R1 t1],[R t]);
                method_list(kk).r(k)= err(1);
                method_list(kk).t(k)= err(2);
            end
            currtime = currtime + 1;
            waitbar(currtime / totaltime, f);
        end
        fprintf('\n');
        
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
    end
    
    % save result if you want
    % save('xxxx','xxxx');
end
delete(f);

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


