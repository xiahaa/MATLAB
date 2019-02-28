% experiment on different size, no outlier, varing noisy level and scrambling
% rate.
clear; 
clc; 
close all;

if ismac
    addpath('/Users/xiaohu/Documents/MATLAB/3rdparty/3dcalib/3dslam');
    addpath('/Users/xiaohu/Documents/MATLAB/3rdparty/3dcalib/3d-dualquat');
    addpath('/Users/xiaohu/Documents/MATLAB/3rdparty/3dcalib/3d-euler');
else
    addpath('../3rdparty/mit3dslam');
end
addpath('../beautiful_plot');

%% Add file dependencies
addpath C:\Users\xiahaa\Documents\MATLAB\hand_eye_calib\solver\
addpath C:\Users\xiahaa\Documents\MATLAB\MatrixLieGroup\barfoot_tro14\
addpath C:\Users\xiahaa\Documents\MATLAB\hand_eye_calib\data_gen\
addpath D:\dtu\sourcecode\hand_eye\axxb_calibration-stable_stochastics_lie\axxb_calibration-stable\matlab\Batch_Method_ED_KL_BS

%% Initialize Parameters
num = 50; % Number of steps
gmean = [0;0;0;0;0;0];	%Gaussian Noise Mean
stds = 1; % Gaussian Noise standard deviation Range
n_trials = 50; %60
perm_rate = 50;

% percentage_of_added_sample = 0:0.1:0.5;
outlier_percentage = 0:0.1:0.5;

convSolver = {
    @batchSolveNew, ...                                   %% batch1
    @batchSolveNew, ...                                   %% batch2
    @batchSolveSoftUseScrew, ...                          %% batch1
};

solver_name = {'B1','B2','BS'};

usedsolver = convSolver;

options = {1,4,[]};

%% Computation Loops
optPDF = 1; % select the distribution for A and B
A_noise = [];
B = [];

for solver_id = 1:size(solver_name,2)
    times = zeros(numel(outlier_percentage),100);
    error_r = zeros(numel(outlier_percentage),n_trials);
    error_t = zeros(numel(outlier_percentage),n_trials);
    valid_id = ones(numel(outlier_percentage),n_trials);
    for id1 = 1:numel(outlier_percentage)
        numout = outlier_percentage(id1)*num;
        disp([solver_id,id1])
        std = stds(1);
        cov = std*eye(6,6);
        for k = 1:n_trials
            x = randn(6,1); x = x./norm(x); X = expm(se3_vec(x)); % Generate a Random X
            while 1
                xout = randn(6,1); 
                xout = xout./norm(xout);
                if norm(xout-x) > 0.8
                    break;
                end
            end
            Xout = expm(se3_vec(xout)); % Generate a Random Xout
            
            [A, B] = generateAB(num, optPDF, X, gmean, cov);
%             [Aout, Bout] = generateAB(numout, optPDF, X, gmean, cov);
            Aout = genRandomPose(numout);
            Bout = genRandomPose(numout);
        
%             if rand(1) > 0.5
%                 A = cat(3,A,Aout);
%                 B = cat(3,B);
%             else
%                 A = cat(3,A);
%                 B = cat(3,B, Bout);
%             end
            A = cat(3,A,Aout);
            B = cat(3,B, Bout);
                
            PA = (1:size(A,3));
            PB = (1:size(B,3));

            for i = 1:length(PA)
                if rand <= 0.01*perm_rate
                    index = randi(num,1);
                    PA([i index]) = PA([index i]);
                end
            end
            A_perm = A(:, :, PA);
            B_perm = B(:, :, PB);
        
            handle_sol = usedsolver{solver_id};
            
            if strcmp(solver_name(solver_id),'B1') || strcmp(solver_name(solver_id),'B2')
                tic
                [TX, ~, ~, ~, ~, ~] = handle_sol(A_perm, B_perm, cell2mat(options(solver_id)));
                time = toc;
            else
                tic
                [TX, ~, ~, ~] = handle_sol(A_perm, B_perm);
                time = toc;
            end
            times(id1,k) = time;
            error_r(id1,k) = roterror(TX, X);
            error_t(id1,k) = tranerror(TX, X);
        end
    end
    fig = figure();
    xlbs = cellstr(string(outlier_percentage));
    box_labels = categorical(xlbs);
    box_labels = reordercats(box_labels,xlbs);
    boxplot(error_r', box_labels);
    res.error_r = error_r;
    res.error_t = error_t;
    res.times = times;
    save(strcat('C:/Users/xiahaa/Documents/MATLAB/hand_eye_calib/result/sto/test_3',solver_name{solver_id},'.mat'),'res');
end

function X = genRandomPose(num) 
    X = zeros(4,4,num);
    for i = 1:num
        x = randn(6,1); x = x./norm(x); X(:,:,i) = expm(se3_vec(x)); % Generate a Random X
    end
end