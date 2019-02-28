% experiment on same size, no outlier, varing noisy level and scrambling
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
stds = 0.3:0.3:1.5; % Gaussian Noise standard deviation Range
n_trials = 50; %60

perm_rate = 100;

n_noise = length(stds);

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
numout = 0.0*num;

for solver_id = 1:size(solver_name,2)
    times = zeros(numel(stds),100);
    error_r = zeros(numel(stds),n_trials);
    error_t = zeros(numel(stds),n_trials);
    valid_id = ones(numel(stds),n_trials);
    for id_std = 1:numel(stds)
        disp([solver_id,id_std])
        std = stds(id_std);
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
            [Aout, Bout] = generateAB(numout, optPDF, Xout, gmean, cov);
        
            A = cat(3,A,Aout);
            B = cat(3,B);
                
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
            times(id_std,k) = time;
            error_r(id_std,k) = roterror(TX, X);
            error_t(id_std,k) = tranerror(TX, X);
        end
    end
    fig = figure();
    xlbs = cellstr(string(stds));
    box_labels = categorical(xlbs);
    box_labels = reordercats(box_labels,xlbs);
    boxplot(error_r', box_labels);
    res.error_r = error_r;
    res.error_t = error_t;
    res.times = times;
    save(strcat('C:/Users/xiahaa/Documents/MATLAB/hand_eye_calib/result/sto/',solver_name{solver_id},'.mat'),'res');
end

function X = genRandomPose(num) 
    X = zeros(4,4,num);
    for i = 1:num
        x = randn(6,1); x = x./norm(x); X(:,:,i) = expm(se3_vec(x)); % Generate a Random X
    end
end