% This main fucntion compares three Batch methods with three 
% corresponding definitions (representations) of mean and the Kronecker
% product method
%
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
num = 100; % Number of measurements
gmean = [0;0;0;0;0;0];	%Gaussian Noise Mean

nstd = 0.1:0.1:1.0; % Gaussian Noise standard deviation Range

n_trials = 50; %60

x = randn(6,1); x = x./norm(x); X = expm(se3_vec(x)); % Generate a Random X

while 1
    xout = randn(6,1); 
    xout = xout./norm(xout);
    if norm(xout-x) > 0.8
        break;
    end
end
Xout = expm(se3_vec(xout)); % Generate a Random Xout

skip = 10; %10;
perm_rate = 100;

n_rate = length(perm_rate);
n_num = length(nstd);

%% Initialize error containers
runtime_batch_1 = zeros(n_trials, n_num);
runtime_batch_2 = zeros(n_trials, n_num);
runtime_batch_3 = zeros(n_trials, n_num);


%% Computation Loops
A_noise = [];
B = [];
opt = 2;
optPDF = 1; % select the distribution for A and B

for j = 1:n_num   
    cov = nstd(j).*eye(6,6);
    for k = 1:n_trials
        numout = 0.5*num;
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
        
        tic
        [X_batch_New_1, ~, ~, ~, ~, ~] = batchSolveNew(A_perm, B_perm, 1);
        runtime_batch_1(k,j) = toc;
        tic
        [X_batch_New_2, MA_New, MB_New, ~, ~, t_error] = batchSolveNew(A_perm, B_perm, 4); % 2 is changing step size ,4 is M=M(I+X)
        runtime_batch_2(k,j) = toc;
%         [X_batch_New_3, MA, MB, SA, SB, ~] = batchSolveNew(A_perm, B_perm, 5);
        tic
        [X_batch_New_4, ~, ~, ~] = batchSolveSoftUseScrew(A_perm, B_perm);
        runtime_batch_3(k,j) = toc;
    end
end

%% compute mean time
meantimes = zeros(3, n_num);
meantimes(1,:) = mean(runtime_batch_1);
meantimes(2,:) = mean(runtime_batch_2);
meantimes(3,:) = mean(runtime_batch_3);

figure()
cmap = lines(3);
cmapalpha = [cmap 0.8*ones(size(cmap,1),1)];
xlabels = nstd;
plot_case = {'B1','B2','BS'};
for solver_id = 1:3
    plot(xlabels', meantimes(solver_id,:),'-o', 'Color', cmapalpha(solver_id,:),'LineWidth', 2.5);hold on;
end
fontsize = 16;
xlabel('Standard deviation of additional noise','FontSize',fontsize);
ylabel('Time: (s)','FontSize',fontsize);
legend(plot_case,'FontSize',fontsize,'Location', 'northwest');
title('Runtime Comparison','FontSize',fontsize);
grid on;
yticks([0 0.1 0.2 0.3 0.4 0.5 0.6])
ylim([0 0.6])


function X = genRandomPose(num) 
    X = zeros(4,4,num);
    for i = 1:num
        x = randn(6,1); x = x./norm(x); X(:,:,i) = expm(se3_vec(x)); % Generate a Random X
    end
end