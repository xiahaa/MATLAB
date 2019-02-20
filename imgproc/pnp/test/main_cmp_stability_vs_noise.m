clc;close all;clear all;

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

number_of_points = 100;%:20:50;
noise_level = 0.03:0.01:0.05;

numstd = length(noise_level);
res_error_r = cell(numstd,2);%, 100);
res_error_t = cell(numstd,2);%, 100);

for i = 1:length(number_of_points)
    num = number_of_points(i);
    for j= 1:numel(noise_level)
        nl= noise_level(j);
        clear method_list;
        load(strcat(prefix,'noise_cmp_',num2str(num),'_',num2str(nl),suffix));
%         meanerr(i,1) = mean(method_list(1).r);meanert(i,1) = mean(method_list(1).t);
%         meanerr(i,2) = mean(method_list(2).r);meanert(i,2) = mean(method_list(2).t);
        res_error_r{j,  1} = method_list(1).r;
        res_error_r{j,  2} = method_list(2).r;
        res_error_t{j,  1} = method_list(1).t;   
        res_error_t{j,  2} = method_list(2).t;
    end
end
cmap = lines(2);
box_labels = convertStringsToChars(string(noise_level));
plot_case = {'Proposed','ICP'};
cmapalpha = [cmap 0.3*ones(size(cmap,1),1)];
fig = figure();
multiple_boxplot(res_error_r, box_labels, plot_case, cmapalpha');
fontsize = 20;
xlabel('$Standard\ deviation\ of\ added\ noise$','Interpreter','latex','FontSize',fontsize);
ylabel('$E_{\mathbf{R}}$','Interpreter','latex','FontSize',fontsize);
%     legend(solver_name,'Interpreter','latex','FontSize',8,'Location', 'northeast');
title('Rotation Error\ Comparison','Interpreter','latex','FontSize',fontsize);
ylim([0,2]);
%% t
figure();
multiple_boxplot(res_error_t, box_labels, plot_case, cmapalpha');
fontsize = 20;
xlabel('$Standard\ deviation\ of\ added\ noise$','Interpreter','latex','FontSize',fontsize);
ylabel('$E_{\mathbf{t}}$','Interpreter','latex','FontSize',fontsize);
%     legend(solver_name,'Interpreter','latex','FontSize',8,'Location', 'northeast');
title('Translation Error\ Comparison','Interpreter','latex','FontSize',fontsize);
ylim([0,3]);