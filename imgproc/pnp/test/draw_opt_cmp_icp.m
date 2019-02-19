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

number_of_points = 10:10:100;%:20:50;

meanerr = zeros(length(number_of_points),2);
stderr = zeros(length(number_of_points),2);
meanert = zeros(length(number_of_points),2);
stdert = zeros(length(number_of_points),2);
for i = 1:length(number_of_points)
    clear method_list;
    load(strcat(prefix,'optimality_cmp_',num2str(number_of_points(i)),suffix));
    meanerr(i,1) = mean(method_list(1).r);meanert(i,1) = mean(method_list(1).t);
    meanerr(i,2) = mean(method_list(2).r);meanert(i,2) = mean(method_list(2).t);
    stderr(i,1) = std(method_list(1).r);stdert(i,1) = std(method_list(1).t);
    stderr(i,2) = std(method_list(2).r);stdert(i,2) = std(method_list(2).t);
end

fontsize = 25;
figure;
cmap = lines(2);
cmapalpha = [cmap 0.8*ones(size(cmap,1),1)];
plot(number_of_points, meanerr(:,1)','Color', cmapalpha(1,:),'LineWidth', 2.5);hold on;
plot(number_of_points, meanerr(:,2)','Color', cmapalpha(2,:),'LineWidth', 2.5);
grid on;
title('Comparison of Rotation Error','Interpreter','latex','FontSize',fontsize);
xlabel('$Number\ of\ samples$','Interpreter','latex','FontSize',fontsize);
ylabel('$E_{\mathbf{R}}: (degree)$','Interpreter','latex','FontSize',fontsize);
legend({'Proposed','ICP'},'Interpreter','latex','FontSize',fontsize,'Location', 'northwest');
yticks([0 50 100 150])

figure;
cmap = lines(2);
cmapalpha = [cmap 0.8*ones(size(cmap,1),1)];
plot(number_of_points, meanert(:,1)','Color', cmapalpha(1,:),'LineWidth', 2.5);hold on;
plot(number_of_points, meanert(:,2)','Color', cmapalpha(2,:),'LineWidth', 2.5);
grid on;
title('Comparison of Translation Error','Interpreter','latex','FontSize',fontsize);
xlabel('$Number\ of\ samples$','Interpreter','latex','FontSize',fontsize);
ylabel('$E_{\mathbf{t}}: (m)$','Interpreter','latex','FontSize',fontsize);
legend({'Proposed','ICP'},'Interpreter','latex','FontSize',fontsize,'Location', 'northwest');
% yticks([0 50 100 150])