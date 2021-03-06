%% this file draw the results 4 different solutions:
% 1. safe convex feasibility set on dual quaternion;
% 2. sequencial quadratic programming on dual quaternion;
% 3. safe convex feasibility set on kronecker product;
% 4. sequencial quadratic programming on kronecker product;
% generally speaking, they all perform very well, and their performance
% depends on how the noise is added and motion.
clc;
close all;
clear all;

addpath('../beautiful_plot/');

naaray = [10 20 30 40];
nstd2 = [0.1 0.2 0.3 0.4 0.5 1];
usedstd = nstd2;
prefix = 'data/SCF/scfCmp';

convSols = {'SCFQ', 'SQPQ', 'SCFKRON', 'SQPKRON'};
nsols = size(convSols, 2);
N = 50;
ts = zeros(numel(naaray),nsols);
nnn = {'5', '10', '20', '40'};

rotErrors = cell(numel(usedstd), numel(naaray));
tranErrors = cell(numel(usedstd), numel(naaray));

if 1
    for j = 2:numel(usedstd)
        for i =  1:numel(naaray)
            numPair = naaray(i);
            noisylv = num2str(usedstd(j));
            noisylv = replace(noisylv,'.','_');
            filename = strcat(prefix,'_',num2str(numPair), '_', noisylv, '.mat');
            dat = load(filename);
            
            Xs = dat.Xs;
            flags = dat.flags;
            rotError100 = zeros(N,nsols);
            tranError100 = zeros(N,nsols);
            for kk = 1:N
                for k = 1:nsols
                    if flags{kk}(k) == 1
                        ts(i,k) = ts(i,k) + dat.tsols{kk}(k);
                        rotError100(kk,k) = roterror(Xs{kk}, dat.xsols{kk}(:,:,k));
                        tranError100(kk,k) = tranerror(Xs{kk}, dat.xsols{kk}(:,:,k));
                    else
                        disp('no solution');
                    end
                end
            end
            rotErrors{j,i} = rotError100;
            tranErrors{j,i} = tranError100;
        end
    end
    ts = ts ./ N;
    save('./drawing/scf_cmp.mat', 'ts', 'rotErrors', 'tranErrors');
else
    load('./drawing/scf_cmp.mat', 'ts', 'rotErrors', 'tranErrors');
end

cc = jet(4);
mm = {'-o','-*','-s','-d','-x','-+','-^'};

font_size = 16;
bar_labels = categorical(nnn);
bar_labels = reordercats(bar_labels,nnn);

fig = figure();
set(fig,'defaulttextinterpreter','latex');
for i = 1:nsols
    h(i)=plot(bar_labels, ts(:,i)', mm{i}, 'LineWidth',2,...
    'MarkerEdgeColor',cc(i,:),...
    'MarkerFaceColor',cc(i,:),...
    'MarkerSize',10);hold on;
end
title('Runtime');
grid on;
xlabel('Number of Measurements','FontSize', font_size, 'Interpreter', 'latex');
ylabel('Time: (s)','FontSize', font_size, 'Interpreter', 'latex');
legend('SCFQ', 'SQPQ', 'SCFKRON', 'SQPKRON');
set(gca,'FontSize', font_size, 'TickLabelInterpreter','latex');
%% bar plot of error
font_size = 16;
line_width = 2.5;
subplot_margin = 0.12;
subplot_spacing = 0.1;
%% compare different noise level
fig = figure();
set(fig,'defaulttextinterpreter','latex');
subaxis(6,1,1, 'Margin', subplot_margin, 'Spacing', subplot_spacing);
boxdata = rotErrors{1,end};
box_labels = categorical(convSols);
box_labels = reordercats(box_labels,convSols);
boxplot(boxdata,box_labels);grid on;
title('Gaussian Noise standard deviation: (0.0), Samples: 100')
ylabel('$E_{R_X}$','Interpreter','latex');

subaxis(6,1,2, 'Margin', subplot_margin, 'Spacing', subplot_spacing);
boxdata = rotErrors{2,end};grid on;
boxplot(boxdata,box_labels);
title('Gaussian Noise standard deviation: (0.05), Samples: 100')
ylabel('$E_{R_X}$','Interpreter','latex');

subaxis(6,1,3, 'Margin', subplot_margin, 'Spacing', subplot_spacing);
boxdata = rotErrors{3,end};grid on;
boxplot(boxdata,box_labels);
title('Gaussian Noise standard deviation: (0.1), Samples: 100')
ylabel('$E_{R_X}$','Interpreter','latex');

subaxis(6,1,4, 'Margin', subplot_margin, 'Spacing', subplot_spacing);
boxdata = rotErrors{4,end};grid on;
boxplot(boxdata,box_labels);
title('Gaussian Noise standard deviation: (0.2), Samples: 100')
ylabel('$E_{R_X}$','Interpreter','latex');

subaxis(6,1,5, 'Margin', subplot_margin, 'Spacing', subplot_spacing);
boxdata = rotErrors{5,end};grid on;
boxplot(boxdata,box_labels);
title('Gaussian Noise standard deviation: (0.5), Samples: 100')
ylabel('$E_{R_X}$','Interpreter','latex');

subaxis(6,1,6, 'Margin', subplot_margin, 'Spacing', subplot_spacing);
boxdata = rotErrors{6,end};grid on;
boxplot(boxdata,box_labels);
title('Gaussian Noise standard deviation: (1), Samples: 100')
ylabel('$E_{R_X}$','Interpreter','latex');

%% compare different sample size
fig = figure();
set(fig,'defaulttextinterpreter','latex');
subaxis(6,1,1, 'Margin', subplot_margin, 'Spacing', subplot_spacing);
boxdata = rotErrors{3,1};grid on;
boxplot(boxdata,box_labels);
title('Gaussian Noise standard deviation: (0.1), Samples: 10')
ylabel('$E_{R_X}$','Interpreter','latex');

subaxis(6,1,2, 'Margin', subplot_margin, 'Spacing', subplot_spacing);
boxdata = rotErrors{3,2};grid on;
boxplot(boxdata,box_labels);
title('Gaussian Noise standard deviation: (0.1), Samples: 20')
ylabel('$E_{R_X}$','Interpreter','latex');

subaxis(6,1,3, 'Margin', subplot_margin, 'Spacing', subplot_spacing);
boxdata = rotErrors{3,3};grid on;
boxplot(boxdata,box_labels);
title('Gaussian Noise standard deviation: (0.1), Samples: 40')
ylabel('$E_{R_X}$','Interpreter','latex');

subaxis(6,1,4, 'Margin', subplot_margin, 'Spacing', subplot_spacing);
boxdata = rotErrors{3,4};grid on;
boxplot(boxdata,box_labels);
title('Gaussian Noise standard deviation: (0.1), Samples: 60')
ylabel('$E_{R_X}$','Interpreter','latex');

% subaxis(6,1,5, 'Margin', subplot_margin, 'Spacing', subplot_spacing);
% boxdata = rotErrors{3,5};grid on;
% boxplot(boxdata,box_labels);
% title('Gaussian Noise standard deviation: (0.1), Samples: 80')
% ylabel('$E_{R_X}$','Interpreter','latex');
% 
% subaxis(6,1,6, 'Margin', subplot_margin, 'Spacing', subplot_spacing);
% boxdata = rotErrors{3,6};grid on;
% boxplot(boxdata,box_labels);
% title('Gaussian Noise standard deviation: (0.1), Samples: 100')
% ylabel('$E_{R_X}$','Interpreter','latex');


%% t error
font_size = 16;
line_width = 2.5;
subplot_margin = 0.12;
subplot_spacing = 0.1;
fig = figure();
set(fig,'defaulttextinterpreter','latex');
subaxis(6,1,1, 'Margin', subplot_margin, 'Spacing', subplot_spacing);
boxdata = tranErrors{1,end};grid on;
boxplot(boxdata,box_labels);
title('Gaussian Noise standard deviation: (0.0)')
ylabel('$E_{t_X}$: (m)','Interpreter','latex');

subaxis(6,1,2, 'Margin', subplot_margin, 'Spacing', subplot_spacing);
boxdata = tranErrors{2,end};grid on;
boxplot(boxdata,box_labels);
title('Gaussian Noise standard deviation: (0.05)')
ylabel('$E_{t_X}$: (m)','Interpreter','latex');

subaxis(6,1,3, 'Margin', subplot_margin, 'Spacing', subplot_spacing);
boxdata = tranErrors{3,end};grid on;
boxplot(boxdata,box_labels);
title('Gaussian Noise standard deviation: (0.1)')
ylabel('$E_{t_X}$: (m)','Interpreter','latex');

subaxis(6,1,4, 'Margin', subplot_margin, 'Spacing', subplot_spacing);
boxdata = tranErrors{4,end};grid on;
boxplot(boxdata,box_labels);
title('Gaussian Noise standard deviation: (0.2)')
ylabel('$E_{t_X}$: (m)','Interpreter','latex');

subaxis(6,1,5, 'Margin', subplot_margin, 'Spacing', subplot_spacing);
boxdata = tranErrors{5,end};grid on;
boxplot(boxdata,box_labels);
title('Gaussian Noise standard deviation: (0.5)')
ylabel('$E_{t_X}$: (m)','Interpreter','latex');

subaxis(6,1,6, 'Margin', subplot_margin, 'Spacing', subplot_spacing);
boxdata = tranErrors{6,end};grid on;
boxplot(boxdata,box_labels);
title('Gaussian Noise standard deviation: (0.8)')
ylabel('$E_{t_X}$: (m)','Interpreter','latex');

%% compare different sample size
fig = figure();
set(fig,'defaulttextinterpreter','latex');
subaxis(6,1,1, 'Margin', subplot_margin, 'Spacing', subplot_spacing);
boxdata = tranErrors{3,1};grid on;
boxplot(boxdata,box_labels);
title('Gaussian Noise standard deviation: (0.1), Samples: 10')
ylabel('$E_{R_X}$','Interpreter','latex');

subaxis(6,1,2, 'Margin', subplot_margin, 'Spacing', subplot_spacing);
boxdata = tranErrors{3,2};grid on;
boxplot(boxdata,box_labels);
title('Gaussian Noise standard deviation: (0.1), Samples: 20')
ylabel('$E_{R_X}$','Interpreter','latex');

subaxis(6,1,3, 'Margin', subplot_margin, 'Spacing', subplot_spacing);
boxdata = tranErrors{3,3};grid on;
boxplot(boxdata,box_labels);
title('Gaussian Noise standard deviation: (0.1), Samples: 40')
ylabel('$E_{R_X}$','Interpreter','latex');

subaxis(6,1,4, 'Margin', subplot_margin, 'Spacing', subplot_spacing);
boxdata = tranErrors{3,4};grid on;
boxplot(boxdata,box_labels);
title('Gaussian Noise standard deviation: (0.1), Samples: 60')
ylabel('$E_{R_X}$','Interpreter','latex');

% subaxis(6,1,5, 'Margin', subplot_margin, 'Spacing', subplot_spacing);
% boxdata = tranErrors{3,5};grid on;
% boxplot(boxdata,box_labels);
% title('Gaussian Noise standard deviation: (0.1), Samples: 80')
% ylabel('$E_{R_X}$','Interpreter','latex');
% 
% subaxis(6,1,6, 'Margin', subplot_margin, 'Spacing', subplot_spacing);
% boxdata = tranErrors{3,6};grid on;
% boxplot(boxdata,box_labels);
% title('Gaussian Noise standard deviation: (0.1), Samples: 100')
% ylabel('$E_{R_X}$','Interpreter','latex');

