% draw plots for comparing advanced hand eye calibration methods.
clc;
close all;
clear all;

addpath('../beautiful_plot/');

% naaray = [10 20 30 40 50 60 70 80 90 100];
% naaray = [10 20 40 60 80 100];
naaray = [10 20 30 40];


nstd = [0 0.01 0.05 0.1 0.2 0.5 0.8 1];  %Gaussian Noise standard deviation Range
nstd1 = [0 0.01 0.05 0.1 0.2 0.5 0.8 1];  %Gaussian Noise standard deviation Range
% nstd2 = [0 0.05 0.1 0.2 0.5 1];
nstd2 = [0.1 0.2 0.3 0.4 0.5 1];

usedstd = nstd2;

prefix = 'data/adv/';
folders = {'nlopt','socp','gpoly','dual','scf','se3'};
nnn = {'5', '10', '20', '40'};

convSols = {'NLOPT', 'SOCP', 'GPOLY', 'DUAL', 'SCF', 'SE3OPT'};
nsols = size(convSols, 2);
N = 50; % Times of simulation

ts = zeros(numel(naaray),nsols);

rotErrors = cell(numel(usedstd), numel(naaray));
tranErrors = cell(numel(usedstd), numel(naaray));
fErrors = cell(numel(usedstd), numel(naaray));

if 1
    for j = 1:numel(usedstd)
        for i =  1:numel(naaray)
            numPair = naaray(i);
            noisylv = num2str(usedstd(j));
            noisylv = replace(noisylv,'.','_');

            rotError100 = zeros(N,nsols);
            tranError100 = zeros(N,nsols);
            fError100 = zeros(N,nsols);

            for k = 1:nsols
                clear dat;
                filename = strcat(prefix, folders{k},'/convCmp_',num2str(numPair), '_', noisylv, '.mat');
                dat = load(filename);
                Xs = dat.Xs;
                flags = dat.flags;
                for kk = 1:N
                    if flags{kk} == 1
                        ts(i,k) = ts(i,k) + dat.tsols{kk};
                        rotError100(kk,k) = roterror(Xs{kk}, dat.xsols{kk}(:,:));
                        tranError100(kk,k) = tranerror(Xs{kk}, dat.xsols{kk}(:,:));
                        fError100(kk,k) = f_rot_error(Xs{kk}, dat.xsols{kk}(:,:));
                    else
                        disp('no solution');
                    end
                end
            end

            rotErrors{j,i} = rotError100;
            tranErrors{j,i} = tranError100;
            fErrors{j,i} = fError100;
        end
    end
    ts = ts ./ N;
    save('./drawing/adv.mat', 'ts', 'rotErrors', 'tranErrors');
else
    load('./drawing/adv.mat', 'ts', 'rotErrors', 'tranErrors');
end

cc = jet(6);
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
    'MarkerSize',5);hold on;
end
title('Runtime');
grid on;
xlabel('Number of Measurements','FontSize', font_size, 'Interpreter', 'latex');
ylabel('Time: (s)','FontSize', font_size, 'Interpreter', 'latex');
legend('NL', 'SOCP', 'GPOLY', 'DUAL', 'SCF', 'SE3','Interpreter','latex');
set(gca,'FontSize', font_size, 'TickLabelInterpreter','latex');
%% bar plot of error
font_size = 16;
line_width = 2.5;
subplot_margin = 0.12;
subplot_spacing = 0.1;
fig = figure();
set(fig,'defaulttextinterpreter','latex');
subaxis(3,2,1, 'Margin', subplot_margin, 'Spacing', subplot_spacing);
boxdata = rotErrors{1,end};
box_labels = categorical(convSols);
box_labels = reordercats(box_labels,convSols);
boxplot(boxdata,box_labels);
title('$\sigma$: (0.0), N = 100')
ylabel('$E_{R_X}$','Interpreter','latex');
grid on;

subaxis(3,2,2, 'Margin', subplot_margin, 'Spacing', subplot_spacing);
boxdata = rotErrors{2,end};
boxplot(boxdata,box_labels);
title('$\sigma$: (0.01), N = 100')
ylabel('$E_{R_X}$','Interpreter','latex');
grid on;

subaxis(3,2,3, 'Margin', subplot_margin, 'Spacing', subplot_spacing);
boxdata = rotErrors{3,end};
boxplot(boxdata,box_labels);
title('$\sigma$: (0.1), N = 100')
ylabel('$E_{R_X}$','Interpreter','latex');
grid on;

subaxis(3,2,4, 'Margin', subplot_margin, 'Spacing', subplot_spacing);
boxdata = rotErrors{4,end};
boxplot(boxdata,box_labels);
title('$\sigma$: (0.2), N = 100')
ylabel('$E_{R_X}$','Interpreter','latex');
grid on;

subaxis(3,2,5, 'Margin', subplot_margin, 'Spacing', subplot_spacing);
boxdata = rotErrors{5,end};
boxplot(boxdata,box_labels);
title('$\sigma$: (0.5), N = 100')
ylabel('$E_{R_X}$','Interpreter','latex');
grid on;

subaxis(3,2,6, 'Margin', subplot_margin, 'Spacing', subplot_spacing);
boxdata = rotErrors{6,end};
boxplot(boxdata,box_labels);
title('$\sigma$: (1), N = 100')
ylabel('$E_{R_X}$','Interpreter','latex');
grid on;

% legend('TSAI', 'LIE', 'QSEP', 'KRSEP', 'DQ', 'IDQ','Interpreter','latex');
fig = figure();
set(fig,'defaulttextinterpreter','latex');
subaxis(3,2,1, 'Margin', subplot_margin, 'Spacing', subplot_spacing);
boxdata = rotErrors{3,1};
boxplot(boxdata,box_labels);
title('$\sigma$: (0.1), N = 10')
ylabel('$E_{R_X}$','Interpreter','latex');
grid on;

subaxis(3,2,2, 'Margin', subplot_margin, 'Spacing', subplot_spacing);
boxdata = rotErrors{3,2};
boxplot(boxdata,box_labels);
title('$\sigma$: (0.1), N = 20')
ylabel('$E_{R_X}$','Interpreter','latex');
grid on;

subaxis(3,2,3, 'Margin', subplot_margin, 'Spacing', subplot_spacing);
boxdata = rotErrors{3,3};
boxplot(boxdata,box_labels);
title('$\sigma$: (0.1), N = 40')
ylabel('$E_{R_X}$','Interpreter','latex');
grid on;

subaxis(3,2,4, 'Margin', subplot_margin, 'Spacing', subplot_spacing);
boxdata = rotErrors{3,4};
boxplot(boxdata,box_labels);
title('$\sigma$: (0.1), N = 60')
ylabel('$E_{R_X}$','Interpreter','latex');
grid on;

if size(rotErrors,2) > 4
    subaxis(3,2,5, 'Margin', subplot_margin, 'Spacing', subplot_spacing);
    boxdata = rotErrors{3,5};
    boxplot(boxdata,box_labels);
    title('$\sigma$: (0.1), N = 80')
    ylabel('$E_{R_X}$','Interpreter','latex');
    grid on;

    subaxis(3,2,6, 'Margin', subplot_margin, 'Spacing', subplot_spacing);
    boxdata = rotErrors{3,6};
    boxplot(boxdata,box_labels);
    title('$\sigma$: (0.1), N = 100')
    ylabel('$E_{R_X}$','Interpreter','latex');
    grid on;
end

%% t error
font_size = 16;
line_width = 2.5;
subplot_margin = 0.12;
subplot_spacing = 0.1;
fig = figure();
set(fig,'defaulttextinterpreter','latex');
subaxis(3,2,1, 'Margin', subplot_margin, 'Spacing', subplot_spacing);
boxdata = tranErrors{1,end};
boxplot(boxdata,box_labels);
title('$\sigma$: (0.0), N = 100')
ylabel('$E_{t_X}$: (m)','Interpreter','latex');
grid on;

subaxis(3,2,2, 'Margin', subplot_margin, 'Spacing', subplot_spacing);
boxdata = tranErrors{2,end};
boxplot(boxdata,box_labels);
title('$\sigma$: (0.05), N = 100')
ylabel('$E_{t_X}$: (m)','Interpreter','latex');
grid on;

subaxis(3,2,3, 'Margin', subplot_margin, 'Spacing', subplot_spacing);
boxdata = tranErrors{3,end};
boxplot(boxdata,box_labels);
title('$\sigma$: (0.1), N = 100')
ylabel('$E_{t_X}$: (m)','Interpreter','latex');
grid on;

subaxis(3,2,4, 'Margin', subplot_margin, 'Spacing', subplot_spacing);
boxdata = tranErrors{4,end};
boxplot(boxdata,box_labels);
title('$\sigma$: (0.2), N = 100')
ylabel('$E_{t_X}$: (m)','Interpreter','latex');
grid on;

subaxis(3,2,5, 'Margin', subplot_margin, 'Spacing', subplot_spacing);
boxdata = tranErrors{5,end};
boxplot(boxdata,box_labels);
title('$\sigma$: (0.5), N = 100')
ylabel('$E_{t_X}$: (m)','Interpreter','latex');
grid on;

subaxis(3,2,6, 'Margin', subplot_margin, 'Spacing', subplot_spacing);
boxdata = tranErrors{6,end};
boxplot(boxdata,box_labels);
title('$\sigma$: (1), N = 100')
ylabel('$E_{t_X}$: (m)','Interpreter','latex');
grid on;

fig = figure();
set(fig,'defaulttextinterpreter','latex');
subaxis(3,2,1, 'Margin', subplot_margin, 'Spacing', subplot_spacing);
boxdata = tranErrors{3,1};
boxplot(boxdata,box_labels);
title('$\sigma$: (0.1), N = 10')
ylabel('$E_{t_X}$','Interpreter','latex');
grid on;

subaxis(3,2,2, 'Margin', subplot_margin, 'Spacing', subplot_spacing);
boxdata = tranErrors{3,2};
boxplot(boxdata,box_labels);
title('$\sigma$: (0.1), N = 20')
ylabel('$E_{t_X}$','Interpreter','latex');
grid on;

subaxis(3,2,3, 'Margin', subplot_margin, 'Spacing', subplot_spacing);
boxdata = tranErrors{3,3};
boxplot(boxdata,box_labels);
title('$\sigma$: (0.1), N = 40')
ylabel('$E_{t_X}$','Interpreter','latex');
grid on;

subaxis(3,2,4, 'Margin', subplot_margin, 'Spacing', subplot_spacing);
boxdata = tranErrors{3,4};
boxplot(boxdata,box_labels);
title('$\sigma$: (0.1), N = 60')
ylabel('$E_{t_X}$','Interpreter','latex');
grid on;
if size(rotErrors,2) > 4
    subaxis(3,2,5, 'Margin', subplot_margin, 'Spacing', subplot_spacing);
    boxdata = tranErrors{3,5};
    boxplot(boxdata,box_labels);
    title('$\sigma$: (0.1), N = 80')
    ylabel('$E_{t_X}$','Interpreter','latex');
    grid on;

    subaxis(3,2,6, 'Margin', subplot_margin, 'Spacing', subplot_spacing);
    boxdata = tranErrors{3,6};
    boxplot(boxdata,box_labels);
    title('$\sigma$: (0.1), N = 100')
    ylabel('$E_{t_X}$','Interpreter','latex');
    grid on;
end

%% f error
font_size = 16;
line_width = 2.5;
subplot_margin = 0.12;
subplot_spacing = 0.1;
fig = figure();
set(fig,'defaulttextinterpreter','latex');
subaxis(3,2,1, 'Margin', subplot_margin, 'Spacing', subplot_spacing);
boxdata = fErrors{1,end};
boxplot(boxdata,box_labels);
title('$\sigma$: (0.0), N = 100')
ylabel('$E_{Frob}$','Interpreter','latex');
grid on;

subaxis(3,2,2, 'Margin', subplot_margin, 'Spacing', subplot_spacing);
boxdata = fErrors{2,end};
boxplot(boxdata,box_labels);
title('$\sigma$: (0.05), N = 100')
ylabel('$E_{Frob}$','Interpreter','latex');
grid on;

subaxis(3,2,3, 'Margin', subplot_margin, 'Spacing', subplot_spacing);
boxdata = fErrors{3,end};
boxplot(boxdata,box_labels);
title('$\sigma$: (0.1), N = 100')
ylabel('$E_{Frob}$','Interpreter','latex');
grid on;

subaxis(3,2,4, 'Margin', subplot_margin, 'Spacing', subplot_spacing);
boxdata = fErrors{4,end};
boxplot(boxdata,box_labels);
title('$\sigma$: (0.2), N = 100')
ylabel('$E_{Frob}$','Interpreter','latex');
grid on;

subaxis(3,2,5, 'Margin', subplot_margin, 'Spacing', subplot_spacing);
boxdata = fErrors{5,end};
boxplot(boxdata,box_labels);
title('$\sigma$: (0.5), N = 100')
ylabel('$E_{Frob}$','Interpreter','latex');
grid on;

subaxis(3,2,6, 'Margin', subplot_margin, 'Spacing', subplot_spacing);
boxdata = fErrors{6,end};
boxplot(boxdata,box_labels);
title('$\sigma$: (1), N = 100')
ylabel('$E_{Frob}$','Interpreter','latex');
grid on;

fig = figure();
set(fig,'defaulttextinterpreter','latex');
subaxis(3,2,1, 'Margin', subplot_margin, 'Spacing', subplot_spacing);
boxdata = fErrors{3,1};
boxplot(boxdata,box_labels);
title('$\sigma$: (0.1), N = 10')
ylabel('$E_{Frob}$','Interpreter','latex');
grid on;

subaxis(3,2,2, 'Margin', subplot_margin, 'Spacing', subplot_spacing);
boxdata = fErrors{3,2};
boxplot(boxdata,box_labels);
title('$\sigma$: (0.1), N = 20')
ylabel('$E_{Frob}$','Interpreter','latex');
grid on;

subaxis(3,2,3, 'Margin', subplot_margin, 'Spacing', subplot_spacing);
boxdata = fErrors{3,3};
boxplot(boxdata,box_labels);
title('$\sigma$: (0.1), N = 40')
ylabel('$E_{Frob}$','Interpreter','latex');
grid on;

subaxis(3,2,4, 'Margin', subplot_margin, 'Spacing', subplot_spacing);
boxdata = fErrors{3,4};
boxplot(boxdata,box_labels);
title('$\sigma$: (0.1), N = 60')
ylabel('$E_{Frob}$','Interpreter','latex');
grid on;

if size(rotErrors,2) > 4
    subaxis(3,2,5, 'Margin', subplot_margin, 'Spacing', subplot_spacing);
    boxdata = fErrors{3,5};
    boxplot(boxdata,box_labels);
    title('$\sigma$: (0.1), N = 80')
    ylabel('$E_{Frob}$','Interpreter','latex');
    grid on;

    subaxis(3,2,6, 'Margin', subplot_margin, 'Spacing', subplot_spacing);
    boxdata = fErrors{3,6};
    boxplot(boxdata,box_labels);
    title('$\sigma$: (0.1), N = 100')
    ylabel('$E_{Frob}$','Interpreter','latex');
    grid on;
end
