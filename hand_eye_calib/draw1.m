clc;
close all;
clear all;

addpath('../beautiful_plot/');

naaray = [10 20 30 40 50 60 70 80 90 100];

nstd = [0 0.01 0.05 0.1 0.2 0.5 0.8 1];  %Gaussian Noise standard deviation Range
nstd1 = [0 0.01 0.05 0.1 0.2 0.5 0.8 1];  %Gaussian Noise standard deviation Range

prefix = 'data/convCmp';

convSols = {'TSAI', 'LIE', 'QSEP', 'KR', 'DQ', 'IDQ'};
nsols = size(convSols, 2);
N = 100; % Times of simulation

ts = zeros(numel(naaray),nsols);

rotErrors = cell(numel(nstd1), numel(naaray));
tranErrors = cell(numel(nstd1), numel(naaray));

if 0
    for j = 1:numel(nstd1)
        for i =  1:numel(naaray)
            numPair = naaray(i);
            noisylv = num2str(nstd(j));
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
    save('./drawing/conv.mat', 'ts', 'rotErrors', 'tranErrors');
else
    load('./drawing/conv.mat', 'ts', 'rotErrors', 'tranErrors');
end

cc = jet(6);
mm = ['-o','-*','-s','-d','-x','-+','-^'];

font_size = 16;
bar_labels = categorical({'$10$', '$20$', '$30$', '$40$', '$50$', '$60$', '$70$', '$80$', '$90$', '$100$'});
bar_labels = reordercats(bar_labels,{'$10$', '$20$', '$30$', '$40$', '$50$', '$60$', '$70$', '$80$', '$90$', '$100$'});
line_width = 2.5;
subplot_margin = 0.12;
subplot_spacing = 0.16;
fig = figure();
set(fig,'defaulttextinterpreter','latex');

title('Runtime');
for i = 1:nsols
    h(i)=plot(bar_labels, ts(:,i)', 'LineWidth',2,...
    'MarkerEdgeColor',cc(i,:),...
    'MarkerFaceColor',cc(i,:),...
    'MarkerSize',10);hold on;
end
grid on;
xlabel('Number of Measurements','FontSize', font_size, 'Interpreter', 'latex');
ylabel('Time: (s)','FontSize', font_size, 'Interpreter', 'latex');
legend('TSAI', 'LIE', 'QSEP', 'KRSEP', 'DQ', 'IDQ','Interpreter','latex');
set(gca,'FontSize', font_size, 'TickLabelInterpreter','latex');
%% bar plot of error

bar_labels = categorical({'$10$', '$20$', '$30$', '$40$', '$50$', '$60$', '$70$', '$80$', '$90$', '$100$'});
bar_labels = reordercats(bar_labels,{'$10$', '$20$', '$30$', '$40$', '$50$', '$60$', '$70$', '$80$', '$90$', '$100$'});
font_size = 16;
line_width = 2.5;
subplot_margin = 0.12;
subplot_spacing = 0.05;
fig = figure();
set(fig,'defaulttextinterpreter','latex');
subaxis(8,1,1, 'Margin', subplot_margin, 'Spacing', subplot_spacing);
boxdata = rotErrors{1,10};
box_labels = categorical({'TSAI', 'LIE', 'QSEP', 'KRSEP', 'DQ', 'IDQ'});
box_labels = reordercats(box_labels,{'TSAI', 'LIE', 'QSEP', 'KRSEP', 'DQ', 'IDQ'});
boxplot(boxdata,box_labels);
title('Gaussian Noise standard deviation: (0.0)')
ylabel('$E_{R_X}$','Interpreter','latex');

subaxis(8,1,2, 'Margin', subplot_margin, 'Spacing', subplot_spacing);
boxdata = rotErrors{2,10};
box_labels = categorical({'TSAI', 'LIE', 'QSEP', 'KRSEP', 'DQ', 'IDQ'});
box_labels = reordercats(box_labels,{'TSAI', 'LIE', 'QSEP', 'KRSEP', 'DQ', 'IDQ'});
boxplot(boxdata,box_labels);
title('Gaussian Noise standard deviation: (0.01)')
ylabel('$E_{R_X}$','Interpreter','latex');

subaxis(8,1,3, 'Margin', subplot_margin, 'Spacing', subplot_spacing);
boxdata = rotErrors{3,10};
box_labels = categorical({'TSAI', 'LIE', 'QSEP', 'KRSEP', 'DQ', 'IDQ'});
box_labels = reordercats(box_labels,{'TSAI', 'LIE', 'QSEP', 'KRSEP', 'DQ', 'IDQ'});
boxplot(boxdata,box_labels);
title('Gaussian Noise standard deviation: (0.05)')
ylabel('$E_{R_X}$','Interpreter','latex');

subaxis(8,1,4, 'Margin', subplot_margin, 'Spacing', subplot_spacing);
boxdata = rotErrors{4,10};
box_labels = categorical({'TSAI', 'LIE', 'QSEP', 'KRSEP', 'DQ', 'IDQ'});
box_labels = reordercats(box_labels,{'TSAI', 'LIE', 'QSEP', 'KRSEP', 'DQ', 'IDQ'});
boxplot(boxdata,box_labels);
title('Gaussian Noise standard deviation: (0.1)')
ylabel('$E_{R_X}$','Interpreter','latex');

subaxis(8,1,5, 'Margin', subplot_margin, 'Spacing', subplot_spacing);
boxdata = rotErrors{5,10};
box_labels = categorical({'TSAI', 'LIE', 'QSEP', 'KRSEP', 'DQ', 'IDQ'});
box_labels = reordercats(box_labels,{'TSAI', 'LIE', 'QSEP', 'KRSEP', 'DQ', 'IDQ'});
boxplot(boxdata,box_labels);
title('Gaussian Noise standard deviation: (0.2)')
ylabel('$E_{R_X}$','Interpreter','latex');

subaxis(8,1,6, 'Margin', subplot_margin, 'Spacing', subplot_spacing);
boxdata = rotErrors{6,10};
box_labels = categorical({'TSAI', 'LIE', 'QSEP', 'KRSEP', 'DQ', 'IDQ'});
box_labels = reordercats(box_labels,{'TSAI', 'LIE', 'QSEP', 'KRSEP', 'DQ', 'IDQ'});
boxplot(boxdata,box_labels);
title('Gaussian Noise standard deviation: (0.5)')
ylabel('$E_{R_X}$','Interpreter','latex');

subaxis(8,1,7, 'Margin', subplot_margin, 'Spacing', subplot_spacing);
boxdata = rotErrors{7,10};
box_labels = categorical({'TSAI', 'LIE', 'QSEP', 'KRSEP', 'DQ', 'IDQ'});
box_labels = reordercats(box_labels,{'TSAI', 'LIE', 'QSEP', 'KRSEP', 'DQ', 'IDQ'});
boxplot(boxdata,box_labels);
title('Gaussian Noise standard deviation: (0.8)')
ylabel('$E_{R_X}$','Interpreter','latex');

subaxis(8,1,8, 'Margin', subplot_margin, 'Spacing', subplot_spacing);
boxdata = rotErrors{8,10};
box_labels = categorical({'TSAI', 'LIE', 'QSEP', 'KRSEP', 'DQ', 'IDQ'});
box_labels = reordercats(box_labels,{'TSAI', 'LIE', 'QSEP', 'KRSEP', 'DQ', 'IDQ'});
boxplot(boxdata,box_labels);
title('Gaussian Noise standard deviation: (1)')
ylabel('$E_{R_X}$','Interpreter','latex');

% legend('TSAI', 'LIE', 'QSEP', 'KRSEP', 'DQ', 'IDQ','Interpreter','latex');


%% t error
bar_labels = categorical({'$10$', '$20$', '$30$', '$40$', '$50$', '$60$', '$70$', '$80$', '$90$', '$100$'});
bar_labels = reordercats(bar_labels,{'$10$', '$20$', '$30$', '$40$', '$50$', '$60$', '$70$', '$80$', '$90$', '$100$'});
font_size = 16;
line_width = 2.5;
subplot_margin = 0.12;
subplot_spacing = 0.05;
fig = figure();
set(fig,'defaulttextinterpreter','latex');
subaxis(8,1,1, 'Margin', subplot_margin, 'Spacing', subplot_spacing);
boxdata = tranErrors{1,10};
box_labels = categorical({'TSAI', 'LIE', 'QSEP', 'KRSEP', 'DQ', 'IDQ'});
box_labels = reordercats(box_labels,{'TSAI', 'LIE', 'QSEP', 'KRSEP', 'DQ', 'IDQ'});
boxplot(boxdata,box_labels);
title('Gaussian Noise standard deviation: (0.0)')
ylabel('$E_{t_X}$: (m)','Interpreter','latex');

subaxis(8,1,2, 'Margin', subplot_margin, 'Spacing', subplot_spacing);
boxdata = tranErrors{2,10};
box_labels = categorical({'TSAI', 'LIE', 'QSEP', 'KRSEP', 'DQ', 'IDQ'});
box_labels = reordercats(box_labels,{'TSAI', 'LIE', 'QSEP', 'KRSEP', 'DQ', 'IDQ'});
boxplot(boxdata,box_labels);
title('Gaussian Noise standard deviation: (0.01)')
ylabel('$E_{t_X}$: (m)','Interpreter','latex');

subaxis(8,1,3, 'Margin', subplot_margin, 'Spacing', subplot_spacing);
boxdata = tranErrors{3,10};
box_labels = categorical({'TSAI', 'LIE', 'QSEP', 'KRSEP', 'DQ', 'IDQ'});
box_labels = reordercats(box_labels,{'TSAI', 'LIE', 'QSEP', 'KRSEP', 'DQ', 'IDQ'});
boxplot(boxdata,box_labels);
title('Gaussian Noise standard deviation: (0.05)')
ylabel('$E_{t_X}$: (m)','Interpreter','latex');

subaxis(8,1,4, 'Margin', subplot_margin, 'Spacing', subplot_spacing);
boxdata = tranErrors{4,10};
box_labels = categorical({'TSAI', 'LIE', 'QSEP', 'KRSEP', 'DQ', 'IDQ'});
box_labels = reordercats(box_labels,{'TSAI', 'LIE', 'QSEP', 'KRSEP', 'DQ', 'IDQ'});
boxplot(boxdata,box_labels);
title('Gaussian Noise standard deviation: (0.1)')
ylabel('$E_{t_X}$: (m)','Interpreter','latex');

subaxis(8,1,5, 'Margin', subplot_margin, 'Spacing', subplot_spacing);
boxdata = tranErrors{5,10};
box_labels = categorical({'TSAI', 'LIE', 'QSEP', 'KRSEP', 'DQ', 'IDQ'});
box_labels = reordercats(box_labels,{'TSAI', 'LIE', 'QSEP', 'KRSEP', 'DQ', 'IDQ'});
boxplot(boxdata,box_labels);
title('Gaussian Noise standard deviation: (0.2)')
ylabel('$E_{t_X}$: (m)','Interpreter','latex');


subaxis(8,1,6, 'Margin', subplot_margin, 'Spacing', subplot_spacing);
boxdata = tranErrors{6,10};
box_labels = categorical({'TSAI', 'LIE', 'QSEP', 'KRSEP', 'DQ', 'IDQ'});
box_labels = reordercats(box_labels,{'TSAI', 'LIE', 'QSEP', 'KRSEP', 'DQ', 'IDQ'});
boxplot(boxdata,box_labels);
title('Gaussian Noise standard deviation: (0.5)')
ylabel('$E_{t_X}$: (m)','Interpreter','latex');


subaxis(8,1,7, 'Margin', subplot_margin, 'Spacing', subplot_spacing);
boxdata = tranErrors{7,10};
box_labels = categorical({'TSAI', 'LIE', 'QSEP', 'KRSEP', 'DQ', 'IDQ'});
box_labels = reordercats(box_labels,{'TSAI', 'LIE', 'QSEP', 'KRSEP', 'DQ', 'IDQ'});
boxplot(boxdata,box_labels);
title('Gaussian Noise standard deviation: (0.8)')
ylabel('$E_{t_X}$: (m)','Interpreter','latex');


subaxis(8,1,8, 'Margin', subplot_margin, 'Spacing', subplot_spacing);
boxdata = tranErrors{8,10};
box_labels = categorical({'TSAI', 'LIE', 'QSEP', 'KRSEP', 'DQ', 'IDQ'});
box_labels = reordercats(box_labels,{'TSAI', 'LIE', 'QSEP', 'KRSEP', 'DQ', 'IDQ'});
boxplot(boxdata,box_labels);
title('Gaussian Noise standard deviation: (1)')
ylabel('$E_{t_X}$: (m)','Interpreter','latex');
