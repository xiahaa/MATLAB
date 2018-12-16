% load the same data and draw plots for comparing traditional and advanced
% hand eye calibration methods.
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

advSols = {'NLOPT', 'SOCP', 'GPOLY', 'DUAL', 'SCF', 'SE3OPT'};
nsols = size(advSols, 2);
N = 50; % Times of simulation

advts = zeros(numel(naaray),nsols);

advrotErrors = cell(numel(usedstd), numel(naaray));
advtranErrors = cell(numel(usedstd), numel(naaray));
advfErrors = cell(numel(usedstd), numel(naaray));

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
                        advts(i,k) = advts(i,k) + dat.tsols{kk};
                        rotError100(kk,k) = roterror(Xs{kk}, dat.xsols{kk}(:,:));
                        tranError100(kk,k) = tranerror(Xs{kk}, dat.xsols{kk}(:,:));
                        fError100(kk,k) = ferror(Xs{kk}, dat.xsols{kk}(:,:));
                    else
                        disp('no solution');
                    end
                end
            end
            
            advrotErrors{j,i} = rotError100;
            advtranErrors{j,i} = tranError100;
            advfErrors{j,i} = fError100;
        end
    end
    advts = advts ./ N;
    save('./drawing/adv.mat', 'advts', 'advrotErrors', 'advtranErrors','advfErrors');
else
    load('./drawing/adv.mat', 'advts', 'advrotErrors', 'advtranErrors','advfErrors');
end

convSols = {'TSAI', 'LIE', 'QSEP', 'KR', 'DQ', 'IDQ'};
cnsols = size(convSols, 2);
convrotErrors = cell(numel(usedstd), numel(naaray));
convtranErrors = cell(numel(usedstd), numel(naaray));
convfErrors = cell(numel(usedstd), numel(naaray));
prefix = 'data/conv/convCmp';
convts = zeros(numel(naaray),cnsols);

if 1
    for j = 1:numel(usedstd)
        for i =  1:numel(naaray)
            numPair = naaray(i);
            noisylv = num2str(usedstd(j));
            noisylv = replace(noisylv,'.','_');
            filename = strcat(prefix,'_',num2str(numPair), '_', noisylv, '.mat');
            clear dat;
            dat = load(filename);
            Xs = dat.Xs;
            flags = dat.flags;
            rotError100 = zeros(N,cnsols);
            tranError100 = zeros(N,cnsols);
            fError100 = zeros(N,cnsols);
            for kk = 1:N
                for k = 1:cnsols
                    if flags{kk}(k) == 1
                        convts(i,k) = convts(i,k) + dat.tsols{kk}(k);
                        rotError100(kk,k) = roterror(Xs{kk}, dat.xsols{kk}(:,:,k));
                        tranError100(kk,k) = tranerror(Xs{kk}, dat.xsols{kk}(:,:,k));
                        fError100(kk,k) = ferror(Xs{kk}, dat.xsols{kk}(:,:,k));
                    else
                        disp('no solution');
                    end
                end
            end
            convrotErrors{j,i} = rotError100;
            convtranErrors{j,i} = tranError100;
            convfErrors{j,i} = fError100;
        end
    end
    convts = convts ./ N;
    save('./drawing/conv.mat', 'convts', 'convrotErrors', 'convtranErrors','convfErrors');
else
    load('./drawing/conv.mat', 'convts', 'convrotErrors', 'convtranErrors','convfErrors');
end

cc = jet(12);
mm = {'-o','-*','-s','-d','-x','-+','-^'};

font_size = 10;
bar_labels = categorical(nnn);
bar_labels = reordercats(bar_labels,nnn);
fig = figure();
set(fig,'defaulttextinterpreter','latex');
for i = 1:nsols
    plot(bar_labels, advts(:,i)', mm{i}, 'LineWidth',2,...
        'MarkerEdgeColor',cc(i,:),...
        'MarkerFaceColor',cc(i,:),...
        'MarkerSize',3);hold on;
end
for i = 1:cnsols
    plot(bar_labels, convts(:,i)', mm{i}, 'LineWidth',2,...
        'MarkerEdgeColor',cc(i+nsols,:),...
        'MarkerFaceColor',cc(i+nsols,:),...
        'MarkerSize',3);hold on;
end
title('Runtime');
grid on;
xlabel('Number of Measurements','FontSize', font_size, 'Interpreter', 'latex');
ylabel('Time: (s)','FontSize', font_size, 'Interpreter', 'latex');
legend('NL', 'SOCP', 'GPOLY', 'DUAL', 'SCF', 'SE3', 'TSAI', 'LIE', 'QSEP', 'KR', 'DQ', 'IDQ', 'Interpreter','latex');
set(gca,'FontSize', font_size, 'TickLabelInterpreter','latex');


%% bar plot of error
font_size = 16;
line_width = 2.5;
subplot_margin = 0.12;
subplot_spacing = 0.1;
fig = figure();
set(fig,'defaulttextinterpreter','latex');
% subaxis(2,2,1, 'Margin', subplot_margin, 'Spacing', subplot_spacing);
x = categorical({'$0.1$', '$0.2$', '$0.3$', '$0.4$', '$0.5$', '$1$'});
x = reordercats(x,{'$0.1$', '$0.2$', '$0.3$', '$0.4$', '$0.5$', '$1$'});

ya = zeros(numel(usedstd),nsols);
yc = zeros(numel(usedstd),cnsols);
for i = 1:numel(usedstd)
    yy = advrotErrors{i,1};
    ya(i,:) = mean(yy);
    yy = convrotErrors{i,1};
    yc(i,:) = mean(yy);
end
ccc = {'r','g','b','m','c','k'};
for i = 1:nsols
    errorbar(x, ya(:,i)',std(ya(:,i)')*ones(size(x)),'Color',ccc{i},'LineWidth',2);hold on;
end
for i = 1:cnsols
    errorbar(x, yc(:,i)',std(yc(:,i)')*ones(size(x)),'--','Color',cc(i+nsols,:), 'LineWidth',2);hold on;
end
title(['N=',num2str(naaray(1))]);
set(fig,'defaulttextinterpreter','latex');
legend('NL', 'SOCP', 'GPOLY', 'DUAL', 'SCF', 'SE3', 'TSAI', 'LIE', 'QSEP', 'KR', 'DQ', 'IDQ', 'Interpreter','latex');



fig = figure();
set(fig,'defaulttextinterpreter','latex');
% subaxis(2,2,1, 'Margin', subplot_margin, 'Spacing', subplot_spacing);
x = categorical({'$0.1$', '$0.2$', '$0.3$', '$0.4$', '$0.5$', '$1$'});
x = reordercats(x,{'$0.1$', '$0.2$', '$0.3$', '$0.4$', '$0.5$', '$1$'});

ya = zeros(numel(usedstd),nsols);
yc = zeros(numel(usedstd),cnsols);
for i = 1:numel(usedstd)
    yy = advrotErrors{i,2};
    ya(i,:) = mean(yy);
    yy = convrotErrors{i,2};
    yc(i,:) = mean(yy);
end
for i = 1:nsols
    errorbar(x, ya(:,i)',std(ya(:,i)')*ones(size(x)),'Color',ccc{i},'LineWidth',2);hold on;
end
for i = 1:cnsols
    errorbar(x, yc(:,i)',std(yc(:,i)')*ones(size(x)),'--','Color',cc(i+nsols,:), 'LineWidth',2);hold on;
end
title(['N=',num2str(naaray(2))]);
set(fig,'defaulttextinterpreter','latex');
legend('NL', 'SOCP', 'GPOLY', 'DUAL', 'SCF', 'SE3', 'TSAI', 'LIE', 'QSEP', 'KR', 'DQ', 'IDQ', 'Interpreter','latex');

fig = figure();
set(fig,'defaulttextinterpreter','latex');
% subaxis(2,2,1, 'Margin', subplot_margin, 'Spacing', subplot_spacing);
x = categorical({'$0.1$', '$0.2$', '$0.3$', '$0.4$', '$0.5$', '$1$'});
x = reordercats(x,{'$0.1$', '$0.2$', '$0.3$', '$0.4$', '$0.5$', '$1$'});

ya = zeros(numel(usedstd),nsols);
yc = zeros(numel(usedstd),cnsols);
for i = 1:numel(usedstd)
    yy = advrotErrors{i,3};
    ya(i,:) = mean(yy);
    yy = convrotErrors{i,3};
    yc(i,:) = mean(yy);
end
for i = 1:nsols
    errorbar(x, ya(:,i)',std(ya(:,i)')*ones(size(x)),'Color',ccc{i},'LineWidth',2);hold on;
end
for i = 1:cnsols
    errorbar(x, yc(:,i)',std(yc(:,i)')*ones(size(x)),'--','Color',cc(i+nsols,:), 'LineWidth',2);hold on;
end
title(['N=',num2str(naaray(3))]);
set(fig,'defaulttextinterpreter','latex');
legend('NL', 'SOCP', 'GPOLY', 'DUAL', 'SCF', 'SE3', 'TSAI', 'LIE', 'QSEP', 'KR', 'DQ', 'IDQ', 'Interpreter','latex');

fig = figure();
set(fig,'defaulttextinterpreter','latex');
% subaxis(2,2,1, 'Margin', subplot_margin, 'Spacing', subplot_spacing);
x = categorical({'$0.1$', '$0.2$', '$0.3$', '$0.4$', '$0.5$', '$1$'});
x = reordercats(x,{'$0.1$', '$0.2$', '$0.3$', '$0.4$', '$0.5$', '$1$'});

ya = zeros(numel(usedstd),nsols);
yc = zeros(numel(usedstd),cnsols);
for i = 1:numel(usedstd)
    yy = advrotErrors{i,4};
    ya(i,:) = mean(yy);
    yy = convrotErrors{i,4};
    yc(i,:) = mean(yy);
end
for i = 1:nsols
    errorbar(x, ya(:,i)',std(ya(:,i)')*ones(size(x)),'Color',ccc{i},'LineWidth',2);hold on;
end
for i = 1:cnsols
    errorbar(x, yc(:,i)',std(yc(:,i)')*ones(size(x)),'--','Color',cc(i+nsols,:), 'LineWidth',2);hold on;
end
title(['N=',num2str(naaray(4))]);
set(fig,'defaulttextinterpreter','latex');
legend('NL', 'SOCP', 'GPOLY', 'DUAL', 'SCF', 'SE3', 'TSAI', 'LIE', 'QSEP', 'KR', 'DQ', 'IDQ', 'Interpreter','latex');





fig = figure();
set(fig,'defaulttextinterpreter','latex');
% subaxis(2,2,1, 'Margin', subplot_margin, 'Spacing', subplot_spacing);
x = categorical({'$0.1$', '$0.2$', '$0.3$', '$0.4$', '$0.5$', '$1$'});
x = reordercats(x,{'$0.1$', '$0.2$', '$0.3$', '$0.4$', '$0.5$', '$1$'});

ya = zeros(numel(usedstd),nsols);
yc = zeros(numel(usedstd),cnsols);
for i = 1:numel(usedstd)
    yy = advtranErrors{i,2};
    ya(i,:) = mean(yy);
    yy = convtranErrors{i,2};
    yc(i,:) = mean(yy);
end
for i = 1:nsols
    errorbar(x, ya(:,i)',std(ya(:,i)')*ones(size(x)),'Color',ccc{i},'LineWidth',2);hold on;
end
for i = 1:cnsols
    errorbar(x, yc(:,i)',std(yc(:,i)')*ones(size(x)),'--','Color',cc(i+nsols,:), 'LineWidth',2);hold on;
end
title(['N=',num2str(naaray(2))]);
set(fig,'defaulttextinterpreter','latex');
legend('NL', 'SOCP', 'GPOLY', 'DUAL', 'SCF', 'SE3', 'TSAI', 'LIE', 'QSEP', 'KR', 'DQ', 'IDQ', 'Interpreter','latex');

fig = figure();
set(fig,'defaulttextinterpreter','latex');
% subaxis(2,2,1, 'Margin', subplot_margin, 'Spacing', subplot_spacing);
x = categorical({'$0.1$', '$0.2$', '$0.3$', '$0.4$', '$0.5$', '$1$'});
x = reordercats(x,{'$0.1$', '$0.2$', '$0.3$', '$0.4$', '$0.5$', '$1$'});

ya = zeros(numel(usedstd),nsols);
yc = zeros(numel(usedstd),cnsols);
for i = 1:numel(usedstd)
    yy = advtranErrors{i,3};
    ya(i,:) = mean(yy);
    yy = convtranErrors{i,3};
    yc(i,:) = mean(yy);
end
for i = 1:nsols
    errorbar(x, ya(:,i)',std(ya(:,i)')*ones(size(x)),'Color',ccc{i},'LineWidth',2);hold on;
end
for i = 1:cnsols
    errorbar(x, yc(:,i)',std(yc(:,i)')*ones(size(x)),'--','Color',cc(i+nsols,:), 'LineWidth',2);hold on;
end
title(['N=',num2str(naaray(3))]);
set(fig,'defaulttextinterpreter','latex');
legend('NL', 'SOCP', 'GPOLY', 'DUAL', 'SCF', 'SE3', 'TSAI', 'LIE', 'QSEP', 'KR', 'DQ', 'IDQ', 'Interpreter','latex');

fig = figure();
set(fig,'defaulttextinterpreter','latex');
% subaxis(2,2,1, 'Margin', subplot_margin, 'Spacing', subplot_spacing);
x = categorical({'$0.1$', '$0.2$', '$0.3$', '$0.4$', '$0.5$', '$1$'});
x = reordercats(x,{'$0.1$', '$0.2$', '$0.3$', '$0.4$', '$0.5$', '$1$'});

ya = zeros(numel(usedstd),nsols);
yc = zeros(numel(usedstd),cnsols);
for i = 1:numel(usedstd)
    yy = advtranErrors{i,4};
    ya(i,:) = mean(yy);
    yy = convtranErrors{i,4};
    yc(i,:) = mean(yy);
end
for i = 1:nsols
    errorbar(x, ya(:,i)',std(ya(:,i)')*ones(size(x)),'Color',ccc{i},'LineWidth',2);hold on;
end
for i = 1:cnsols
    errorbar(x, yc(:,i)',std(yc(:,i)')*ones(size(x)),'--','Color',cc(i+nsols,:), 'LineWidth',2);hold on;
end
title(['N=',num2str(naaray(4))]);
set(fig,'defaulttextinterpreter','latex');
legend('NL', 'SOCP', 'GPOLY', 'DUAL', 'SCF', 'SE3', 'TSAI', 'LIE', 'QSEP', 'KR', 'DQ', 'IDQ', 'Interpreter','latex');





fig = figure();
set(fig,'defaulttextinterpreter','latex');
% subaxis(2,2,1, 'Margin', subplot_margin, 'Spacing', subplot_spacing);
x = categorical({'$0.1$', '$0.2$', '$0.3$', '$0.4$', '$0.5$', '$1$'});
x = reordercats(x,{'$0.1$', '$0.2$', '$0.3$', '$0.4$', '$0.5$', '$1$'});

ya = zeros(numel(usedstd),nsols);
yc = zeros(numel(usedstd),cnsols);
for i = 1:numel(usedstd)
    yy = advfErrors{i,2};
    ya(i,:) = mean(yy);
    yy = convfErrors{i,2};
    yc(i,:) = mean(yy);
end
for i = 1:nsols
    errorbar(x, ya(:,i)',std(ya(:,i)')*ones(size(x)),'Color',ccc{i},'LineWidth',2);hold on;
end
for i = 1:cnsols
    errorbar(x, yc(:,i)',std(yc(:,i)')*ones(size(x)),'--','Color',cc(i+nsols,:), 'LineWidth',2);hold on;
end
title(['N=',num2str(naaray(2))]);
set(fig,'defaulttextinterpreter','latex');
legend('NL', 'SOCP', 'GPOLY', 'DUAL', 'SCF', 'SE3', 'TSAI', 'LIE', 'QSEP', 'KR', 'DQ', 'IDQ', 'Interpreter','latex');

fig = figure();
set(fig,'defaulttextinterpreter','latex');
% subaxis(2,2,1, 'Margin', subplot_margin, 'Spacing', subplot_spacing);
x = categorical({'$0.1$', '$0.2$', '$0.3$', '$0.4$', '$0.5$', '$1$'});
x = reordercats(x,{'$0.1$', '$0.2$', '$0.3$', '$0.4$', '$0.5$', '$1$'});

ya = zeros(numel(usedstd),nsols);
yc = zeros(numel(usedstd),cnsols);
for i = 1:numel(usedstd)
    yy = advfErrors{i,3};
    ya(i,:) = mean(yy);
    yy = convfErrors{i,3};
    yc(i,:) = mean(yy);
end
for i = 1:nsols
    errorbar(x, ya(:,i)',std(ya(:,i)')*ones(size(x)),'Color',ccc{i},'LineWidth',2);hold on;
end
for i = 1:cnsols
    errorbar(x, yc(:,i)',std(yc(:,i)')*ones(size(x)),'--','Color',cc(i+nsols,:), 'LineWidth',2);hold on;
end
title(['N=',num2str(naaray(3))]);
set(fig,'defaulttextinterpreter','latex');
legend('NL', 'SOCP', 'GPOLY', 'DUAL', 'SCF', 'SE3', 'TSAI', 'LIE', 'QSEP', 'KR', 'DQ', 'IDQ', 'Interpreter','latex');

fig = figure();
set(fig,'defaulttextinterpreter','latex');
% subaxis(2,2,1, 'Margin', subplot_margin, 'Spacing', subplot_spacing);
x = categorical({'$0.1$', '$0.2$', '$0.3$', '$0.4$', '$0.5$', '$1$'});
x = reordercats(x,{'$0.1$', '$0.2$', '$0.3$', '$0.4$', '$0.5$', '$1$'});

ya = zeros(numel(usedstd),nsols);
yc = zeros(numel(usedstd),cnsols);
for i = 1:numel(usedstd)
    yy = advfErrors{i,4};
    ya(i,:) = mean(yy);
    yy = convfErrors{i,4};
    yc(i,:) = mean(yy);
end
for i = 1:nsols
    errorbar(x, ya(:,i)',std(ya(:,i)')*ones(size(x)),'Color',ccc{i},'LineWidth',2);hold on;
end
for i = 1:cnsols
    errorbar(x, yc(:,i)',std(yc(:,i)')*ones(size(x)),'--','Color',cc(i+nsols,:), 'LineWidth',2);hold on;
end
title(['N=',num2str(naaray(4))]);
set(fig,'defaulttextinterpreter','latex');
legend('NL', 'SOCP', 'GPOLY', 'DUAL', 'SCF', 'SE3', 'TSAI', 'LIE', 'QSEP', 'KR', 'DQ', 'IDQ', 'Interpreter','latex');
