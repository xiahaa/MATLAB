%% this file draw the results sdp solution on ETH asl dataset.
% box plot compares 6 hand eye calibration methods on 4 datasets:
% 1. primesense_1; 
% 2. primesense_2;
% 3. robotarm;
% 4. robotarm_sim

clc;clear all;close all;
addpath('../3rdparty/mit3dslam');
addpath('../beautiful_plot');
addpath('./solver/');
addpath('../MatrixLieGroup/barfoot_tro14');
addpath('../quaternion');
addpath('./solver/atadq');
addpath('../beautiful_plot/aboxplot');

% prefix = './data/sdp/';
prefix = './data/SE3/';

for id = 1:2
%     clear dat;
    dat(id) = load(strcat(prefix,'prime_', num2str(id),'_res','.mat'));
    numCls = numel(dat(id).tsol);
    convSols = dat(id).convSols;
    
end

cmap = lines(numCls);
cmapalpha = [cmap 0.2*ones(size(cmap,1),1)];
    
%% 1
font_size = 12;
line_width = 2.5;
subplot_margin = 0.12;
subplot_spacing = 0.03;
fig = figure();
% grid on;
set(fig,'defaulttextinterpreter','latex');
% subplot(1,2,1);
subaxis(1,2,1, 'Margin', subplot_margin, 'Spacing', subplot_spacing);
box_labels = categorical(convSols);
box_labels = reordercats(box_labels,convSols);
boxplot((dat(1).errcs1)',box_labels);grid on;
% Apply colors
h = findobj(gca,'Tag','Box');
for jj=1:length(h)
   patch(get(h(jj),'XData'),get(h(jj),'YData'),cmapalpha(jj,1:3),'FaceAlpha',cmapalpha(jj,4));
end
decorateBox();
ylabel('RMSE Rotation','FontSize',font_size);
xticklabel_rotate([],45,[],'Fontsize',font_size);

% subplot(1,2,2);
subaxis(1,2,2, 'Margin', subplot_margin, 'Spacing', subplot_spacing);
box_labels = categorical(convSols);
box_labels = reordercats(box_labels,convSols);
boxplot((dat(1).errcs2)',box_labels);grid on;
% Apply colors
h = findobj(gca,'Tag','Box');
for jj=1:length(h)
   patch(get(h(jj),'XData'),get(h(jj),'YData'),cmapalpha(jj,1:3),'FaceAlpha',cmapalpha(jj,4));
end
decorateBox();
ylabel('RMSE Translation','FontSize',font_size);
xticklabel_rotate([],45,[],'Fontsize',font_size);
[~,h1]=suplabel('Prime\_1' ,'t'); 
set(h1,'FontSize',font_size);
% ha = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0  1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
% text(0.5, 0.98,'Title')

%% 2
fig = figure();
% subplot(1,2,1);
subaxis(1,2,1, 'Margin', subplot_margin, 'Spacing', subplot_spacing);
box_labels = categorical(convSols);
box_labels = reordercats(box_labels,convSols);
boxplot(dat(2).errcs1',box_labels);grid on;
% Apply colors
h = findobj(gca,'Tag','Box');
for jj=1:length(h)
   patch(get(h(jj),'XData'),get(h(jj),'YData'),cmapalpha(jj,1:3),'FaceAlpha',cmapalpha(jj,4));
end
decorateBox();
ylabel('RMSE Rotation','FontSize',font_size);
xticklabel_rotate([],45,[],'Fontsize',font_size);

% subplot(1,2,2);
subaxis(1,2,2, 'Margin', subplot_margin, 'Spacing', subplot_spacing);
box_labels = categorical(convSols);
box_labels = reordercats(box_labels,convSols);
boxplot(dat(2).errcs2',box_labels);grid on;
% Apply colors
h = findobj(gca,'Tag','Box');
for jj=1:length(h)
   patch(get(h(jj),'XData'),get(h(jj),'YData'),cmapalpha(jj,1:3),'FaceAlpha',cmapalpha(jj,4));
end
decorateBox();
ylabel('RMSE Translation','FontSize',font_size);
xticklabel_rotate([],45,[],'Fontsize',font_size);
[~,h1]=suplabel('Prime_2' ,'t'); 
set(h1,'Interpreter','latex','FontSize',font_size);

%% 3
fig = figure();
% subplot(1,2,1);
subaxis(1,2,1, 'Margin', subplot_margin, 'Spacing', subplot_spacing);
box_labels = categorical(convSols);
box_labels = reordercats(box_labels,convSols);
boxplot(dat(3).errcs1',box_labels);grid on;
% Apply colors
h = findobj(gca,'Tag','Box');
for jj=1:length(h)
   patch(get(h(jj),'XData'),get(h(jj),'YData'),cmapalpha(jj,1:3),'FaceAlpha',cmapalpha(jj,4));
end
decorateBox();
ylabel('$RMSE\ Rotation$','Interpreter','latex','FontSize',font_size);
xticklabel_rotate([],45,[],'Fontsize',font_size);

% subplot(1,2,2);
subaxis(1,2,2, 'Margin', subplot_margin, 'Spacing', subplot_spacing);
box_labels = categorical(convSols);
box_labels = reordercats(box_labels,convSols);
boxplot(dat(3).errcs2',box_labels);grid on;
% Apply colors
h = findobj(gca,'Tag','Box');
for jj=1:length(h)
   patch(get(h(jj),'XData'),get(h(jj),'YData'),cmapalpha(jj,1:3),'FaceAlpha',cmapalpha(jj,4));
end
decorateBox();
ylabel('$RMSE\ Translation$','Interpreter','latex','FontSize',font_size);xticklabel_rotate([],45,[],'Fontsize',font_size);
[~,h1]=suplabel('$robotarm$' ,'t'); 
set(h1,'Interpreter','latex','FontSize',font_size);

%% 4
fig = figure();
% subplot(1,2,1);
subaxis(1,2,1, 'Margin', subplot_margin, 'Spacing', subplot_spacing);
box_labels = categorical(convSols);
box_labels = reordercats(box_labels,convSols);
boxplot(dat(4).errcs1',box_labels);grid on;
% Apply colors
h = findobj(gca,'Tag','Box');
for jj=1:length(h)
   patch(get(h(jj),'XData'),get(h(jj),'YData'),cmapalpha(jj,1:3),'FaceAlpha',cmapalpha(jj,4));
end
decorateBox();
ylabel('$RMSE\ Rotation$','Interpreter','latex','FontSize',font_size);
xticklabel_rotate([],45,[],'Fontsize',font_size);

% subplot(1,2,2);
subaxis(1,2,2, 'Margin', subplot_margin, 'Spacing', subplot_spacing);
box_labels = categorical(convSols);
box_labels = reordercats(box_labels,convSols);
boxplot(dat(4).errcs2',box_labels);grid on;
% Apply colors
h = findobj(gca,'Tag','Box');
for jj=1:length(h)
   patch(get(h(jj),'XData'),get(h(jj),'YData'),cmapalpha(jj,1:3),'FaceAlpha',cmapalpha(jj,4));
end
decorateBox();
ylabel('$RMSE\ Translation$','Interpreter','latex','FontSize',font_size);xticklabel_rotate([],45,[],'Fontsize',font_size);
[~,h1]=suplabel('$robotarm\_sim$' ,'t'); 
set(h1,'Interpreter','latex','FontSize',font_size);

function decorateBox()
    hlw = findall(gca,'tag','Lower Whisker');
    huw = findall(gca,'tag','Upper Whisker');
    set(hlw,'linestyle','-');
    set(huw,'linestyle','-');
    hout = findall(gca,'tag','Outliers');
    % set(hout,'marker','.');
    delete(hout);
end