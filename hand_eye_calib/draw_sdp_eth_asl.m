clc;clear all;close all;
addpath('../3rdparty/mit3dslam');
addpath('../beautiful_plot');
addpath('./solver/');
addpath('../MatrixLieGroup/barfoot_tro14');
addpath('../quaternion');
addpath('./solver/atadq');
addpath('../beautiful_plot/aboxplot');

for id = 1:4
    clear dat;
    dat(id) = load(strcat('./data/sdp/','prime_', num2str(id),'_res','.mat'));
    numCls = numel(dat(id).tsol);
    convSols = dat(id).convSols;
    
end
% res_error_r = cell(1,numCls);
% res_error_t = cell(1,numCls);
% res_error_T = cell(1,numCls);
% 
% for i = 1:numCls
%     res_error_r{1,i} = errcs1(i,:);
%     res_error_t{1,i} = errcs2(i,:);
%     res_error_T{1,i} = errcs3(i,:);
% end
cmap = lines(numCls);
cmapalpha = [cmap 0.2*ones(size(cmap,1),1)];
    
%% 1
font_size = 15;
line_width = 2.5;
subplot_margin = 0.12;
subplot_spacing = 0.1;
fig = figure();
% grid on;
set(fig,'defaulttextinterpreter','latex');
subplot(4,2,1);
% subaxis(4,2,1, 'Margin', subplot_margin, 'Spacing', subplot_spacing);
box_labels = categorical(convSols);
box_labels = reordercats(box_labels,convSols);
boxplot((dat(1).errcs1)',box_labels);grid on;
% Apply colors
h = findobj(gca,'Tag','Box');
for jj=1:length(h)
   patch(get(h(jj),'XData'),get(h(jj),'YData'),cmapalpha(jj,1:3),'FaceAlpha',cmapalpha(jj,4));
end
decorateBox();
ylabel('$RMSE\ Rotation$','Interpreter','latex','FontSize',font_size);

subplot(4,2,2);
% subaxis(4,2,2, 'Margin', subplot_margin, 'Spacing', subplot_spacing);
box_labels = categorical(convSols);
box_labels = reordercats(box_labels,convSols);
boxplot((dat(1).errcs2)',box_labels);grid on;
% Apply colors
h = findobj(gca,'Tag','Box');
for jj=1:length(h)
   patch(get(h(jj),'XData'),get(h(jj),'YData'),cmapalpha(jj,1:3),'FaceAlpha',cmapalpha(jj,4));
end
decorateBox();
ylabel('$RMSE\ Translation$','Interpreter','latex','FontSize',font_size);

% ha = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0  1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
% text(0.5, 0.98,'Title')

%% 2
subplot(4,2,3);
% subaxis(4,2,3, 'Margin', subplot_margin, 'Spacing', subplot_spacing);
box_labels = categorical(convSols);
box_labels = reordercats(box_labels,convSols);
boxplot(dat(2).errcs1',box_labels);grid on;
% Apply colors
h = findobj(gca,'Tag','Box');
for jj=1:length(h)
   patch(get(h(jj),'XData'),get(h(jj),'YData'),cmapalpha(jj,1:3),'FaceAlpha',cmapalpha(jj,4));
end
decorateBox();
ylabel('$RMSE\ Rotation$','Interpreter','latex','FontSize',font_size);

subplot(4,2,4);
% subaxis(4,2,4, 'Margin', subplot_margin, 'Spacing', subplot_spacing);
box_labels = categorical(convSols);
box_labels = reordercats(box_labels,convSols);
boxplot(dat(2).errcs2',box_labels);grid on;
% Apply colors
h = findobj(gca,'Tag','Box');
for jj=1:length(h)
   patch(get(h(jj),'XData'),get(h(jj),'YData'),cmapalpha(jj,1:3),'FaceAlpha',cmapalpha(jj,4));
end
decorateBox();
ylabel('$RMSE\ Translation$','Interpreter','latex','FontSize',font_size);

%% 3
% subplot(4,2,3);
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
subaxis(4,2,4, 'Margin', subplot_margin, 'Spacing', subplot_spacing);
box_labels = categorical(convSols);
box_labels = reordercats(box_labels,convSols);
boxplot(dat(3).errcs2',box_labels);grid on;
% Apply colors
h = findobj(gca,'Tag','Box');
for jj=1:length(h)
   patch(get(h(jj),'XData'),get(h(jj),'YData'),cmapalpha(jj,1:3),'FaceAlpha',cmapalpha(jj,4));
end
decorateBox();
ylabel('$RMSE\ Translation$','Interpreter','latex','FontSize',font_size);

%% 4
subplot(4,1,4);
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
ylabel('$RMSE\ Translation$','Interpreter','latex','FontSize',font_size);

function decorateBox()
    hlw = findall(gca,'tag','Lower Whisker');
    huw = findall(gca,'tag','Upper Whisker');
    set(hlw,'linestyle','-');
    set(huw,'linestyle','-');
    hout = findall(gca,'tag','Outliers');
    % set(hout,'marker','.');
    delete(hout);
end