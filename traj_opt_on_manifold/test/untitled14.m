clf;clear all;close all;
figure(1);
addpath ../utils/

% for i = 1:100
%     a = (rand(3,1)); a = a ./ norm(a);
%     Rreg(:,:,i) = expSO3(a);
% end
% Rdata = Rreg(:,:,[1,30,60,100]);
% showSO3(Rdata,Rreg);

for i = 1:100
    a = expSO3(rand(3,1));
    quat(i,:) = rot2quat(a)';
    pos(i,:) = (i-1).*[0.01,0.01,0.01];
end
quatd = quat([1,30,60,100],:);
posd = pos([1,30,60,100],:);

posPlot = pos;
quatPlot = quat;

% Extend final sample to delay end of animation
extraTime = 1;
samplePeriod = 1/100;
onesVector = ones(extraTime*(1/samplePeriod), 1);
posPlot = [posPlot; [posPlot(end, 1)*onesVector, posPlot(end, 2)*onesVector, posPlot(end, 3)*onesVector]];
quatPlot = [quatPlot; [quatPlot(end, 1)*onesVector, quatPlot(end, 2)*onesVector, quatPlot(end, 3)*onesVector, quatPlot(end, 4)*onesVector]];

% Create 6 DOF animation
SamplePlotFreq = 8;
Spin = 120;
filename = 'gait.mp4';
SixDOFanimation(posPlot, quat2rotm(quatPlot), posd, quat2rotm(quatd),...
                'SamplePlotFreq', SamplePlotFreq, 'Trail', 'All', ...
                'Position', [9 39 1280 768], 'View', [(100:(Spin/(length(posPlot)-1)):(100+Spin))', 10*ones(length(posPlot), 1)], ...
                'AxisLength', 0.1, 'ShowArrowHead', false, ...
                'Xlabel', 'X (m)', 'Ylabel', 'Y (m)', 'Zlabel', 'Z (m)', 'ShowLegend', false, ...
                'CreateAVI', false, 'AVIfileNameEnum', filename, 'AVIfps', ((1/samplePeriod) / SamplePlotFreq));
            
% Use hold on and hold off to plot multiple cubes
hold on;
% Call the function to plot a cube with dimension of X, Y, Z, at point [x,y,z].
cube_plot([1,1,1],1,1,1,'r');
% Figure configurations
% Define the range of x-axis, y-axis, and z-axis in form of
% [xmin,xmax,ymin,ymax,zmin,zmax].
% axis([0,1,0,1,0,1]);
% Set the axis with equal unit.
axis equal;
% Show grids on the plot
grid on;
% Set the lable and the font size
xlabel('X','FontSize',18);
ylabel('Y','FontSize',18)
zlabel('Z','FontSize',18)
% Control the ticks on the axises
h = gca; % Get the handle of the figure
% h.XTick = 0:0.5:1;
% h.YTick = 0:0.5:1;
% h.ZTick = 0:0.5:1;
% Set the color as transparient
material metal
% alpha('color');
% alphamap('rampup');
% Set the view point
view(30,30);
hold off;
% plot the figure in the form of eps with 600 ppi named 'filename'
% print(gcf,'-depsc2','-r600','filename.eps')