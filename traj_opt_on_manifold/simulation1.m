if ismac
    addpath './SOn_regression-master/'
    addpath './SOn_regression-master/STL'
    addpath './utils/'
    addpath './libso3'
else
    addpath './SOn_regression-master\SOn_regression-master/'
    addpath './jiachao/'
    addpath './utils/'
end
clc;close all;clear all;
warning off;
    
lat = linspace(-pi/4,pi/3,10);
lon = linspace(-pi,pi,10);
r = zeros(3,10);
data.n = 3;
data.N = 10;
data.p = zeros(3,3,10);
for i = 1:10
    r(:,i) = [cos(lat(i))*cos(lon(i));cos(lat(i))*sin(lon(i));sin(lat(i))];
    data.p(:,:,i) = expSO3(r(:,i));
end
    
n = data.n;
N = data.N;
p = data.p;
w = ones(N, 1);
Nd = 100;
s = round(linspace(1, Nd, N));
delta_tau = 1/(Nd-1);
lambda = 0;
mu = 1e-2;

%% Pack all data defining the regression problem in a problem structure.
problem.n = n;
problem.N = N;
problem.Nd = Nd;
problem.p = p;
problem.s = s;
problem.w = w;
problem.delta_tau = delta_tau;
problem.lambda = lambda;
problem.mu = mu;

%% Call the optimization procedure to compute the regression curve.
% Compute an initial guess for the curve. If this step is omitted, digress
% (below) will compute one itself. X0 is a 3D matrix of size n x n x Nd,
% such that each slice X0(:, :, k) is a rotation matrix.
%
X0 = initguess(problem);

%% my part
N1 = problem.N;
N2 = problem.Nd;
indices =  problem.s;
tau = problem.delta_tau;
miu = problem.mu;

Rdata = zeros(3,3*N1);
Rreg = zeros(3,3*N2);

% fill in data to Rdata
for i = 1:N1
    Rdata(:,i*3-2:i*3) = X0(:,:,indices(i));
end
for i = 1:N2
    Rreg(:,i*3-2:i*3) = X0(:,:,i);
end

[Rreg,newcosts] = non_optimization_on_so3(Rdata, Rreg, miu, lambda, indices, tau);
    
figure(7);
plot(newcosts,'r-o','LineWidth',2);
    

for i = 1:N2
    X1(:,:,i) = Rreg(:,i*3-2:i*3);
end

showSO3(Rdata,Rreg);

[speed0, acc0] = compute_profiles(problem, X0);
[speed1, acc1] = compute_profiles(problem, X1);

% Passage time of each point on the discrete curves.
time = problem.delta_tau*( 0 : (problem.Nd-1) );

figure(5);

subplot(1, 2, 1);
plot(1:N2,speed0,1:N2,speed1);
%     plot(time, speed0, time, speed1);
title('Speed of initial curve and optimized curve');
xlabel('Time');
ylabel('Speed');
legend('Initial curve', 'Optimized curve', 'Location', 'SouthEast');
pbaspect([1.6, 1, 1]);

subplot(1, 2, 2);
plot(1:N2,acc0,1:N2,acc1);
%     plot(time, acc0, time, acc1);
title('Acceleration of initial curve and optimized curve');
xlabel('Time');
ylabel('Acceleration');
legend('Initial curve', 'Optimized curve', 'Location', 'NorthWest');
pbaspect([1.6, 1, 1]);

ylim([0, 100]);

% Extend final sample to delay end of animation
Rdata = reshape(Rdata,3,3,[]);
Rreg = reshape(Rreg,3,3,[]);

for i = 1:size(Rreg,3)
    quat(i,:) = rot2quat(Rreg(:,:,i))';
    pos(i,:) = (i-1).*[0.01,0.01,0.01];
end

for i = 1:size(Rdata,3)
    quatd(i,:) = rot2quat(Rdata(:,:,i))';
end
posd = pos(indices,:);


posPlot = pos;
quatPlot = quat;

extraTime = 1;
samplePeriod = 1/100;
onesVector = ones(extraTime*(1/samplePeriod), 1);
posPlot = [posPlot; [posPlot(end, 1)*onesVector, posPlot(end, 2)*onesVector, posPlot(end, 3)*onesVector]];
quatPlot = [quatPlot; [quatPlot(end, 1)*onesVector, quatPlot(end, 2)*onesVector, quatPlot(end, 3)*onesVector, quatPlot(end, 4)*onesVector]];

% Create 6 DOF animation
SamplePlotFreq = 2;
Spin = 120;
filename = 'gait.mp4';
SixDOFanimation(posPlot, quat2rotm(quatPlot), posd, quat2rotm(quatd),...
                'SamplePlotFreq', SamplePlotFreq, 'Trail', 'All', ...
                'Position', [9 39 1280 768], 'View', [(100:(Spin/(length(posPlot)-1)):(100+Spin))', 10*ones(length(posPlot), 1)], ...
                'AxisLength', 0.1, 'ShowArrowHead', false, ...
                'Xlabel', 'X (m)', 'Ylabel', 'Y (m)', 'Zlabel', 'Z (m)', 'ShowLegend', false, ...
                'CreateAVI', false, 'AVIfileNameEnum', filename, 'AVIfps', ((1/samplePeriod) / SamplePlotFreq));


for i = 1:size(X0,3)
    quat0(i,:) = rot2quat(X0(:,:,i))';
end
ts = linspace(0,1,size(Rreg,3));
td = ts(indices);
t0 = ts;

figure
subplot(2,2,1);
plot(ts, quat(:,1),'r-','LineWidth',1.5);
hold on;grid on;
plot(td, quatd(:,1),'bs','MarkerSize',5);
plot(t0, quat0(:,1),'g--','LineWidth',1.5);
legend({'Regression','Anchor','Geodesics'},'FontName','Arial','FontSize',15);
title('q0','FontName','Arial','FontSize',15);

subplot(2,2,2);
plot(ts, quat(:,2),'r-','LineWidth',1.5);
hold on;grid on;
plot(td, quatd(:,2),'bs','MarkerSize',5);
plot(t0, quat0(:,2),'g--','LineWidth',1.5);
title('q1','FontName','Arial','FontSize',15);


subplot(2,2,3);
plot(ts, quat(:,3),'r-','LineWidth',1.5);
hold on;grid on;
plot(td, quatd(:,3),'bs','MarkerSize',5);
plot(t0, quat0(:,3),'g--','LineWidth',1.5);
title('q2','FontName','Arial','FontSize',15);


subplot(2,2,4);
plot(ts, quat(:,4),'r-','LineWidth',1.5);
hold on;grid on;
plot(td, quatd(:,4),'bs','MarkerSize',5);
plot(t0, quat0(:,4),'g--','LineWidth',1.5);
title('q3','FontName','Arial','FontSize',15);

warning on;