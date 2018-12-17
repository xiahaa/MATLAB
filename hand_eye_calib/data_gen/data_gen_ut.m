function varargout = data_gen_ut()
    addpath('../3rdparty/mit3dslam');
    addpath('../beautiful_plot');

    close all; clear all; clc;
    rand('seed', 1);

    %% from rot vec to quar, then to euler angles
    [r,p,w]=bot_quat_to_roll_pitch_yaw(bot_angle_axis_to_quat(pi/5,[.4 .5 .6]));
    trueCalib_euler = [0.6 0.3 0.4 r p w ];
    calib_dq = EulerStateToDualQuat(trueCalib_euler);
    %% transformation matrix from s1 to s2
    sensor2_expressedIn_sensor1 = GetHomoTransform(trueCalib_euler);

    motion = 3;
    N = 51; % Num poses
    N_interp = 300;% interpolation for surface
    if motion == 1
        N_modes = ceil(N/10);% perturbate modes
        rad1 = 10;%% radius
        rad2 = 5;
        amp_range = [0.0, 0.8];% amplitude
        freq_range = [0.5, 1.5];% frequency
        %% sim data, return values: position, derivative dz/dx, dz/dy, grid for surface
        [x, y, z,dz_dx,dz_dy, XX, YY, Z] = random_smooth_traj(N,rad1, rad2, N_interp, N_modes, amp_range, freq_range);
        xaxis = [diff(x); diff(y); diff(z)];%% always forward
    elseif motion == 2
        %% line motion
        len = 20; theta = deg2rad(25);
        [x, y, z, dz_dx, dz_dy, XX, YY, Z] = random_smooth_traj1(N, len, N_interp, theta);
    elseif motion == 3
        %% pure rotation
        x = zeros(1,N);
        y = zeros(1,N);
        z = zeros(1,N);
        X = linspace(min(x)-1.0, max(x)+1.0, N_interp);
        Y = linspace(min(y)-1.0, max(y)+1.0, N_interp);
        [XX,YY] = meshgrid(X,Y);
        Z = zeros(size(XX));
    elseif motion == 4
        N_modes = ceil(N/10);% perturbate modes
        rad1 = 10;%% radius
        amp_range = [0.0, 0.8];% amplitude
        freq_range = [0.5, 1.5];% frequency
       [x, y, z,dz_dx,dz_dy, XX, YY, Z] = random_smooth_traj2(N,rad1, N_interp, N_modes, ...
                                        amp_range, freq_range);
    elseif motion == 5
        %% line motion
        len = 20; theta = deg2rad(25);
        [x, y, z, dz_dx, dz_dy, XX, YY, Z] = random_smooth_traj1(N, len, N_interp, theta);
    end

    if motion == 1 || motion == 2 || motion == 4
        xaxis = [diff(x); diff(y); diff(z)];%% always forward
        yaxis = zeros(size(xaxis));
        zaxis = zeros(size(xaxis));
        normal = zeros(size(xaxis));
        for i = 1:length(x)-1
            normal(:,i) = cross( [1 0 dz_dx(i)], [0 1 dz_dy(i) ] );
            normal(:,i) = normal(:,i) / norm(normal(:,i));
            xaxis(:,i) = xaxis(:,i) /  norm(xaxis(:,i));
            yaxis(:,i) = cross ( normal(:,i), xaxis(:,i) );
            yaxis(:,i) = yaxis(:,i) /  norm(yaxis(:,i));
            zaxis(:,i) = cross ( xaxis(:,i), yaxis(:,i) );
            zaxis(:,i) = zaxis(:,i) /  norm(zaxis(:,i));
        end
    elseif motion == 3
        xaxis = [diff(x); diff(y); diff(z)];%% always forward
        yaxis = zeros(size(xaxis));
        zaxis = zeros(size(xaxis));
        normal = zeros(size(xaxis));
        for i = 1:length(x)-1
            u = rand(1,3);
            v = rand(1,3);
            normal(:,i) = cross( u, v);
            normal(:,i) = normal(:,i) / norm(normal(:,i));
            k = rand(1);
            xaxis(:,i) = k.*u + (1-k).*v;
            xaxis(:,i) = xaxis(:,i) / norm(xaxis(:,i));
            yaxis(:,i) = cross ( normal(:,i), xaxis(:,i) );
            yaxis(:,i) = yaxis(:,i) /  norm(yaxis(:,i));
            zaxis(:,i) = cross ( xaxis(:,i), yaxis(:,i) );
            zaxis(:,i) = zaxis(:,i) /  norm(zaxis(:,i));
        end
    else
        xaxis = [diff(x); diff(y); diff(z)];%% always forward
        yaxis = zeros(size(xaxis));
        zaxis = zeros(size(xaxis));
        normal = zeros(size(xaxis));
        scale = 0.1;
        for i = 1:length(x)-1
            normal(:,i) = cross( [1 0 dz_dx(i)], [0 1 dz_dy(i) ] );
            normal(:,i) = normal(:,i) / norm(normal(:,i));
            normal(:,i) = normal(:,i) + rand(3,1).*scale;normal(:,i) = normal(:,i) /  norm(normal(:,i));

            xaxis(:,i) = xaxis(:,i) /  norm(xaxis(:,i));
            xaxis(:,i) = xaxis(:,i) + rand(3,1).*scale;xaxis(:,i) = xaxis(:,i) /  norm(xaxis(:,i));

            yaxis(:,i) = cross ( normal(:,i), xaxis(:,i) );
            yaxis(:,i) = yaxis(:,i) /  norm(yaxis(:,i));
            zaxis(:,i) = cross ( xaxis(:,i), yaxis(:,i) );
            zaxis(:,i) = zaxis(:,i) /  norm(zaxis(:,i));

        end
    end

    sensor1_expressedIn_world = zeros(length(x),6);
    sensor2_expressedIn_world = zeros(length(x),6);
    T = eye(4);
    for i = 1:length(x)-1
        T(1:3,1:3) = [xaxis(:,i) yaxis(:,i) zaxis(:,i)];%% from body-1 to world
        T(1:3,4) = [x(i), y(i), z(i)];%%
        hold on;
        DrawAxis(T, 0.5, 'r', 'g', 'b');
        p1 = T(1:3,4);
        B = T * sensor2_expressedIn_sensor1;
        DrawAxis(B, 0.5, 'r', 'g', 'b');%% body-2 to world
        p2 = B(1:3,4);
        p = [p1 p2];
        disp(["baseline: ",num2str(norm(p1-p2))]);
        plot3(p(1,:),p(2,:),p(3,:),'k-','LineWidth',2);
        sensor1_expressedIn_world(i,:) = GetState(T);
        sensor2_expressedIn_world(i,:) = GetState(T * sensor2_expressedIn_sensor1);
    end
    sensor2_expressedIn_world(end,:) = sensor2_expressedIn_world(1,:);

    x2 = zeros(size(x));
    y2 = zeros(size(y));
    z2 = zeros(size(z));
    for idx = 1:length(x)
        x2(idx) = sensor2_expressedIn_world(idx,1);
        y2(idx) = sensor2_expressedIn_world(idx,2);
        z2(idx) = sensor2_expressedIn_world(idx,3);
    end

    font_size = 14;
    blue_color = [0 116 186]/255;
    orange_color = [223 80 35]/255;
    fig = gca;%figure();
    set(fig,'defaulttextinterpreter','latex');
    h1=plot3(x,y,z+0.1, 'LineWidth', 2.5, 'Color', blue_color);
    hold on;
    h2=plot3(x2,y2,z2+0.1, 'LineWidth', 2.5, 'Color', orange_color, 'LineStyle', '-.');
    s = surf(XX, YY, Z);
    lgnd = legend([h1,h2,s],{'Sensor $a$','Sensor $b$', 'Terrain'}, 'Location', 'NorthWest');
    set(lgnd, 'Interpreter', 'Latex','FontSize', font_size);
    colormap summer
    s.EdgeColor = 'none';
    s.FaceAlpha = 0.7;
    set(gca,'TickLabelInterpreter','latex');
    xlabel('$x$ (m)','FontSize', font_size, 'Interpreter', 'latex');
    ylabel('$y$ (m)','FontSize', font_size, 'Interpreter', 'latex');
    zlabel('$z$ (m)','FontSize', font_size, 'Interpreter', 'latex');
    hold off;
    axis equal;
    grid on;
    title('Simulation Data in 3D','Interpreter', 'latex','FontSize', font_size);
%     print(['./docs/figures/mit_ut_data_gen/', 'circle_terrain'],'-dpdf','-bestfit','-r300');

    sensor1_expressedIn_prevSensor1 = MakeRelState(sensor1_expressedIn_world);
    sensor2_expressedIn_prevSensor2 = TransformDiffHomo( sensor2_expressedIn_sensor1, sensor1_expressedIn_prevSensor1 );

    sensor2Initial_expressedIn_world = GetState(GetHomoTransform(sensor1_expressedIn_world(1,:)) * sensor2_expressedIn_sensor1);
    [~,sensor1check_expressedIn_world] = MakeAbsStates(sensor1_expressedIn_world(1,:), [], sensor1_expressedIn_prevSensor1, []);
    [~,sensor2check_expressedIn_world] = MakeAbsStates(sensor2Initial_expressedIn_world, [], sensor2_expressedIn_prevSensor2, []);

    font_size = 14;
    blue_color = [0 116 186]/255;
    red_color = [153 0 0]/255;
    figure(2);

    h4=plot3(sensor1check_expressedIn_world(1:end-1,1),sensor1check_expressedIn_world(1:end-1,2),sensor1check_expressedIn_world(1:end-1,3),'LineWidth', 2.5, 'Color', blue_color);hold on;
    h5=plot3(sensor2check_expressedIn_world(1:end-1,1),sensor2check_expressedIn_world(1:end-1,2),sensor2check_expressedIn_world(1:end-1,3),'LineWidth', 2.5, 'Color', red_color);
    axis equal;
    grid on;
    hold on;
    for i = 1:size(sensor1check_expressedIn_world,1)-1
        T1 = GetHomoTransform(sensor1check_expressedIn_world(i,:));
        DrawAxis(T1, 0.2, 'r', 'g', 'b');
        T2 = GetHomoTransform(sensor2check_expressedIn_world(i,:));
        DrawAxis(T2, 0.2, 'r', 'g', 'b');
        p1 = T1(1:3,4);
        p2 = T2(1:3,4);
        p = [p1 p2];
        disp(["baseline: ",num2str(norm(p1-p2))]);
        plot3(p(1,:),p(2,:),p(3,:),'k-','LineWidth',2);
    end
    hold off;
    lgnd2 = legend([h4,h5],{'Sensor $a$','Sensor $b$'}, 'Location', 'NorthWest');
    set(lgnd2, 'Interpreter', 'Latex','FontSize', font_size);
    colormap summer
    xlabel('$x$ (m)','FontSize', font_size, 'Interpreter', 'latex');
    ylabel('$y$ (m)','FontSize', font_size, 'Interpreter', 'latex');
    zlabel('$z$ (m)','FontSize', font_size, 'Interpreter', 'latex');
    set(gca,'TickLabelInterpreter','latex');

    if motion == 1
        title('Circle','Interpreter', 'latex','FontSize', font_size);
        print(['./docs/figures/mit_ut_data_gen/', 'circle'],'-dpdf','-bestfit','-r300');
        save('C:/Users/xiahaa/Documents/MATLAB/hand_eye_calib/data/circle.mat','sensor1check_expressedIn_world','sensor2check_expressedIn_world');
    elseif motion == 2
        title('Line','Interpreter', 'latex','FontSize', font_size);
        print(['./docs/figures/mit_ut_data_gen/', 'line'],'-dpdf','-bestfit','-r300');
        save('C:/Users/xiahaa/Documents/MATLAB/hand_eye_calib/data/line.mat','sensor1check_expressedIn_world','sensor2check_expressedIn_world');
    elseif motion == 3
        title('Pure Rotation','Interpreter', 'latex','FontSize', font_size);
        print(['./docs/figures/mit_ut_data_gen/', 'rotation'],'-dpdf','-bestfit','-r300');
        save('C:/Users/xiahaa/Documents/MATLAB/hand_eye_calib/data/rotation.mat','sensor1check_expressedIn_world','sensor2check_expressedIn_world');
    elseif motion == 4
        title('$8-Shape$','Interpreter', 'latex','FontSize', font_size);
        print(['./docs/figures/mit_ut_data_gen/', 'shape8'],'-dpdf','-bestfit','-r300');
        save('C:/Users/xiahaa/Documents/MATLAB/hand_eye_calib/data/shape8.mat','sensor1check_expressedIn_world','sensor2check_expressedIn_world');
    elseif motion == 5
        title('$Small Rotation$','Interpreter', 'latex','FontSize', font_size);
        print(['./docs/figures/mit_ut_data_gen/', 'smallr'],'-dpdf','-bestfit','-r300');
        save('C:/Users/xiahaa/Documents/MATLAB/hand_eye_calib/data/samllr.mat','sensor1check_expressedIn_world','sensor2check_expressedIn_world');
    end

    % generate random covariances
    %RandStream.setDefaultStream(RandStream('mt19937ar','seed',0)); % set the random seed so we always get the same random values
    std = 0.55;%0.25;%0.5 0.75 1
    tstd = std*ones(1,3);
    tstr = std*ones(1,3);

    std_sensor1_per_unit = [tstd tstr];
    std_sensor2_per_unit = [tstd tstr];

    cov1 = zeros(size(sensor1_expressedIn_prevSensor1,1), 36);
    cov2 = zeros(size(sensor1_expressedIn_prevSensor1,1), 36);
    for i = 1:size(sensor1_expressedIn_prevSensor1,1)
        motion1 = abs(sensor1_expressedIn_prevSensor1(i,:));
        motion2 = abs(sensor2_expressedIn_prevSensor2(i,:));
        S1 = diag((motion1 .* std_sensor1_per_unit).^2+(1e-20)^2);% + randn(6)*(1e-3)^2;
        S2 = diag((motion2 .* std_sensor2_per_unit).^2+(1e-20)^2);% + randn(6)*(1e-3)^2;
        S1 = triu(S1)+triu(S1,1)';
        S2 = triu(S2)+triu(S2,1)';
        cov1(i,:) = reshape(S1,1,[]);
        cov2(i,:) = reshape(S2,1,[]);
    end

    numberOfRuns = 100;
    run.true.sensor1_expressedIn_prevSensor1 = sensor1_expressedIn_prevSensor1;
    run.true.sensor2_expressedIn_prevSensor2 = sensor2_expressedIn_prevSensor2;
    run.true.sensor1_expressedIn_world = sensor1_expressedIn_world;
    run.true.sensor2_expressedIn_world = sensor2_expressedIn_world;
    run.true.cov1 = cov1;
    run.true.cov2 = cov2;
    run.true.calib = trueCalib_euler;
    for i = 1:numberOfRuns
        run.observations{i}.sensor1_expressedIn_prevSensor1 = SampleVelocitiesWithCovariance(run.true.sensor1_expressedIn_prevSensor1, run.true.cov1);
        run.observations{i}.sensor2_expressedIn_prevSensor2 = SampleVelocitiesWithCovariance(run.true.sensor2_expressedIn_prevSensor2, run.true.cov2);
    end

    if motion == 1
        save(['C:/Users/xiahaa/Documents/MATLAB/hand_eye_calib/data/circle100_',num2str(std),'.mat'], 'run');
    elseif motion == 2
        save(['C:/Users/xiahaa/Documents/MATLAB/hand_eye_calib/data/line100_',num2str(std),'.mat'],'run');
    elseif motion == 3
        save(['C:/Users/xiahaa/Documents/MATLAB/hand_eye_calib/data/rotation100_',num2str(std),'.mat'], 'run');
    elseif motion == 4
        save(['C:/Users/xiahaa/Documents/MATLAB/hand_eye_calib/data/shape8100_',num2str(std),'.mat'], 'run');
    elseif motion == 5
        save(['C:/Users/xiahaa/Documents/MATLAB/hand_eye_calib/data/smallr100_',num2str(std),'.mat'], 'run');
    end

end

function observed_pos_expressedIn_prevPos = SampleVelocitiesWithCovariance(true_pos_expressedIn_prevPos, true_pos_expressedIn_prevPos_cov)
    observed_pos_expressedIn_prevPos = zeros(size(true_pos_expressedIn_prevPos));
    for i = 1:size(observed_pos_expressedIn_prevPos,1)
        c = reshape(true_pos_expressedIn_prevPos_cov(i,:),6,6);
        observed_pos_expressedIn_prevPos(i,:) = mgd(1, 6, true_pos_expressedIn_prevPos(i,:), c);
    end
end

function [x, y, z, dz_dx, dz_dy, XX, YY, Z] = random_smooth_traj1(N, len, N_interp, theta)
%random_smooth_traj

xg = linspace(0,len,N);
yg = zeros(1,N);

pg = [xg;yg];

x = [cos(theta) sin(theta)]*pg;
y = [-sin(theta) cos(theta)]*pg;

X = linspace(min(x)-1.0, max(x)+1.0, N_interp);
Y = linspace(min(y)-1.0, max(y)+1.0, N_interp);
[XX,YY] = meshgrid(X,Y);
Z = zeros(size(XX));

dz_dx = zeros(size(x));
dz_dy = zeros(size(y));

z = interp2(XX,YY,Z,x,y);

end

function [x, y, z, dz_dx, dz_dy, XX, YY, Z] = random_smooth_traj(N,rad1, rad2, N_interp, N_modes, ...
                                        amp_range, freq_range)
%random_smooth_traj

theta = linspace(0, 2*pi, N);
x = cos(theta)*rad1;
y = sin(theta)*rad2;

X = linspace(min(x)-1.0, max(x)+1.0, N_interp);
Y = linspace(min(y)-1.0, max(y)+1.0, N_interp);
[XX,YY] = meshgrid(X,Y);
Z = zeros(size(XX));

dz_dx = zeros(size(x));
dz_dy = zeros(size(y));

for idx=1:N_modes
    phase = 2*rand*pi -pi;
    amp = rand*(amp_range(2) - amp_range(1)) + amp_range(1);
    freq = rand*(freq_range(2) - freq_range(1)) + freq_range(1);
    mix = rand;
    Z = Z + amp*cos((mix*XX + (1-mix)*YY)*freq + phase);

    dz_dx = dz_dx - amp*sin((mix*x + (1-mix)*y)*freq + phase)*freq*mix;
    dz_dy = dz_dy - amp*sin((mix*x + (1-mix)*y)*freq + phase)*(1-mix)*freq;

    % Do it again for Y
%     phase = 2*rand*pi -pi;
%     amp = rand*(amp_range(2) - amp_range(1)) + amp_range(1);
%     freq = rand*(freq_range(2) - freq_range(1)) + freq_range(1);
%     Z = Z + amp*cos(YY*freq + phase);

end

z = interp2(XX,YY,Z,x,y);

% figure;
% surf(XX,YY,Z);
end


function [x, y, z, dz_dx, dz_dy, XX, YY, Z] = random_smooth_traj2(N,rad1, N_interp, N_modes, ...
                                        amp_range, freq_range)
%random_smooth_traj
hN1 = round(N/2);
hN2 = N - hN1;
theta1 = linspace(0, pi, hN1);
theta2 = linspace(0, -pi, hN2);
theta = [theta1 theta2];

hN3 = round(hN1*0.5);
hN4 = hN1 - hN3;
hN5 = round(hN2*0.5);
hN6 = hN2 - hN5;

rho = [linspace(0.1, rad1, hN3) linspace(rad1, 0.1, hN4) linspace(0.1, rad1, hN5) linspace(rad1, 0.1, hN6)];
x = cos(theta).*rho;
y = sin(theta).*rho;

X = linspace(min(x)-1.0, max(x)+1.0, N_interp);
Y = linspace(min(y)-1.0, max(y)+1.0, N_interp);
[XX,YY] = meshgrid(X,Y);
Z = zeros(size(XX));

dz_dx = zeros(size(x));
dz_dy = zeros(size(y));

for idx=1:N_modes
    phase = 2*rand*pi -pi;
    amp = rand*(amp_range(2) - amp_range(1)) + amp_range(1);
    freq = rand*(freq_range(2) - freq_range(1)) + freq_range(1);
    mix = rand;
    Z = Z + amp*cos((mix*XX + (1-mix)*YY)*freq + phase);

    dz_dx = dz_dx - amp*sin((mix*x + (1-mix)*y)*freq + phase)*freq*mix;
    dz_dy = dz_dy - amp*sin((mix*x + (1-mix)*y)*freq + phase)*(1-mix)*freq;

    % Do it again for Y
%     phase = 2*rand*pi -pi;
%     amp = rand*(amp_range(2) - amp_range(1)) + amp_range(1);
%     freq = rand*(freq_range(2) - freq_range(1)) + freq_range(1);
%     Z = Z + amp*cos(YY*freq + phase);

end

z = interp2(XX,YY,Z,x,y);

% figure;
% surf(XX,YY,Z);
end