% run as series of hand eye calibration methods on the RGB dataset from MIT
% and save the results.
clc;
close all;
clear all;

addpath('../beautiful_plot/');
addpath('./solver/');
addpath('../MatrixLieGroup/barfoot_tro14');
addpath('../quaternion');
addpath('./solver/atadq');

%% dataset path
addpath('D:/dtu/sourcecode/hand_eye/3dcalib/3d-experiment/');

id = 1;

load(strcat('prime_', num2str(id),'.mat'));

eye = dat.dq1;
hand = dat.dq2;

if (size(eye,1) ~= size(hand,1))
    error('incomplete data: size not equal!');
end

truth = [0.000000000000000  -0.000000000000000  -0.034899496702501   0.999390827019096   0.000436243708781   0.012492385337739  0.019987816540382   0.000697989934050];
dq = truth(1:4)';
dqd = truth(5:8)';
R = q2R(dq);
t = 2.*qprod(dqd, conjugateq(dq));
t = t(2:4);
Tt = [R t;[0 0 0 1]];

estMIT = [0.003692920499665  -0.008068642457500  -0.032828870096415   0.999421595041486   0.002952433853242   0.008160798751105   0.017042622011385   0.000614789084959];
dq = estMIT(1:4)';
dqd = estMIT(5:8)';
R = q2R(dq);
t = 2.*qprod(dqd, conjugateq(dq));
t = t(2:4);
Tmit = [R t;[0 0 0 1]];

num = size(eye,1);
T1 = zeros(num,4,4);
T2 = zeros(num,4,4);

for i = 1:num
    dq = eye(i,1:4)';
    dqd = eye(i,5:8)';
    R = q2R(dq);
    t = 2.*qprod(dqd, conjugateq(dq));
    t = t(2:4);
    T1(i,1:4,1:4) = [R t;[0 0 0 1]];
    
    dq = hand(i,1:4)';
    dqd = hand(i,5:8)';
    R = q2R(dq);
    t = 2.*qprod(dqd, conjugateq(dq));
    t = t(2:4);
    T2(i,1:4,1:4) = [R t;[0 0 0 1]];
end

%     @sol_horaud, ...                                    %% QSEP
%     @sol_horaud_nlopt, ...                     %% NLOPT
%     @sol_andreff, ...                                   %% KR


convSolver = { ...                       %% TSAI
    @sol_park_martin, ...                               %% LIE
    @sol_dual_quaternion, ...                           %% DQ
    @sol_chou, ...                                      %% IDQ
    @sol_adjoint_transformation_algo, ...
    @sol_dphec, ...                            %% GLOBAL_POLY
    @sol_dual_sdp_cvx, ...                     %% DUAL_SDP
    @sol_cvx2, ...                             %% SCF
    @sol_manifold_opt_SE3, ...                    %% SE3OPT
    @sol_cvx_sdp, ...
    };
 
 convSolver1 = {@sol_tsai_lenz, ...                       %% TSAI
    @sol_park_martin, ...                               %% LIE
    @sol_horaud, ...                                    %% QSEP
    @sol_andreff, ...                                   %% KR
    @HandEye_DQ, ...                                    %% DQ
    @sol_chou, ...                                      %% CHOU
    @sol_improved_dual_quaternion, ...                   %% IDQ
    };

advSolver = {@sol_adjoint_transformation_algo, ...      %% ATA
             @sol_horaud_nlopt, ...                     %% NLOPT
             @sol_cvx1, ...                             %% SOCP
             @sol_dphec, ...                            %% GLOBAL_POLY
             @sol_dual_sdp_cvx, ...                     %% DUAL_SDP
             @sol_cvx2, ...                             %% SCF
             @sol_cvx_sdp};                    %% SE3OPT

     convSolver1 = {
        @sol_andreff, ...                                   %% KR
        @sol_cvx1, ...                                      %% SOCP
        @sol_adjoint_transformation_algo, ...               %% ATA
        @sol_dphec, ...                                     %% GLOBAL_POLY
        @sol_dual_sdp_cvx, ...                              %% DUAL_SDP
        @sol_cvx_sdp, ...                                   %% SDP
        };
    
    convSolver1 = {
        @sol_andreff, ...                                   %% KR
        @sol_horaud_nlopt, ...                              %% DUAL_SDP
        @sol_cvx1, ...                                      %% SOCP
        @sol_adjoint_transformation_algo, ...               %% ATA
        @sol_dphec, ...                                     %% GLOBAL_POLY
        @sol_manifold_opt_SE3, ...                                   %% SDP
        };
         
 usedsolver = convSolver1;
     

for kk = 1:size(usedsolver, 2)
    handle_sol = usedsolver{kk};
    tic
    TX = handle_sol(T1,T2,num);
    tsol(kk) = toc;
    if isempty(TX) == false
        xsol(1:4,1:4,kk) = TX;
        flag(kk) = 1;
    else
        flag(kk) = 0;
    end
end
ts = [3.797894 tsol];
%% compared with truth
err_with_truth = zeros(size(usedsolver, 2)+1,3);
[err1,err2,err3] = errors2(Tmit, Tt);
err_with_truth(1,:) = [err1,err2,err3];
for kk = 1:size(usedsolver, 2)
    X = xsol(1:4,1:4,kk);
    [err1,err2,err3] = errors2(X, Tt);
    err_with_truth(kk+1,:) = [err1,err2,err3];
end
% 'KR',
% convSols = {'MIT', 'TSAI', 'LIE', 'QSEP', 'KR', 'DQ', 'CHOU', 'IDQ'};
% convSols = {'MIT', 'LIE',  'DQ', 'CHOU', 'ATA', 'GPOLY', 'DUAL', 'SCF', 'SE3OPT', 'SDP'};
% convSols = {'BL','KR','SOCP','ATA','GPOLY','DUAL','SDR'};
    convSols = {'BL','KR','NLQ','SOCP','ATA','GPOLY','SE3'};


red_color = [153 0 0]/255;
red_color = [red_color 0.5];
font_size = 15;

fig = figure();
set(fig,'defaulttextinterpreter','latex');
box_labels = categorical(convSols);
box_labels = reordercats(box_labels,convSols);
subplot(3,1,1);plot(box_labels, err_with_truth(:,1)', '-o', 'LineWidth', 2.5, 'Color', red_color);grid on;ylabel('$E_{\mathbf{R}}$','Interpreter','latex');
title('$RGB-D\ Camera\ Dataset$','FontSize', font_size,'Interpreter','latex');set(gca,'FontSize', font_size, 'TickLabelInterpreter','latex');

subplot(3,1,2);plot(box_labels, err_with_truth(:,2)', '-o', 'LineWidth', 2.5, 'Color', red_color);grid on;ylabel('$E_{\mathbf{t}}$','Interpreter','latex');
set(gca,'FontSize', font_size, 'TickLabelInterpreter','latex');

subplot(3,1,3);plot(box_labels, ts, '-d', 'LineWidth', 2.5, 'Color', red_color);grid on;ylabel('$Time:\ (s)$','Interpreter','latex');
set(gca,'FontSize', font_size, 'TickLabelInterpreter','latex');
% print(['./docs/figures/eth_asl/advResult_', 'MIT_Dataset'],'-dpdf','-bestfit','-r300');

fig = figure();
set(fig,'defaulttextinterpreter','latex');
box_labels = categorical(convSols);
box_labels = reordercats(box_labels,convSols);
subplot(3,1,1);plot(box_labels, err_with_truth(:,1)', '-o', 'LineWidth', 3, 'Color', red_color);grid on;ylabel('$E_{\mathbf{R}}$','Interpreter','latex');
title('$RGB-D\ Camera\ Dataset$','FontSize', font_size,'Interpreter','latex');set(gca,'FontSize', font_size, 'TickLabelInterpreter','latex');

subplot(3,1,2);plot(box_labels, err_with_truth(:,2)', '-o', 'LineWidth', 3, 'Color', red_color);grid on;ylabel('$E_{\mathbf{t}}$','Interpreter','latex');
set(gca,'FontSize', font_size, 'TickLabelInterpreter','latex');

subplot(3,1,3);plot(box_labels, err_with_truth(:,3)', '-o', 'LineWidth', 3, 'Color', red_color);grid on;ylabel('$E_{\mathbf{T}}$','Interpreter','latex');
set(gca,'FontSize', font_size, 'TickLabelInterpreter','latex');
% print(['./docs/figures/eth_asl/advResult_', 'MIT_Dataset'],'-dpdf','-bestfit','-r300');

%% compute RMSE
errcs1 = zeros(size(usedsolver, 2),num);
errcs2 = zeros(size(usedsolver, 2),num);
errcs3 = zeros(size(usedsolver, 2),num);

for kk = 1:size(usedsolver, 2)
    X = xsol(1:4,1:4,kk);
    errs = zeros(num,3);
    for i = 1:num
        A = T2(i,:,:);A = reshape(A,4,4,1);
        B = T1(i,:,:);B = reshape(B,4,4,1);
        [err1,err2,err3] = errors(A,B,X);
        errs(i,:) = [err1,err2,err3];
    end
    errcs1(kk,:) = errs(:,1)';
    errcs2(kk,:) = errs(:,2)';
    errcs3(kk,:) = errs(:,3)';
end
convSols = {'KR','SOCP','ATA','GPOLY','DUAL','SDR'};

RMSE_R = mean(errcs1,2);
RMSE_t = mean(errcs2,2);
RMSE_T = mean(errcs3,2);
box_labels = categorical(convSols);
box_labels = reordercats(box_labels,convSols);
fig = figure();
set(fig,'defaulttextinterpreter','latex');
subplot(3,1,1);plot(box_labels, RMSE_R');
subplot(3,1,2);plot(box_labels, RMSE_t');
subplot(3,1,3);plot(box_labels, RMSE_T');

%% bar plot of error
font_size = 16;
line_width = 2.5;
subplot_margin = 0.12;
subplot_spacing = 0.1;
%% compare different noise level
fig = figure();
set(fig,'defaulttextinterpreter','latex');
box_labels = categorical(convSols);
box_labels = reordercats(box_labels,convSols);
boxplot(errcs1',box_labels);grid on;
title('Rotational Error')
ylabel('$E_{R}$','Interpreter','latex');

fig = figure();
set(fig,'defaulttextinterpreter','latex');
box_labels = categorical(convSols);
box_labels = reordercats(box_labels,convSols);
boxplot(errcs2',box_labels);grid on;
title('Translation Error')
ylabel('$E_{t}$','Interpreter','latex');

fig = figure();
set(fig,'defaulttextinterpreter','latex');
box_labels = categorical(convSols);
box_labels = reordercats(box_labels,convSols);
boxplot(errcs3',box_labels);grid on;
title('Total Error')
ylabel('$E_{T}$','Interpreter','latex');

function varargout = errors(A,B,X)
    R = A(1:3,1:3)*X(1:3,1:3) - X(1:3,1:3)*B(1:3,1:3);
    err1 = norm(R,'fro');
    t = A(1:3,1:3) * X(1:3,4) + A(1:3,4) - X(1:3,1:3) * B(1:3,4) - X(1:3,4);
    err2 = norm(t,2);
    T = A*X-X*B;
    err3 = norm(T,'fro');
    varargout{1} = err1;
    varargout{2} = err2;
    varargout{3} = err3;
end

function varargout = errors2(T, Tt)
    R = Tt(1:3,1:3)'*T(1:3,1:3) - eye(3);
    err1 = norm(R,'fro');
    t = Tt(1:3,4) - T(1:3,4);
    err2 = norm(t,2);
    Te = inv(Tt)*T - eye(4);
    err3 = norm(Te,'fro');
    varargout{1} = err1;
    varargout{2} = err2;
    varargout{3} = err3;
end

