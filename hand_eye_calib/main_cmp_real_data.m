% run a series of hand eye calibration methods on the ETH-ASL datasets and
% save the result.
clc;
close all;
clear all;

addpath('../beautiful_plot/');
addpath('./solver/');
addpath('../MatrixLieGroup/barfoot_tro14');
addpath('../quaternion');
addpath('./solver/atadq');

%% dataset path
addpath('C:/Users/xiahaa/Documents/MATLAB/hand_eye_calib/data/eth_asl_real/');

id = 2;

load(strcat('prime_', num2str(id),'.mat'));

eye = dat.eye;
hand = dat.hand;

idnan1 = isnan(eye(:,2)) | isnan(eye(:,3)) | isnan(eye(:,4)) | isnan(eye(:,5)) | isnan(eye(:,6)) | isnan(eye(:,7)) | isnan(eye(:,8));
idnan2 = isnan(hand(:,2)) | isnan(hand(:,3)) | isnan(hand(:,4)) | isnan(hand(:,5)) | isnan(hand(:,6)) | isnan(hand(:,7)) | isnan(hand(:,8));
idfilter = idnan1 | idnan2;
eye = eye(~idfilter,:);
hand = hand(~idfilter,:);

if (size(eye,1) ~= size(hand,1))
    error('incomplete data: size not equal!');
end

num = size(eye,1);
nsam = num - 1;
T1 = zeros(nsam,4,4);
T2 = zeros(nsam,4,4);

% num = 10;
% dR = euler2rot(rand(1)*pi-pi*0.5,rand(1)*pi-pi*0.5,rand(1)*2*pi-pi)
% dt = rand(3,1)*5-2.5
% dT = [dR dt;[0 0 0 1]];
% dTt = [dR' -dR'*dt;[0 0 0 1]];
% for i = 1:num
%     R1 = euler2rot(rand(1)*pi-pi*0.5,rand(1)*pi-pi*0.5,rand(1)*2*pi-pi);
%     t1 = rand(3,1)*20-10;
%     TSS = [R1 t1;[0 0 0 1]];
%     TB = TSS;
%     TA = dT * TSS * dTt;
%     a = rot2quat(TB(1:3,1:3))';
%     hand(i,2:end) = [TB(1:3,4)' [a(2) a(3) a(4) a(1)]];
%     a = rot2quat(TA(1:3,1:3))';
%     eye(i,2:end) = [TA(1:3,4)' a(2) a(3) a(4) a(1)];
% end

t = eye(1,2:4);
q = [eye(1,8);eye(1,5:7)'];
R = q2R(q);
Te0 = [R t';[0 0 0 1]];
t = hand(1,2:4);
q = [hand(1,8);hand(1,5:7)'];
R = q2R(q);
Th0 = [R t';[0 0 0 1]];

for i = 2:num
    t = eye(i,2:4);
    q = [eye(i,8);eye(i,5:7)'];
    R = q2R(q);
    T1(i-1,1:4,1:4) = [R t';[0 0 0 1]] * inv(Te0);
    t = hand(i,2:4);
    q = [hand(i,8);hand(i,5:7)'];
    R = q2R(q);
    T2(i-1,1:4,1:4) = [R t';[0 0 0 1]] * inv(Th0);
end



convSolver = {@sol_tsai_lenz, ...                       %% TSAI
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
             @sol_cvx_sdp};                             %% SE3OPT

    convSolver1 = {
        @sol_andreff, ...                                   %% KR
        @sol_horaud_nlopt, ...                              %% DUAL_SDP
        @sol_adjoint_transformation_algo, ...               %% ATA
        @sol_cvx1, ...                                      %% SOCP
        @sol_dphec, ...                                     %% GLOBAL_POLY
        @sol_dual_sdp_cvx, ...
        @batchSolveSoftUseScrew, ...
        };
%       @sol_cvx1, ...                                      %% SOCP
%       @sol_manifold_opt_SE3, ...                                   %% SDP

    solver_name = {'KR','NLQ','ATA','SOCP', 'GPOLY', 'DUAL','BS'};%'SOCP',
%     solver_name = {'KR','NLQ','SOCP','ATA','GPOLY','SE3'};
         
         
 usedsolver = convSolver1;

for kk = 1:size(usedsolver, 2)
    handle_sol = usedsolver{kk};
    if strcmp(solver_name(kk),'SDP')
        tic
        TX = handle_sol(T1,T2,nsam, Th0(1:3,1));
        tsol(kk) = toc;
    elseif strcmp(solver_name(kk),'BS')
        A = zeros(4,4,nsam);
        B = zeros(4,4,nsam);
        for i = 1:nsam
            Ts = T2(i,:,:);Ts = reshape(Ts,4,4);
            A(:,:,i) = Ts;
            Ts = T1(i,:,:);Ts = reshape(Ts,4,4);
            B(:,:,i) = Ts;
        end 
        tic
        TX = handle_sol(A,B);
        tsol(kk) = toc;
    else
        tic
        TX = handle_sol(T1,T2,nsam);
        tsol(kk) = toc;
    end
    if isempty(TX) == false
        xsol(1:4,1:4,kk) = TX;
        flag(kk) = 1;
    else
        flag(kk) = 0;
    end
end

%% compute RMSE
errcs1 = zeros(size(usedsolver, 2),nsam);
errcs2 = zeros(size(usedsolver, 2),nsam);
errcs3 = zeros(size(usedsolver, 2),nsam);

for kk = 1:size(usedsolver, 2)    
    X = xsol(1:4,1:4,kk);
    errs = zeros(nsam,3);
    for i = 1:nsam
        A = T2(i,:,:);A = reshape(A,4,4,1);
        B = T1(i,:,:);B = reshape(B,4,4,1);
        [err1,err2,err3] = errors(A,B,X);
        errs(i,:) = [err1,err2,err3];
    end
    errcs1(kk,:) = errs(:,1)';
    errcs2(kk,:) = errs(:,2)';
    errcs3(kk,:) = errs(:,3)';
end
% convSols = {'TSAI', 'LIE', 'QSEP', 'KR', 'DQ', 'CHOU', 'IDQ'};
convSols = solver_name;

red_color = [153 0 0]/255;
font_size = 12;

RMSE_R = mean(errcs1,2);
RMSE_t = mean(errcs2,2);
RMSE_T = mean(errcs3,2);
box_labels = categorical(convSols);
box_labels = reordercats(box_labels,convSols);
fig = figure();
set(fig,'defaulttextinterpreter','latex');
subplot(3,1,1);plot(box_labels, RMSE_R', '-o', 'LineWidth', 3, 'Color', red_color);grid on;ylabel('$E_{\mathbf{R}}$','Interpreter','latex');
title(['Dataset: ', 'Prime_', num2str(id)]);set(gca,'FontSize', font_size, 'TickLabelInterpreter','latex');
subplot(3,1,2);plot(box_labels, RMSE_t', '-o', 'LineWidth', 3, 'Color', red_color);grid on;ylabel('$E_{\mathbf{t}}$','Interpreter','latex');
set(gca,'FontSize', font_size, 'TickLabelInterpreter','latex');
subplot(3,1,3);plot(box_labels, RMSE_T', '-o', 'LineWidth', 3, 'Color', red_color);grid on;ylabel('$E_{\mathbf{T}}$','Interpreter','latex');
set(gca,'FontSize', font_size, 'TickLabelInterpreter','latex');
% print(['./docs/figures/eth_asl/adv_Result_', 'Prime_', num2str(id)],'-dpdf','-bestfit','-r300');

save(strcat('./data/SE3/','prime_', num2str(id),'_res','.mat'),'errcs1','errcs2','errcs3','convSols','tsol');

%% bar plot of error
fontsize = 12;
xlabel('$Standard\ deviation\ of\ added\ noise$','Interpreter','latex','FontSize',fontsize);
ylabel('$E_{\mathbf{R}}$','Interpreter','latex','FontSize',fontsize);
%     legend(solver_name,'Interpreter','latex','FontSize',8,'Location', 'northeast');
title('Rotation Error\ Comparison','Interpreter','latex','FontSize',fontsize);
ylim([0,2]);

font_size = 12;
%% compare different noise level
fig = figure();
set(fig,'defaulttextinterpreter','latex');
box_labels = categorical(convSols);
box_labels = reordercats(box_labels,convSols);
boxplot(errcs1',box_labels);grid on;
title('Rotational Error')
ylabel('$E_{R}$','Interpreter','latex');
set(gca,'FontSize', font_size, 'TickLabelInterpreter','latex');

fig = figure();
set(fig,'defaulttextinterpreter','latex');
box_labels = categorical(convSols);
box_labels = reordercats(box_labels,convSols);
boxplot(errcs2',box_labels);grid on;
title('Translation Error')
ylabel('$E_{t}$','Interpreter','latex');
set(gca,'FontSize', font_size, 'TickLabelInterpreter','latex');

fig = figure();
set(fig,'defaulttextinterpreter','latex');
box_labels = categorical(convSols);
box_labels = reordercats(box_labels,convSols);
boxplot(errcs3',box_labels);grid on;
title('Total Error')
ylabel('$E_{T}$','Interpreter','latex');
set(gca,'FontSize', font_size, 'TickLabelInterpreter','latex');

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

