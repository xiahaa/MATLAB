function main_sim_cmp_runtime_ut_mit
% main entry for comparing the runtimes for a series of hand eye
% calibration methods by varing the number of available motion samples.
    clc;clear all;close all;
    addpath('../3rdparty/mit3dslam');
    addpath('../beautiful_plot');
    addpath('./solver/');
    addpath('../MatrixLieGroup/barfoot_tro14');
    addpath('../quaternion');
    addpath('./solver/atadq');

    id = 1;
    prefix = 'C:/Users/xiahaa/Documents/MATLAB/hand_eye_calib/data/';
    name = {'circle100','line100','rotation100','shape8100','smallr100'};
    suffix = '.mat';
%     stds = [0 0.1 0.25 0.5 0.75 1];
    stds = [0.05 0.15 0.25 0.35 0.45 0.55];
    
    truth = [0.9511    0.1409    0.1761    0.2113   -0.1109    0.2818    0.1074    0.2219];
    dq = truth(1:4)';
    dqd = truth(5:8)';
    R = q2R(dq);
    t = 2.*qprod(dqd, conjugateq(dq));
    t = t(2:4);
    Tt = [R t;[0 0 0 1]];
    
    
%     convSolver = {@sol_tsai_lenz, ...                       %% TSAI
%         @sol_park_martin, ...                               %% LIE
%         @sol_horaud, ...                                    %% QSEP
%         @sol_andreff, ...                                   %% KR
%         @sol_dual_quaternion, ...                           %% DQ
%         @sol_chou, ...                                      %% IDQ
%         @sol_adjoint_transformation_algo, ...               %% ATA
%         @sol_dphec, ...                                     %% GLOBAL_POLY
%         @sol_dual_sdp_cvx, ...                              %% DUAL_SDP
%         @sol_cvx2, ...                                      %% SCF
%         @sol_manifold_opt_SE3, ...                          %% SE3OPT
%         @sol_cvx_sdp, ...                                   %% SDP
%         @sol_horaud_nlopt, ...                              %% NLOPT
%         @sol_cvx1, ...                                      %% SOCP
%         };
    %     @HandEye_DQ, ...                                    %% DQ
    %     @sol_improved_dual_quaternion, ...                   %% IDQ
         
%     solver_name = {'TS','LIE','QS','KR','DQ','CHOU','ATA','GPOLY','DUAL','SCF','SE3','SDR','QNL','SOCP'};
    
%     convSolver = {
%         @sol_andreff, ...                                   %% KR
%         @sol_adjoint_transformation_algo, ...               %% ATA
%         @sol_dphec, ...                                     %% GLOBAL_POLY
%         @sol_dual_sdp_cvx, ...                              %% DUAL_SDP
%         @sol_cvx_sdp, ...                                   %% SDP
%         @sol_cvx1, ...                                      %% SOCP
%         sol_manifold_opt_SE3, ...                           %% SE3
%         };
% 
%     solver_name = {'KR','ATA','GPOLY','DUAL','SDR','SOCP','SE3'};
%   mapping_id = [1,2,6,3,4,5];

    convSolver = {
        @sol_andreff, ...                                   %% KR
        @sol_horaud_nlopt, ...                              %% NLOPT
        @sol_cvx1, ...                                      %% SOCP
        @sol_adjoint_transformation_algo, ...               %% ATA
        @sol_dphec, ...                                     %% GLOBAL_POLY
        @sol_manifold_opt_SE3, ...                           %% SE3
        @batchSolveSoftUseScrew, ...
        };

%     solver_name = {'KR','SOCP','ATA','GPOLY','DUAL','SDR'};%,'DUAL','SDR'
    solver_name = {'KR','NLQ','SOCP','ATA','GPOLY','SE3','batchnew'};%,'DUAL','SDR'
%     mapping_id = [1,6,4,3,4,5];
    mapping_id = [1,2,6,4,3,6,7];
    
    usedsolver = convSolver;
    numSam = [20 40 60 80 100 120 140 160 180];
    
    std = stds(2);filename = strcat(prefix,name{id},'_',num2str(std),suffix);load(filename);
    times = zeros(numel(numSam), size(solver_name,2), 100);
    valid_id = ones(numel(numSam),size(solver_name,2),100);
    if 0
        for k = 1:numel(numSam)
            ns = numSam(k);
            randid = randperm(200,ns);
            for j = 1:100
                sensor1_expressedIn_prevSensor1 = run.observations{j}.sensor1_expressedIn_prevSensor1;
                sensor2_expressedIn_prevSensor2 = run.observations{j}.sensor2_expressedIn_prevSensor2;
                hand = EulerStateToDualQuat(sensor1_expressedIn_prevSensor1);
                eye = EulerStateToDualQuat(sensor2_expressedIn_prevSensor2);

                num = numel(randid);
                T1 = zeros(num,4,4);
                T2 = zeros(num,4,4);
                for kk = 1:numel(randid)
                    i = randid(kk);
                    dq = eye(i,1:4)';
                    dqd = eye(i,5:8)';
                    R = q2R(dq);
                    t = 2.*qprod(dqd, conjugateq(dq));
                    t = t(2:4);
                    T1(kk,1:4,1:4) = [R t;[0 0 0 1]];

                    dq = hand(i,1:4)';
                    dqd = hand(i,5:8)';
                    R = q2R(dq);
                    t = 2.*qprod(dqd, conjugateq(dq));
                    t = t(2:4);
                    T2(kk,1:4,1:4) = [R t;[0 0 0 1]];
                end

                solver_id = 6;
    %             for solver_id = 1:size(solver_name,2)
                    handle_sol = usedsolver{solver_id};
                    disp([solver_id k,j])
                    tic
                    TX = handle_sol(T1,T2,num);
                    time = toc;
                    if isempty(TX)
                        valid_id(k,j) = 0;
                        continue;
                    end
                    times(k, solver_id, j) = time;
    %             end
            end
        end
        save(strcat('./data/sdp/','runtime_',solver_name{solver_id},'_res','.mat'),'times','valid_id');
    else
        for solver_id = 1:size(solver_name,2)
            dat(solver_id) = load(strcat('./data/sdp/','runtime_',solver_name{solver_id},'_res','.mat'));
            times(:, solver_id, :) = dat(solver_id).times(:,mapping_id(solver_id),:);
            valid_id(:, solver_id, :) = dat(solver_id).valid_id(:,solver_id,:);
        end
        
%         tmp = times(:,6,:);
%         times(:,3:6,:) = times(:,2:5,:);
%         times(:,2,:) = tmp;
%         tmp1 = solver_name{6};
%         for i = 6:-1:3
%             solver_name{i} = solver_name{i-1};
%         end
%         solver_name{2} = tmp1;
        
        %% compute mean time
        meantimes = zeros(numel(numSam), size(solver_name,2));
        for k = 1:numel(numSam)
            for solver_id = 1:size(solver_name,2)
                valids = valid_id(k,solver_id,:)==1;
                num = numel(find(valids));
                meantimes(k, solver_id) = sum(times(k,solver_id,valids))/num;
            end
        end
        figure()
        cmap = lines(size(solver_name,2));
        cmapalpha = [cmap 0.8*ones(size(cmap,1),1)];
        xlabels = numSam;
        plot_case = solver_name;
        for solver_id = 1:size(solver_name,2) 
            plot(xlabels', meantimes(:,solver_id),'-o', 'Color', cmapalpha(solver_id,:),'LineWidth', 2.5);hold on;
        end
        fontsize = 12;
        xlabel('$Number\ of\ samples$','Interpreter','latex','FontSize',fontsize);
        ylabel('$Time: (s)$','Interpreter','latex','FontSize',fontsize);
        legend(plot_case,'Interpreter','latex','FontSize',fontsize,'Location', 'northwest');
        title('Runtime\ Comparison','Interpreter','latex','FontSize',fontsize);
        grid on;
    end
end
