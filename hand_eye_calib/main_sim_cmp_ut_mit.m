function main_sim_cmp_ut_mit
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
    
    convSolver = {
        @sol_andreff, ...                                   %% KR
        @sol_adjoint_transformation_algo, ...               %% ATA
        @sol_dphec, ...                                     %% GLOBAL_POLY
        @sol_dual_sdp_cvx, ...                              %% DUAL_SDP
        @sol_cvx_sdp, ...                                   %% SDP
        @sol_cvx1, ...                                      %% SOCP
        };

    solver_name = {'KR','ATA','GPOLY','DUAL','SDR','SOCP'};
    
    usedsolver = convSolver;
    
    for solver_id = 1:size(solver_name,2)
        times = zeros(numel(stds),100);
        error_r = zeros(numel(stds),100);
        error_t = zeros(numel(stds),100);
        error_T = zeros(numel(stds),100);
        valid_id = ones(numel(stds),100);
        for k = 1:numel(stds)
            clear run;
            std = stds(k);
            filename = strcat(prefix,name{id},'_',num2str(std),suffix);
            load(filename);
            N = min(100, size(run.observations,2));%100
            calibrationEstimates = zeros(N,8);

            handle_sol = usedsolver{solver_id};
            
            for j = 1:N
                disp([solver_id k,j])
                sensor1_expressedIn_prevSensor1 = run.observations{j}.sensor1_expressedIn_prevSensor1;
                sensor2_expressedIn_prevSensor2 = run.observations{j}.sensor2_expressedIn_prevSensor2;
                hand = EulerStateToDualQuat(sensor1_expressedIn_prevSensor1);
                eye = EulerStateToDualQuat(sensor2_expressedIn_prevSensor2);

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
                tic
                TX = handle_sol(T1,T2,num);
                time = toc;
                
                if isempty(TX)
                    valid_id(k,j) = 0;
                    continue;
                end
                times(k,j) = time;
                [err1,err2,err3]  = errors2(TX, Tt);
                error_r(k,j) = err1;
                error_t(k,j) = err2;
                error_T(k,j) = err3;
            end
    %         res{k} = struct('error_r', error_r,'error_t', error_t,'error_T', error_T, 'time', times);
        end
        fig = figure();
        box_labels = categorical({'0','0.1','0.25','0.5','0.75','1'});
        box_labels = reordercats(box_labels,{'0','0.1','0.25','0.5','0.75','1'});
        boxplot(error_r', box_labels);
        res.error_r = error_r;
        res.error_t = error_t;
        res.error_T = error_T;
        res.times = times;
        save(strcat('C:/Users/xiahaa/Documents/MATLAB/hand_eye_calib/result/ut_mit_sim/',name{id},'_',solver_name{solver_id},'.mat'),'res');
    end
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