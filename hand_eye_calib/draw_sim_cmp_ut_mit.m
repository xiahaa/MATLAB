function draw_sim_cmp_ut_mit
% this functionality will draw results on simulated data by using the
% generator given in the MIT folder.
    clc;clear all;close all;
    addpath('../3rdparty/mit3dslam');
    addpath('../beautiful_plot');
    addpath('./solver/');
    addpath('../MatrixLieGroup/barfoot_tro14');
    addpath('../quaternion');
    addpath('./solver/atadq');
    addpath('../beautiful_plot/aboxplot');

    id = 4;
    name = {'circle100','line100','rotation100','shape8100','smallr100'};
%     stds = [0 0.1 0.25 0.5 0.75 1];    
%     stds = [0 0.1 0.25 0.5 0.75 1];   
%     stds = [0 0.1 0.25 0.5];   
    stds = [0.05 0.15 0.25 0.35 0.45 0.55];
    xlabels = categorical(convertStringsToChars(string(stds)));
    box_labels = convertStringsToChars(string(stds));

%     truth = [0.9511    0.1409    0.1761    0.2113   -0.1109    0.2818    0.1074    0.2219];
%     dq = truth(1:4)';
%     dqd = truth(5:8)';
%     R = q2R(dq);
%     t = 2.*qprod(dqd, conjugateq(dq));
%     t = t(2:4);
%     Tt = [R t;[0 0 0 1]];
    
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
%     %     @HandEye_DQ, ...                                    %% DQ
%     %     @sol_improved_dual_quaternion, ...                   %% IDQ
         
    solver_name = {'TS','LIE','QS','KR','DQ','CHOU','ATA','GPOLY','DUAL','SCF','SDR','QNL','SOCP','SE3'};
    solver_name1 = {'QS','KR','DQ','ATA','GPOLY','QNL','SOCP','SE3'};
    solver_name2 = {'KR','SOCP','ATA','GPOLY','DUAL','SDR'};
    solver_name2 = {'KR','NLQ','SOCP','ATA','GPOLY','SE3'};
    plot_case = solver_name2;
    
    numstd = size(stds,2);
    numsolver = size(plot_case,2);
    
    res_error_r = cell(numstd,numsolver);%, 100);
    res_error_t = cell(numstd,numsolver);%, 100);
    res_error_T = cell(numstd,numsolver);%, 100);
    
%     res_error_r = zeros(numstd*numsolver*100,1);%, 100);
%     res_error_t = zeros(numstd*numsolver*100,1);%, 100);
%     res_error_T = zeros(numstd*numsolver*100,1);%, 100);    
    
    holder = ones(numstd*numsolver*100,1);
    
    runtimes = zeros(numsolver, numstd);
%     lblcl1 = nominal(holder,{'1'});
%     lblcl2 = nominal(holder,{'2'});
    for solver_id = 1:numsolver
        clear res;
        res = load(strcat('C:/Users/xiahaa/Documents/MATLAB/hand_eye_calib/result/ut_mit_sim/',name{id},'_',plot_case{solver_id},'.mat'),'res');
        meantime = mean(res.res.times,2)';
        runtimes(solver_id,:) = meantime(1:numstd);
        for j = 1:numstd
%             ids = (solver_id-1)*(numstd*100) + (j-1)*100;
%             res_error_r(ids+1:ids+100) = res.res.error_r(j,:);
%             res_error_t(ids+1:ids+100) = res.res.error_r(j,:);
%             res_error_T(ids+1:ids+100) = res.res.error_r(j,:);
            res_error_r{j,  solver_id} = res.res.error_r(j,:);
            res_error_t{j,  solver_id} = res.res.error_t(j,:);
            res_error_T{j ,solver_id} = res.res.error_T(j,:);
%             lblcl1(ids+1:ids+100) = solver_name{solver_id};
%             lblcl2(ids+1:ids+100) = num2str(stds(j));
        end
    end
    
    fig = figure();
    dtugrey = [0.94 0.94 0.94];
%     set(gca,'Color',[dtugrey 0.2]);
    cmap = lines(numsolver);
    cmapalpha = [cmap 0.8*ones(size(cmap,1),1)];
    for solver_id = 1:numsolver 
        plot(xlabels, runtimes(solver_id,:),'-o', 'Color', cmapalpha(solver_id,:),'LineWidth', 2.5);hold on;
    end
    fontsize = 12;
    xlabel('$Standard\ deviation\ of\ added\ noise$','Interpreter','latex','FontSize',fontsize);
    ylabel('$Time: (s)$','Interpreter','latex','FontSize',fontsize);
    legend(plot_case,'Interpreter','latex','FontSize',fontsize,'Location', 'northwest');
    title('Runtime\ Comparison','Interpreter','latex','FontSize',fontsize);
    grid on;
%     set(gca,'Color',[dtugrey 0.2]);
    
    %% R
    cmapalpha = [cmap 0.3*ones(size(cmap,1),1)];
    fig = figure();
%      lbl = [lblcl1,lblcl2];
%      hierarchicalBoxplotSimple(res_error_r,lbl);
    multiple_boxplot(res_error_r, box_labels, plot_case, cmapalpha');
    fontsize = 12;
    xlabel('$Standard\ deviation\ of\ added\ noise$','Interpreter','latex','FontSize',fontsize);
    ylabel('$E_{\mathbf{R}}$','Interpreter','latex','FontSize',fontsize);
%     legend(solver_name,'Interpreter','latex','FontSize',8,'Location', 'northeast');
    title('Rotation Error\ Comparison','Interpreter','latex','FontSize',fontsize);
    ylim([0,2]);
    %% t
    figure();
    multiple_boxplot(res_error_t, box_labels, plot_case, cmapalpha');
    fontsize = 12;
    xlabel('$Standard\ deviation\ of\ added\ noise$','Interpreter','latex','FontSize',fontsize);
    ylabel('$E_{\mathbf{t}}$','Interpreter','latex','FontSize',fontsize);
%     legend(solver_name,'Interpreter','latex','FontSize',8,'Location', 'northeast');
    title('Translation Error\ Comparison','Interpreter','latex','FontSize',fontsize);
    ylim([0,3]);
    %% T
    figure();
    multiple_boxplot(res_error_T, box_labels, plot_case, cmapalpha');
    fontsize = 12;
    xlabel('$Standard\ deviation\ of\ added\ noise$','Interpreter','latex','FontSize',fontsize);
    ylabel('$E_{\mathbf{T}}$','Interpreter','latex','FontSize',fontsize);
%     legend(solver_name,'Interpreter','latex','FontSize',8,'Location', 'northeast');
    title('Total Error\ Comparison','Interpreter','latex','FontSize',fontsize);
    ylim([0,3]);
    
%     fig = figure();
%     h = boxplot2(res_error_T, x);
%     for ii = 1:numsolver
%         structfun(@(x) set(x(ii,:), 'color', cmap(ii,:), ...
%             'markeredgecolor', cmap(ii,:)), h);
%     end
%     set([h.lwhis h.uwhis], 'linestyle', '-');
%     set(h.out, 'marker', '.');
%     fontsize = 12;
%     xlabel('$std\ of\ added\ noise$','Interpreter','latex','FontSize',fontsize);
%     ylabel('$E_{\mathbf{T}}$','Interpreter','latex','FontSize',fontsize);
%     legend(solver_name,'Interpreter','latex','FontSize',8,'Location', 'northeast');
%     title('Total Error\ Comparison','Interpreter','latex','FontSize',fontsize);
%     ylim([0,10]);
end