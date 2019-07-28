function draw_test_1
% this functionality will draw results on simulated data by using the
% generator given in the MIT folder.
    clc;clear all;close all;
    basepath = 'C:/Users/xiahaa/Documents/MATLAB/hand_eye_calib/';
    addpath(strcat(basepath,'../beautiful_plot'));
    addpath(strcat(basepath,'./solver/'));
    addpath(strcat(basepath,'../MatrixLieGroup/barfoot_tro14'));
    addpath(strcat(basepath,'../quaternion'));
    addpath(strcat(basepath,'./solver/atadq'));
    addpath(strcat(basepath,'../beautiful_plot/aboxplot'));

    %% Initialize Parameters
    
    % test 1
%     stds = 0.3:0.3:1.5; % Gaussian Noise standard deviation Range
%     base_name = 'test_1';
%     xlbs = 'Standard deviation';
    
    % test 2
%     stds = 0.3:0.3:1.5; % Gaussian Noise standard deviation Range
%     base_name = 'test_2';
%     xlbs = 'Standard deviation';
    
    % test 3
    stds = 0.02:0.02:0.10; % Gaussian Noise standard deviation Range
    base_name = 'test_3';
    xlbs = 'Standard deviation of noise';
    
    % test 2
%     percentage_of_added_sample = 0:0.1:0.5;
%     base_name = 'test_2';
%     xlbs = 'Percentage of additional samples';

    % test 3
%     outlier_percentage = 0:0.1:0.5;
%     base_name = 'test_3';
%     xlbs = 'Percentage of additional outliers';
    
    test_case = stds;
    
    solver_name = {'B1','B2','BS'};

    xlabels = categorical(convertStringsToChars(string(test_case)));
    box_labels = convertStringsToChars(string(test_case));

    plot_case = solver_name;
    
    numcases = size(test_case,2);
    numsolver = size(plot_case,2);
    
    res_error_r = cell(numcases,numsolver);%, 100);
    res_error_t = cell(numcases,numsolver);%, 100);
        
    runtimes = zeros(numsolver, numcases);
    cmap = lines(numsolver);
    cmapalpha = [cmap 0.8*ones(size(cmap,1),1)];
    
    showtypes = {'-r','-g','-b'};
    
    m_error_r = zeros(numcases,numsolver);%, 100);
    m_error_t = zeros(numcases,numsolver);%, 100);

    for solver_id = 1:numsolver
        clear res;
        res = load(strcat('C:/Users/xiahaa/Documents/MATLAB/hand_eye_calib/result/sto/',base_name,plot_case{solver_id},'.mat'),'res');
        meantime = mean(res.res.times,2)';
        runtimes(solver_id,:) = meantime(1:numcases);
        
%         x = 1:numstd;
%         y = zeros(50,numstd);
        for j = 1:numcases
            res_error_r{j,  solver_id} = res.res.error_r(j,:);
            res_error_t{j,  solver_id} = res.res.error_t(j,:);
            
%             y(:,j) = res.res.error_t(j,:)';
            m_error_r(j,solver_id) = mean(res.res.error_r(j,:));
            m_error_t(j,solver_id) = mean(res.res.error_t(j,:));
        end
%         shadedErrorBar(x, y, {@mean,@std}, 'lineprops', showtypes(solver_id));hold on;
    end
%     fig = figure();
%     cmap = lines(numsolver);
%     cmapalpha = [cmap 0.8*ones(size(cmap,1),1)];
%     for solver_id = 1:numsolver 
%         plot(xlabels, runtimes(solver_id,:),'-o', 'Color', cmapalpha(solver_id,:),'LineWidth', 2.5);hold on;
%     end
%     fontsize = 12;
%     xlabel('$Standard\ deviation\ of\ added\ noise$','Interpreter','latex','FontSize',fontsize);
%     ylabel('$Time: (s)$','Interpreter','latex','FontSize',fontsize);
%     legend(plot_case,'Interpreter','latex','FontSize',fontsize,'Location', 'northwest');
%     title('Runtime\ Comparison','Interpreter','latex','FontSize',fontsize);
%     grid on;
% %     
%     %% R
%     cmapalpha = [cmap 0.3*ones(size(cmap,1),1)];
%     fig = figure();
%     multiple_boxplot(res_error_r, box_labels, plot_case, cmapalpha');
%     fontsize = 12;
%     xlabel('$Standard\ deviation\ of\ added\ noise$','Interpreter','latex','FontSize',fontsize);
%     ylabel('$E_{\mathbf{R}}$','Interpreter','latex','FontSize',fontsize);
%     title('Rotation Error\ Comparison','Interpreter','latex','FontSize',fontsize);
% %     ylim([0,2]);
%     %% t
%     figure();
%     multiple_boxplot(res_error_t, box_labels, plot_case, cmapalpha');
%     fontsize = 12;
%     xlabel('$Standard\ deviation\ of\ added\ noise$','Interpreter','latex','FontSize',fontsize);
%     ylabel('$E_{\mathbf{t}}$','Interpreter','latex','FontSize',fontsize);
%     title('Translation Error\ Comparison','Interpreter','latex','FontSize',fontsize);
%     ylim([0,3]);
    
    figure
    scatter(test_case, m_error_r(:,1), 100, 'r', 'd', 'LineWidth',2);
    hold on
    scatter(test_case, m_error_r(:,2),  100, 'b', 'o', 'LineWidth',2);
    hold on
    scatter(test_case, m_error_r(:,3), 100, 'g', 's', 'LineWidth',2);
    hold on
    xlabel(xlbs)
    ylabel('Rotation error')
    legend(solver_name,'Location', 'NorthWest')
    grid on;
    xticks(test_case);
    xticklabels(cellstr(string(test_case)));
    ylim([0, 1])
    % axis([0 100 0 4])

    figure
    scatter(test_case, m_error_t(:,1), 100, 'r', 'd', 'LineWidth',2);
    hold on
    scatter(test_case, m_error_t(:,2),  100, 'b', 'o', 'LineWidth',2);
    hold on
    scatter(test_case, m_error_t(:,3), 100, 'g', 's', 'LineWidth',2);
    hold on
    xlabel(xlbs)
    ylabel('Translation error')
    legend(solver_name,'Location', 'NorthWest')
    grid on;
    xticks(test_case);
    xticklabels(cellstr(string(test_case)));
    % axis([0 100 0 4])

end