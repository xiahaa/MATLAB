clear all; clc;
close all;

addpath ../utils/
addpath(genpath('../3rdparty/PnP_Toolbox-master/PnP_Toolbox-master/code'));

% experimental parameters
nl = 5; % level of noise
npts  = 50; % number of points
pouts = 0.6:0.1:0.80;%:0.1:0.1;%[0.0:0.05:0.80];  % percentage of outliers
num   = 100; % total number of trials

% compared methods
A= zeros(size(npts));
B= zeros(num,1);

name= {'Test','REPPnP', 'R1PPnP'};%'RNSC P3P','RNSC RP4P RPnP','RNSC P3P OPnP','RNSC P3P ASPnP', 'REPPnP', 'R1PPnP'
f = {@epnp_re_ransac2,@REPPnP, @R1PPnP};% kP3P,@RPnP,@kP3P,@kP3P
f2 = {[],[],[]}; %post ransac method: [],@RPnP,@OPnP,@ASPnP,
ransacsamples = {0,0,0};%3,4,3,3,0,0}; %number of samples to apply ransac
marker= {'o','d','>'};%,'d','>','<','^','x'};
color= {'g','r','k'};%,[1,0.5,0],'c','b','r','k'};
markerfacecolor=  {'g','r','k'};%,[1,0.5,0],'c','b','r','k'};

method_list= struct('name', name, 'f', f,'f2', f2, 'ransac',ransacsamples, 'mean_r', A, 'mean_t', A,...
    'med_r', A, 'med_t', A, 'std_r', A, 'std_t', A, 'r', B, 't', B,...
    'marker', marker, 'color', color, 'markerfacecolor', markerfacecolor);

% experiments
for i= 1:length(pouts)
        
        npt  = npts;
        pout = pouts(i);
        
        fprintf('npt = %d (sg = %d px) (%2.0f%%): ', npt, nl, pout*100);
    
        for k= 1:length(method_list)
            method_list(k).c = zeros(1,num);
            method_list(k).e = zeros(1,num);
            method_list(k).r = zeros(1,num);
            method_list(k).t = zeros(1,num);
        end

        index_fail = cell(1,length(name));
        
        for j= 1:num
            % camera's parameters
            f  = 1000;
            
            % generate 3d coordinates in camera space
            Xc  = [xrand(1,npt,[-2 2]); xrand(1,npt,[-2 2]); xrand(1,npt,[4 8])];
            t   = mean(Xc,2);
            R   = rodrigues(randn(3,1));
            XXw = inv(R)*(Xc-repmat(t,1,npt));
            
            % projection
            xx  = [Xc(1,:)./Xc(3,:); Xc(2,:)./Xc(3,:)]*f;
           % randomvals = randn(2,npt);
            xxn = xx + randn(2,npt) * nl;
                         
            %generate outliers (assigning to some 3D points more than one 2D correspondence)
            if (pout ~= 0)
                nout = max(1,round((npt * pout)/(1-pout))); %at least we want 1 outlier

                idx  = randi(npt,1,nout);
                XXwo = XXw(:,idx);
            else
               nout = 0;
               XXwo = []; 
            end
            % assignation of random 2D correspondences
            xxo  = [xrand(1,nout,[min(xxn(1,:)) max(xxn(1,:))]); xrand(1,nout,[min(xxn(2,:)) max(xxn(2,:))])];
            
            % test
            % pose estimation
            R1 = []; t1 = []; inliers = [];
            for k= 1:length(method_list)
                 if strcmp(method_list(k).name, 'R1PPnP')
                    tic;
                    minInlierErrorinPixel = 10.0;
                    [R1,t1]= method_list(k).f([XXw, XXwo],[xxn, xxo]/f, f, minInlierErrorinPixel);
                    tcost = toc;
                 elseif strcmp(method_list(k).name, 'REPPnP')
                    tic;
                    [R1,t1,~,robustiters,terr]= method_list(k).f([XXw, XXwo],[xxn, xxo]/f);
                    tcost = toc;
                 else
                     tic
                     [R1,t1]= method_list(k).f([XXw, XXwo],[xxn, xxo]/f);
                     tcost = toc;
                 end

                %no solution
                if size(t1,2) < 1
                   % disp(['The solver - ',method_list(k).name,' - returns no solution!!!\n']);
                    index_fail{k} = [index_fail{k}, j];
                    continue;
                   % break;
                elseif (sum(sum(sum(imag(R1).^2))>0) == size(R1,3) || sum(sum(imag(t1(:,:,1)).^2)>0) == size(t1,2))
                    index_fail{k} = [index_fail{k}, j];
                    continue;
                end
                %choose the solution with smallest error 
                error = inf;
                for jjj = 1:size(R1,3)
                    if (sum(sum(imag(R1(:,:,jjj)).^2)) + sum(imag(t1(:,jjj)).^2) > 0)
                        break
                    end            
                    
                    tempy = cal_pose_err([R1(:,:,jjj) t1(:,jjj)],[R t]);
                    if sum(tempy) < error
                        cost  = tcost;
                        %L2 error is computed without taing into account the outliers
                        ercorr= mean(sqrt(sum((R1(:,:,jjj) * XXw +  t1(:,jjj) * ones(1,npt) - Xc).^2)));
                        y     = tempy;
                        error = sum(tempy);
                    end
                end

                method_list(k).robustiters(j) = 1;%robustiters;
                method_list(k).c(j)= cost * 1000;
                method_list(k).e(j)= ercorr;
                method_list(k).r(j)= y(1);
                method_list(k).t(j)= y(2);
            end

%             showpercent(j,num);
        end
        fprintf('\n');
    
        % save result
        for k= 1:length(method_list)
            
            %results without deleting solutions
            tmethod_list = method_list(k);
            method_list(k).robustiters(index_fail{k}) = [];
            method_list(k).c(index_fail{k}) = [];
            method_list(k).e(index_fail{k}) = [];
            method_list(k).r(index_fail{k}) = [];
            method_list(k).t(index_fail{k}) = [];

            method_list(k).pfail(i) = 100 * numel(index_fail{k})/num;
            
            method_list(k).mean_rit(i) = mean(method_list(k).robustiters);
            method_list(k).mean_c(i)= mean(method_list(k).c);
            method_list(k).mean_e(i)= mean(method_list(k).e);
            method_list(k).med_rit(i) = median(method_list(k).robustiters);
            method_list(k).med_c(i)= median(method_list(k).c);
            method_list(k).med_e(i)= median(method_list(k).e);
            method_list(k).std_rit(i) = std(method_list(k).robustiters);
            method_list(k).std_c(i)= std(method_list(k).c);
            method_list(k).std_e(i)= std(method_list(k).e);

            method_list(k).mean_r(i)= mean(method_list(k).r);
            method_list(k).mean_t(i)= mean(method_list(k).t);
            method_list(k).med_r(i)= median(method_list(k).r);
            method_list(k).med_t(i)= median(method_list(k).t);
            method_list(k).std_r(i)= std(method_list(k).r);
            method_list(k).std_t(i)= std(method_list(k).t);
           
            %results deleting solutions where not all the methods finds one
            tmethod_list.robustiters(unique([index_fail{:}])) = [];
            tmethod_list.c(unique([index_fail{:}])) = [];
            tmethod_list.e(unique([index_fail{:}])) = [];
            tmethod_list.r(unique([index_fail{:}])) = [];
            tmethod_list.t(unique([index_fail{:}])) = [];
            
            method_list(k).deleted_mean_rit(i) = mean(tmethod_list.robustiters);
            method_list(k).deleted_mean_c(i)= mean(tmethod_list.c);
            method_list(k).deleted_mean_e(i)= mean(tmethod_list.e);
            method_list(k).deleted_med_rit(i) = median(tmethod_list.robustiters);
            method_list(k).deleted_med_c(i)= median(tmethod_list.c);
            method_list(k).deleted_med_e(i)= median(tmethod_list.e);
            method_list(k).deleted_std_rit(i) = std(tmethod_list.robustiters);
            method_list(k).deleted_std_c(i)= std(tmethod_list.c);
            method_list(k).deleted_std_e(i)= std(tmethod_list.e);

            method_list(k).deleted_mean_r(i)= mean(tmethod_list.r);
            method_list(k).deleted_mean_t(i)= mean(tmethod_list.t);
            method_list(k).deleted_med_r(i)= median(tmethod_list.r);
            method_list(k).deleted_med_t(i)= median(tmethod_list.t);
            method_list(k).deleted_std_r(i)= std(tmethod_list.r);
            method_list(k).deleted_std_t(i)= std(tmethod_list.t);
        end
        
        prealout(i) = nout/npt; %exact percent of outliers    
        
%        save ordinary3DresultsOutliersRansac_lastiter method_list npt pouts prealout;
end

% addpath ./../3rdparty/PnP_Toolbox-master/PnP_Toolbox-master/code/plot_funcs/;
% plotOrdinary3DoutliersRansac;


close all;

i= 0; w= 300; h= 300;

yrange= [0 max([method_list(:).deleted_mean_r])];
yrange= [0 3];

figure('color','w','position',[w*i,100,w,h]);%i=i+1;
xdrawgraph(pouts*100,yrange,method_list,'deleted_mean_r','Mean Rotation Error',...
    '% of outliers','Rotation Error (degrees)');

% figure('color','w','position',[w*i,100+h,w,h]);i=i+1;
% xdrawgraph(pouts*100,yrange,method_list,'deleted_med_r','Median Rotation Error',...
%     '% of outliers','Rotation Error (degrees)');

yrange= [0 max([method_list(:).deleted_mean_t])];
yrange= [0 3];

figure('color','w','position',[w*i,100,w,h]);%i=i+1;
xdrawgraph(pouts*100,yrange,method_list,'deleted_mean_t','Mean Translation Error',...
    '% of outliers','Translation Error (%)');

figure('color','w','position',[w*i,100+h,w,h]);i=i+1;
xdrawgraph(pouts*100,yrange,method_list,'deleted_med_t','Median Translation Error',...
    '% of outliers','Translation Error (%)');

yrange= [0 max([method_list(:).deleted_mean_e])];

figure('color','w','position',[w*i,100,w,h]);%i=i+1;
xdrawgraph(pouts*100,yrange,method_list,'deleted_mean_e','Mean L2 Error',...
    '% of outliers','L2 error');

figure('color','w','position',[w*i,100+h,w,h]);
xdrawgraph(pouts*100,yrange,method_list,'deleted_med_e','Median L2 Error',...
    '% of outliers','L2 error');

yrange= [0 min(max(1,2*max([method_list(:).pfail])),100)];
figure('color','w','position',[w*i,100+2*h,w,h]);%i=i+1;
xdrawgraph(pouts*100,yrange,method_list,'pfail','No solution x method',...
    '% of outliers','% method fails');
i=i+1;

% yrange= [0 2*max([method_list(:).med_e])];
% 
% figure('color','w','position',[w*i,100,w,h]);%i=i+1;
% xdrawgraph(pouts*100,yrange,method_list,'deleted_mean_e','Mean L2 Error',...
%     '% of outliers/inliers','L2 error');
% 
% figure('color','w','position',[w*i,100+h,w,h]);i=i+1;
% xdrawgraph(pouts*100,yrange,method_list,'deleted_med_e','Median L2 Error',...
%     '% of outliers/inliers','L2 error');
%  

yrange= [0 max([method_list(:).deleted_mean_c])];

figure('color','w','position',[w*i,100,w,h]);%i=i+1;
xdrawgraph(pouts*100,yrange,method_list,'deleted_mean_c','Mean Cost',...
    '% of outliers','Cost (ms)');

figure('color','w','position',[w*i,100+h,w,h]);i=i+1;
xdrawgraph(pouts*100,yrange,method_list,'deleted_med_c','Median Cost',...
    '% of outliers','Cost (ms)');

yrange= [0 max([method_list(:).deleted_mean_rit])];

figure('color','w','position',[w*i,100,w,h]);%i=i+1;
xdrawgraph(pouts*100,yrange,method_list,'deleted_mean_rit','Mean Robust Iters',...
    '% of outliers','Iterations');

figure('color','w','position',[w*i,100+h,w,h]);i=i+1;
xdrawgraph(pouts*100,yrange,method_list,'deleted_med_rit','Median Robust Iters',...
    '% of outliers','Iterations');
