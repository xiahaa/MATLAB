% generate a series of motion data and run a series of hand eye calibration
% methods on the generated data to compare their performance.
clc
clear all
close all

addpath('./solver/');
addpath('./data_gen/');
addpath('./evaluation/');
addpath('./stochastics/utils/');
addpath('../MatrixLieGroup/barfoot_tro14');
addpath('../quaternion');
addpath('./solver/atadq');
addpath('./beatiful_plot');
%Editable Variables
%------------------------------------------------------
% naaray = [10 20 30 40 50 60 70 80 90 100];
naaray = [10 20 30 40];

gmean = [0;0;0;0;0;0];	%Gaussian Noise Mean

nstd = [0 0.01 0.05 0.1 0.2 0.5 0.8 1];  %Gaussian Noise standard deviation Range
nstd1 = [0 0.01 0.05 0.1 0.2 0.5 0.8 1];  %Gaussian Noise standard deviation Range


nstd2 = [0.0 0.2 0.3 0.4 0.5 1];

sigma_r = [0 0.05 0.1 0.2 0.5 1];
sigma_t = [0 0.05 0.1 0.2 0.5 1];

usedstd = nstd2;

ElipseParam = [10, 10, 10];

%%
index = 0; % index that records the number of pairs
N = 50; % Times of simulation

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
convSolver = {@sol_tsai_lenz, ...                       %% TSAI
    @sol_park_martin, ...                               %% LIE
    @sol_horaud, ...                                    %% QSEP
    @sol_andreff, ...                                   %% KR
    @sol_dual_quaternion, ...                           %% DQ
    @sol_chou, ...                                      %% IDQ
    };

convSolver1 = {
    @HandEye_DQ, ...                                    %% DQ
    @sol_dual_quaternion,
    };

advSolver = {
             @sol_adjoint_transformation_algo, ...      %% ATA
             @sol_horaud_nlopt, ...                     %% NLOPT
             @sol_cvx1, ...                             %% SOCP
             @sol_dphec, ...                            %% GLOBAL_POLY
             @sol_dual_sdp_cvx, ...                     %% DUAL_SDP
             @sol_cvx2, ...                             %% SCF
             @sol_manifold_opt_SE3};                    %% SE3OPT
         
%  advSolver1 = {
%              @sol_horaud_nlopt, ...                     %% NLOPT
%              @sol_cvx1, ...                             %% SOCP
%              @sol_dphec, ...                            %% GLOBAL_POLY
%              @sol_dual_sdp_cvx, ...                     %% DUAL_SDP
%              @sol_cvx2, ...                             %% SCF
%              @sol_manifold_opt_SE3};                    %% SE3OPT

advSolver1 = {
             @sol_horaud_nlopt, ...                     %% NLOPT
             };


 usedsolver = convSolver;
         
 prefix = 'data/conv/convCmp';
         
for j = 1:numel(usedstd)
    for i =  1:numel(naaray)
        numPair = naaray(i);
        %Trajectory Generation
        %------------------------------------------------------
        A = zeros(4,4,numPair);
        B = zeros(4,4,numPair); 
        num = numPair + 2;
        rng();
        for k = 1:N
            disp([i,j,k])
            %%  rand X
            x = [randn(3,1) + 0.5;randn(3,1)]; X = expm(se3_vec(x));
            %% A,B
            trajParam = [.5, .5, .5, 0, 0];
            t2 = (0:(2*pi)/((num)):2*pi);
            twist = 0.0*sin(16*t2);
            [A1, B1] = AB_genTrajTwist(X, ElipseParam, trajParam, num/2, twist);
            trajParam = [.5, .5, .5, 0, 0.5*pi];
            [A2, B2] = AB_genTrajTwist(X, ElipseParam, trajParam, num/2, twist);
            A = cat(3, A1, A2);
            B = cat(3, B1, B2);

%             sigma_r = usedstd(j);
%             sigma_t = usedstd(j);
%             
%             for kx = 1:numPair       
% %                 Tnoise = vec2tran([randn(3,1).*sigma_t;randn(3,1).*sigma_r]);
% %                 Xn = X * Tnoise;
% 
%                 B(:,:,kx) = vec2tran(randn(6,1));
%                 A(:,:,kx) = X*B(:,:,kx)*X^-1;
%                 
% %                 RBn = get_rotation_noise(sigma_r);
% %                 tBn = randn(3,1)*sigma_t;
% %                 B(1:3,1:3,kx) = RBn*B(1:3,1:3,kx);
% %                 B(1:3,4,kx) = B(1:3,4,kx) + tBn;
% %                 
% %                 RAn = get_rotation_noise(sigma_r);
% %                 tAn = randn(3,1)*sigma_t;
% %                 A(1:3,1:3,kx) = RAn*A(1:3,1:3,kx);
% %                 A(1:3,4,kx) = A(1:3,4,kx) + tAn;
% 
%                 %% new noise addition
% %                 Tnoise = vec2tran([randn(3,1).*sigma_t;randn(3,1).*sigma_r]);
% %                 A(:,:,kx) =A(:,:,kx);% Tnoise
% %                 Tnoise = vec2tran([randn(3,1).*sigma_t;randn(3,1).*sigma_r]);
% %                 B(:,:,kx) = B(:,:,kx);%Tnoise * 
%             end
            
            

            %% perturbate B
            B = sensorNoise(B, gmean, usedstd(j), 1);
            
%             figure
%             for ii = 1:size(A,3)
%                 ta = A(1:3,4,ii);tb = B(1:3,4,ii);
%                 Ra = A(1:3,1:3,ii);Rb = B(1:3,1:3,ii);
%                 trplot( A(:,:,ii),'color','r');hold on;
%                 trplot( B(:,:,ii),'color','b');
%             end
%             axis equal;
%             view(3)
            figure
            axis equal;grid on;
            view(3);
            %% re-arrange inputs
            for ii = 1:size(A,3)
                T1(ii,:,:) = B(:,:,ii);
                T2(ii,:,:) = A(:,:,ii);
                T = A(:,:,ii); 
                DrawAxis(T, 0.5, 'r', 'g', 'b');hold on;
                T = B(:,:,ii); 
                DrawAxis(T, 0.5, 'r', 'g', 'b');hold on;
                p1 = A(1:3,4,ii);
                p2 = B(1:3,4,ii);
                p = [p1 p2];
                disp(["baseline: ",num2str(norm(p1-p2))]);
                plot3(p(1,:),p(2,:),p(3,:),'k-','LineWidth',2);
            end
            %% Solution for X
            %-------------------------------------------------------
            numsam = size(T1,1);
            for kk = 1:size(usedsolver, 2)
                handle_sol = usedsolver{kk};
                tic
                TX = handle_sol(T1,T2,numsam);
                tsol(kk) = toc;
                if isempty(TX) == false
                    xsol(1:4,1:4,kk) = TX;
                    flag(kk) = 1;
                else
                    flag(kk) = 0;
                end
            end
            xsols{k} = xsol;
            tsols{k} = tsol;
            Xs{k} = X;
            flags{k} = flag;
            As{k} = A;
            Bs{k} = B;
        end
        noisylv = num2str(usedstd(j));
        noisylv = replace(noisylv,'.','_');
        filename = strcat(prefix,'_',num2str(numPair), '_', noisylv, '.mat');
        save(filename,'xsols','tsols','Xs','flags','As','Bs');
    end
end
disp('finishing....');
