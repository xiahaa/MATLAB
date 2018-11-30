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

%Editable Variables
%------------------------------------------------------
naaray = [10 20 30 40 50 60 70 80 90 100];

gmean = [0;0;0;0;0;0];	%Gaussian Noise Mean

nstd = [0 0.01 0.05 0.1 0.2 0.5 0.8 1];  %Gaussian Noise standard deviation Range
nstd1 = [0 0.01 0.05 0.1 0.2 0.5 0.8 1];  %Gaussian Noise standard deviation Range

nstd2 = [0 0.05 0.1 0.2 0.5 1];

sigma_r = [0 0.05 0.1 0.2 0.5 1];
sigma_t = [0 0.05 0.1 0.2 0.5 1];

usedstd = nstd2;

ElipseParam = [10, 10, 10];

%%
index = 0; % index that records the number of pairs
N = 100; % Times of simulation

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
convSolver = {@sol_tsai_lenz, ...                       %% TSAI
    @sol_park_martin, ...                               %% LIE
    @sol_horaud, ...                                    %% QSEP
    @sol_andreff, ...                                   %% KR
    @HandEye_DQ, ...                                    %% DQ
    @sol_improved_dual_quaternion, ...                  %% IDQ
    };

advSolver = {
             @sol_adjoint_transformation_algo, ...      %% ATA
             @sol_horaud_nlopt, ...                     %% NLOPT
             @sol_cvx1, ...                             %% SOCP
             @sol_dphec, ...                            %% GLOBAL_POLY
             @sol_dual_sdp_cvx, ...                     %% DUAL_SDP
             @sol_cvx2, ...                             %% SCF
             @sol_manifold_opt_SE3};                    %% SE3OPT
         
 advSolver1 = {
             @sol_cvx2, ...                             %% SCF
             @sol_manifold_opt_SE3};                    %% SE3OPT

 usedsolver = advSolver1;
         
 prefix = 'data/SCF/advCmp';
         
for j = 1:numel(usedstd)
    for i =  1:numel(naaray)
        numPair = naaray(i);
        %Trajectory Generation
        %------------------------------------------------------
        A = zeros(4,4,numPair);
        B = zeros(4,4,numPair); 
        num = numPair;
        
        for k = 1:N
            %%  rand X
            x = randn(6,1); X = expm(se3_vec(x));
            %% A,B
%             trajParam = [.5, .5, .5, 0, 0];
%             [A1, B1] = AB_genTraj(X, ElipseParam, trajParam, num/2);
%             trajParam = [.5, .5, .5, 0, 0.5*pi];
%             [A2, B2] = AB_genTraj(X, ElipseParam, trajParam, num/2);
%             A = cat(3, A1, A2);
%             B = cat(3, B1, B2);
            sigma_r = usedstd(j);
            sigma_t = usedstd(j);
            for kx = 1:numPair       
                B(:,:,kx) = vec2tran(randn(6,1));
                A(:,:,kx) = X*B(:,:,kx)*X^-1;
                
                
                RBn = get_rotation_noise(sigma_r);
                tBn = randn(3,1)*sigma_t;
                B(1:3,1:3,kx) = RBn*B(1:3,1:3,kx);
                B(1:3,4,kx) = B(1:3,4,kx) + tBn;
                
                RAn = get_rotation_noise(sigma_r);
                tAn = randn(3,1)*sigma_t;
                A(1:3,1:3,kx) = RAn*A(1:3,1:3,kx);
                A(1:3,4,kx) = A(1:3,4,kx) + tAn;
            end

            %% perturbate B
%             B = sensorNoise(B, gmean, usedstd(j), 1);
            

            %% re-arrange inputs
            for ii = 1:size(A,3)
                T1(ii,:,:) = B(:,:,ii);
                T2(ii,:,:) = A(:,:,ii);
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
        end
        noisylv = num2str(usedstd(j));
        noisylv = replace(noisylv,'.','_');
        filename = strcat(prefix,'_',num2str(numPair), '_', noisylv, '.mat');
        save(filename,'xsols','tsols','Xs','flags');
    end
end
disp('finishing....');
